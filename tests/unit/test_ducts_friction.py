"""Darcy-factor, rough curved-duct and correlation-warning contracts."""

import warnings

import fluids
import numpy as np
import pytest

from pyskyfire import ducts
from pyskyfire.common import friction


def _transition_weight(reynolds, critical):
    end = critical * 3_500.0 / 2_300.0
    fraction = np.clip((reynolds - critical) / (end - critical), 0.0, 1.0)
    return fraction**2 * (3.0 - 2.0 * fraction)


def _straight_reference(reynolds, diameter, roughness):
    relative_roughness = (0.0 if roughness is None else roughness) / diameter
    laminar = fluids.friction.friction_laminar(reynolds)
    weight = _transition_weight(reynolds, 2_040.0)
    if weight == 0.0:
        return laminar
    turbulent = fluids.friction.Clamond(reynolds, relative_roughness)
    return (1.0 - weight) * laminar + weight * turbulent


def _curved_reference(reynolds, diameter, radius, roughness):
    coil_diameter = 2.0 * abs(radius)
    roughness = 0.0 if roughness is None else roughness
    critical = fluids.friction.helical_transition_Re_Schmidt(
        diameter, coil_diameter,
    )
    weight = _transition_weight(reynolds, critical)
    laminar = fluids.friction.helical_laminar_fd_Schmidt(
        reynolds, diameter, coil_diameter,
    )
    if weight == 0.0:
        return laminar
    turbulent = fluids.friction.helical_turbulent_fd_Schmidt(
        reynolds, diameter, coil_diameter, roughness,
    )
    return (1.0 - weight) * laminar + weight * turbulent


@pytest.mark.parametrize(
    "reynolds", [1_000.0, 2_040.0, 2_500.0, 3_500.0, 100_000.0],
)
@pytest.mark.parametrize("diameter", [0.001, 0.02])
@pytest.mark.parametrize("roughness", [None, 0.0, 2e-5])
def test_fluids_straight_factor_blends_direct_branches(
    reynolds, diameter, roughness,
) -> None:
    actual = friction.darcy_friction_factor(
        reynolds, diameter, roughness, model="fluids",
    )
    assert actual == pytest.approx(
        _straight_reference(reynolds, diameter, roughness), rel=1e-12,
    )


def test_transition_weight_is_exported() -> None:
    assert ducts.transition_weight is friction.transition_weight
    assert "transition_weight" in friction.__all__


@pytest.mark.parametrize("critical", [500.0, 2_040.0, 8_000.0])
def test_shared_transition_weight_scales_with_critical(critical) -> None:
    end = critical * 3_500.0 / 2_300.0
    assert friction.transition_weight(0.5 * critical, critical) == 0.0
    assert friction.transition_weight(critical, critical) == 0.0
    assert friction.transition_weight(end, critical) == 1.0
    assert friction.transition_weight(2.0 * end, critical) == 1.0
    for fraction in (0.2, 0.5, 0.8):
        reynolds = critical + fraction * (end - critical)
        expected = fraction**2 * (3.0 - 2.0 * fraction)
        assert friction.transition_weight(reynolds, critical) == pytest.approx(expected)


@pytest.mark.parametrize(
    "reynolds",
    [500.0, 3_000.0, 4_000.0, 5_000.0, 10_000.0, 21_999.0, 22_000.0,
     22_001.0, 100_000.0],
)
@pytest.mark.parametrize("roughness", [None, 0.0, 2e-5])
@pytest.mark.parametrize("radius", [0.1, -0.1])
def test_fluids_curvature_ratio_and_gradient_match_schmidt(
    reynolds, roughness, radius,
) -> None:
    diameter, density, velocity = 0.01, 850.0, 3.0
    straight = _straight_reference(reynolds, diameter, roughness)
    curved = _curved_reference(reynolds, diameter, radius, roughness)

    multiplier = friction.curvature_multiplier(
        reynolds, diameter, radius, roughness=roughness, model="fluids",
    )
    gradient = friction.friction_pressure_gradient(
        reynolds, diameter, density, velocity,
        roughness=roughness, R_curv=radius, model="fluids",
    )

    assert multiplier == pytest.approx(curved / straight, rel=1e-12)
    assert straight * multiplier == pytest.approx(curved, rel=1e-12)
    assert gradient == pytest.approx(
        -curved * density * velocity**2 / (2.0 * diameter), rel=1e-12,
    )


@pytest.mark.parametrize("reynolds", [3_000.0, 4_000.0, 5_000.0])
def test_rough_curved_laminar_flow_cancels_straight_turbulent_roughness(
    reynolds,
) -> None:
    diameter, radius, roughness = 0.01, 0.1, 1e-4
    critical = fluids.friction.helical_transition_Re_Schmidt(
        Di=diameter, Dc=2.0 * radius,
    )
    # Here the coil is still laminar, but the straight reference is turbulent.
    # A smooth denominator leaves a spurious roughness dependence in dp/ds.
    assert 2_300.0 < reynolds < critical
    assert _straight_reference(reynolds, diameter, roughness) > (
        _straight_reference(reynolds, diameter, 0.0)
    )
    expected = _curved_reference(reynolds, diameter, radius, roughness)
    assert expected == _curved_reference(reynolds, diameter, radius, 0.0)
    multiplier = friction.curvature_multiplier(
        reynolds, diameter, radius, roughness=roughness, model="fluids",
    )

    assert multiplier * friction.darcy_friction_factor(
        reynolds, diameter, roughness, model="fluids",
    ) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("join", ["critical", "turbulent"])
@pytest.mark.parametrize("side", [-1, 0, 1])
def test_fluids_uses_both_schmidt_transition_boundaries(join, side) -> None:
    diameter, radius, roughness = 0.01, 0.1, 2e-5
    critical = fluids.friction.helical_transition_Re_Schmidt(
        Di=diameter, Dc=2.0 * radius,
    )
    boundary = critical if join == "critical" else critical * 3_500.0 / 2_300.0
    reynolds = boundary if side == 0 else np.nextafter(boundary, side * np.inf)
    actual = friction.curvature_multiplier(
        reynolds, diameter, radius, roughness=roughness, model="fluids",
    )
    assert actual == pytest.approx(
        _curved_reference(reynolds, diameter, radius, roughness)
        / _straight_reference(reynolds, diameter, roughness),
        rel=1e-12,
    )


def test_fluids_calls_branch_correlations_directly(monkeypatch) -> None:
    calls = []

    def laminar(reynolds):
        calls.append(("straight laminar", reynolds))
        return 0.02

    def clamond(reynolds, relative_roughness):
        calls.append(("straight turbulent", reynolds, relative_roughness))
        return 0.04

    def transition(diameter, coil_diameter):
        calls.append(("critical", diameter, coil_diameter))
        return 3_500.0

    def curved_laminar(reynolds, diameter, coil_diameter):
        calls.append(("curved laminar", reynolds, diameter, coil_diameter))
        return 0.03

    def curved_turbulent(reynolds, diameter, coil_diameter, roughness):
        calls.append((
            "curved turbulent", reynolds, diameter, coil_diameter, roughness,
        ))
        return 0.05

    monkeypatch.setattr(fluids.friction, "friction_laminar", laminar)
    monkeypatch.setattr(fluids.friction, "Clamond", clamond)
    monkeypatch.setattr(
        fluids.friction, "helical_transition_Re_Schmidt", transition,
    )
    monkeypatch.setattr(
        fluids.friction, "helical_laminar_fd_Schmidt", curved_laminar,
    )
    monkeypatch.setattr(
        fluids.friction, "helical_turbulent_fd_Schmidt", curved_turbulent,
    )

    assert friction.darcy_friction_factor(
        4_000.0, 0.01, 2e-5, model="fluids",
    ) == 0.04
    multiplier = friction.curvature_multiplier(
        4_000.0, 0.01, -0.1, roughness=2e-5, model="fluids",
    )
    weight = _transition_weight(4_000.0, 3_500.0)
    assert multiplier == pytest.approx(
        ((1.0 - weight) * 0.03 + weight * 0.05) / 0.04,
    )
    assert ("critical", 0.01, 0.2) in calls
    assert ("curved laminar", 4_000.0, 0.01, 0.2) in calls
    assert ("curved turbulent", 4_000.0, 0.01, 0.2, 2e-5) in calls
    assert calls.count(("straight laminar", 4_000.0)) == 0
    assert calls.count(("straight turbulent", 4_000.0, 0.002)) == 2


@pytest.mark.parametrize("geometry", ["straight", "curved"])
@pytest.mark.parametrize("join", ["critical", "turbulent"])
def test_fluids_transition_is_value_and_first_derivative_continuous(
    geometry, join,
) -> None:
    diameter, radius, roughness = 0.01, 0.1, 2e-5
    if geometry == "straight":
        critical = 2_040.0

        def factor(reynolds):
            return friction.darcy_friction_factor(
                reynolds, diameter, roughness, model="fluids",
            )
    else:
        critical = fluids.friction.helical_transition_Re_Schmidt(
            diameter, 2.0 * radius,
        )

        def factor(reynolds):
            straight = friction.darcy_friction_factor(
                reynolds, diameter, roughness, model="fluids",
            )
            multiplier = friction.curvature_multiplier(
                reynolds, diameter, radius,
                roughness=roughness, model="fluids",
            )
            return straight * multiplier

    boundary = critical if join == "critical" else critical * 3_500.0 / 2_300.0
    value = factor(boundary)
    assert factor(np.nextafter(boundary, 0.0)) == pytest.approx(value, rel=1e-12)
    assert factor(np.nextafter(boundary, np.inf)) == pytest.approx(value, rel=1e-12)

    step = boundary * 1e-5
    derivative_left = (value - factor(boundary - step)) / step
    derivative_right = (factor(boundary + step) - value) / step
    assert derivative_left == pytest.approx(derivative_right, rel=1e-3)


@pytest.mark.parametrize("model", ["native", "fluids"])
@pytest.mark.parametrize("roughness", [None, 0.0, 2e-5])
@pytest.mark.parametrize("radius", [np.inf, -np.inf])
def test_infinite_radius_is_exactly_straight(model, roughness, radius) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        multiplier = friction.curvature_multiplier(
            100_000.0, 0.01, radius, roughness=roughness, model=model,
        )
        actual = friction.friction_pressure_gradient(
            100_000.0, 0.01, 850.0, 3.0,
            roughness=roughness, R_curv=radius, model=model,
        )
    straight = friction.darcy_friction_factor(
        100_000.0, 0.01, roughness, model=model,
    )
    assert multiplier == 1.0
    assert actual == pytest.approx(-straight * 850.0 * 3.0**2 / (2.0 * 0.01))


@pytest.mark.parametrize("model", ["native", "fluids"])
def test_disabled_curvature_is_not_evaluated(monkeypatch, model) -> None:
    def unexpected(*args, **kwargs):
        pytest.fail("disabled curvature correction was evaluated")

    monkeypatch.setattr(friction, "curvature_multiplier", unexpected)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        actual = friction.friction_pressure_gradient(
            100_000.0, 0.01, 850.0, 3.0, roughness=2e-5,
            R_curv=0.01, curvature_correction=False, model=model,
        )
    straight = friction.darcy_friction_factor(
        100_000.0, 0.01, 2e-5, model=model,
    )
    assert actual == pytest.approx(-straight * 850.0 * 3.0**2 / (2.0 * 0.01))


def test_correlation_warning_is_exported_user_warning() -> None:
    assert issubclass(friction.FrictionCorrelationWarning, UserWarning)
    assert ducts.FrictionCorrelationWarning is friction.FrictionCorrelationWarning
    assert "FrictionCorrelationWarning" in friction.__all__


@pytest.mark.parametrize(
    "reynolds, diameter, radius",
    [
        (50.0, 0.01, 0.1),
        (100.0, 0.01, 0.1),
        (np.nextafter(150_000.0, np.inf), 0.01, 0.1),
        (200_000.0, 0.01, 0.1),
        (4_000.0, 0.14, 0.5),
        (4_000.0, 0.14, -0.5),
        (4_000.0, 0.15, 0.5),
    ],
)
@pytest.mark.parametrize("entry_point", ["multiplier", "gradient"])
def test_out_of_range_curved_correlations_warn_but_continue(
    reynolds, diameter, radius, entry_point,
) -> None:
    roughness = 2e-5
    curved = _curved_reference(reynolds, diameter, radius, roughness)
    with pytest.warns(friction.FrictionCorrelationWarning):
        if entry_point == "multiplier":
            actual = friction.curvature_multiplier(
                reynolds, diameter, radius, roughness=roughness, model="fluids",
            )
            expected = curved / _straight_reference(reynolds, diameter, roughness)
        else:
            actual = friction.friction_pressure_gradient(
                reynolds, diameter, 850.0, 3.0,
                roughness=roughness, R_curv=radius, model="fluids",
            )
            expected = -curved * 850.0 * 3.0**2 / (2.0 * diameter)
    assert np.isfinite(actual)
    assert actual == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize(
    "reynolds, diameter, radius",
    [
        (np.nextafter(100.0, np.inf), 0.01, 0.1),
        (150_000.0, 0.01, 0.1),
        (4_000.0, np.nextafter(0.14, 0.0), 0.5),
        (4_000.0, 0.01, -0.1),
        (10_000.0, 0.01, 0.1),
        (100_000.0, 0.01, 0.1),
    ],
)
def test_in_range_curved_correlations_do_not_warn(reynolds, diameter, radius):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        friction.curvature_multiplier(
            reynolds, diameter, radius, roughness=2e-5, model="fluids",
        )
        friction.friction_pressure_gradient(
            reynolds, diameter, 850.0, 3.0,
            roughness=2e-5, R_curv=radius, model="fluids",
        )


@pytest.mark.parametrize("reynolds", [50.0, 100.0, 200_000.0])
@pytest.mark.parametrize("radius, enabled", [(np.inf, True), (-np.inf, True),
                                             (0.01, False)])
def test_curved_validity_warnings_do_not_apply_to_straight_or_disabled_flow(
    reynolds, radius, enabled,
) -> None:
    # These are curved-correlation limits, not straight-pipe validity claims.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        straight = friction.darcy_friction_factor(
            reynolds, 0.01, 2e-5, model="fluids",
        )
        if np.isinf(radius):
            assert friction.curvature_multiplier(
                reynolds, 0.01, radius, roughness=2e-5, model="fluids",
            ) == 1.0
        actual = friction.friction_pressure_gradient(
            reynolds, 0.01, 850.0, 3.0, roughness=2e-5,
            R_curv=radius, curvature_correction=enabled, model="fluids",
        )
    assert actual == pytest.approx(-straight * 850.0 * 3.0**2 / (2.0 * 0.01))


@pytest.mark.parametrize("reynolds", [3_500.0, 1e4, 1e6, 1e9, 1e12])
@pytest.mark.parametrize("diameter", [1e-4, 0.01, 1.0])
@pytest.mark.parametrize("relative_roughness", [1e-12, 1e-8, 1e-4, 0.01, 0.1])
def test_native_rough_colebrook_is_positive_and_satisfies_residual(
    reynolds, diameter, relative_roughness,
) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        factor = friction.darcy_friction_factor(
            reynolds, diameter, relative_roughness * diameter, model="native",
        )
    assert np.isfinite(factor)
    assert factor > 0.0
    inverse_sqrt = 1.0 / np.sqrt(factor)
    residual = inverse_sqrt + 2.0 * np.log10(
        relative_roughness / 3.71 + 2.51 * inverse_sqrt / reynolds
    )
    assert residual == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("roughness", [None, 0.0, 1e-8, 1e-4, 1e-3])
@pytest.mark.parametrize("reynolds", [1_000.0, 2_299.0, 2_300.0, 2_301.0,
                                     2_900.0, 3_499.0, 3_500.0, 100_000.0])
def test_native_laminar_transition_and_smooth_petukhov_are_preserved(
    reynolds, roughness,
) -> None:
    diameter = 0.01
    laminar = 64.0 / reynolds
    if roughness:
        # Clamond uses 3.7; rescale e/D to independently reference native 3.71.
        turbulent = fluids.friction_factor(
            Re=reynolds, eD=roughness / diameter * 3.7 / 3.71,
            Method="Clamond", Darcy=True,
        )
    else:
        turbulent = (0.79 * np.log(reynolds) - 1.64) ** -2
    weight = np.clip((reynolds - 2_300.0) / 1_200.0, 0.0, 1.0)
    expected = (1.0 - weight) * laminar + weight * turbulent
    assert friction.darcy_friction_factor(
        reynolds, diameter, roughness, model="native",
    ) == pytest.approx(expected, rel=1e-10)
    assert friction.darcy_friction_factor(
        reynolds, diameter, roughness,
    ) == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize(
    "reynolds, diameter, radius, group",
    [
        # Dh/4 == |R| makes the group exact at and around the threshold.
        (5.999, 0.04, 0.01, 5.999),
        (6.0, 0.04, 0.01, 6.0),
        (6.001, 0.04, 0.01, 6.001),
        (100_000.0, 0.01, 0.1, 62.5),
        (100_000.0, 0.01, 10.0, 0.00625),
    ],
)
@pytest.mark.parametrize("roughness", [None, 0.0, 2e-5])
@pytest.mark.parametrize("sign", [-1, 1])
def test_native_curvature_remains_ito_and_independent_of_roughness(
    reynolds, diameter, radius, group, roughness, sign,
) -> None:
    radius *= sign
    expected = 1.0 if group <= 6.0 else group**0.05
    assert friction.curvature_multiplier(
        reynolds, diameter, radius, roughness=roughness, model="native",
    ) == pytest.approx(expected, rel=1e-12)
    assert friction.curvature_multiplier(
        reynolds, diameter, radius, roughness=roughness,
    ) == pytest.approx(expected, rel=1e-12)
    straight = friction.darcy_friction_factor(reynolds, diameter, roughness)
    assert friction.friction_pressure_gradient(
        reynolds, diameter, 850.0, 3.0, roughness=roughness, R_curv=radius,
    ) == pytest.approx(-straight * expected * 850.0 * 3.0**2 / (2.0 * diameter))
