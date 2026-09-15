"""Local regen-to-duct friction wiring, without a full thermal solve."""

from types import SimpleNamespace
import warnings

import fluids
import numpy as np
import pytest

from pyskyfire.regen import coupled_solver
from pyskyfire.common import physics
from pyskyfire.regen.coupled_solver import CoupledHeatExchangerPhysics


def _transition_weight(reynolds, critical):
    end = critical * 3_500.0 / 2_300.0
    fraction = np.clip((reynolds - critical) / (end - critical), 0.0, 1.0)
    return fraction**2 * (3.0 - 2.0 * fraction)


def _straight_factor(reynolds, diameter, roughness):
    laminar = fluids.friction.friction_laminar(reynolds)
    weight = _transition_weight(reynolds, 2_040.0)
    if weight == 0.0:
        return laminar
    turbulent = fluids.friction.Clamond(reynolds, roughness / diameter)
    return (1.0 - weight) * laminar + weight * turbulent


def _curved_factor(reynolds, diameter, radius, roughness):
    coil_diameter = 2.0 * abs(radius)
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


@pytest.fixture
def friction_case():
    # Match the minimal physics fixtures in test_curvature_corrections.py.
    circuit = SimpleNamespace(
        n_channels=2,
        coolant_transport=SimpleNamespace(get_mu=lambda temperature, pressure: 1e-3),
        A_coolant=lambda x: 1e-4,
        Dh_coolant=lambda x: 0.01,
        radius_of_curvature=lambda x: -0.1,
        roughness=2e-5,
        ds_dx=lambda x: 1.0 + x,
    )
    helper = object.__new__(CoupledHeatExchangerPhysics)
    helper.thrust_chamber = SimpleNamespace(cooling_circuits=[circuit])
    helper.circuit_index = 0
    helper.boundary_conditions = SimpleNamespace(mdot_coolant=0.08)
    helper.coolant_state = SimpleNamespace(saturation=lambda pressure: None)
    helper.bulk_density = lambda *args: 800.0
    helper.friction_model = "fluids"
    helper.pressure_curvature_correction = True
    return helper, circuit


@pytest.mark.parametrize("reynolds", [4_000.0, 10_000.0, 100_000.0])
@pytest.mark.parametrize("roughness_kind", ["scalar", "callable", "none"])
@pytest.mark.parametrize("enabled", [False, True])
def test_regen_fluids_gradient_uses_local_roughness_and_arc_length(
    friction_case, reynolds, roughness_kind, enabled,
) -> None:
    helper, circuit = friction_case
    helper.boundary_conditions.mdot_coolant = reynolds * 2e-5
    helper.pressure_curvature_correction = enabled
    evaluated_at = []

    def roughness_at(x):
        evaluated_at.append(x)
        return 1e-5 * (1.0 + x)

    circuit.roughness = {
        "scalar": 2e-5, "callable": roughness_at, "none": None,
    }[roughness_kind]
    for x in (0.25, 0.75):
        local_roughness = {
            "scalar": 2e-5, "callable": 1e-5 * (1.0 + x), "none": 0.0,
        }[roughness_kind]
        if enabled:
            factor = _curved_factor(
                reynolds, 0.01, 0.1, local_roughness,
            )
        else:
            factor = _straight_factor(reynolds, 0.01, local_roughness)
        velocity = helper.boundary_conditions.mdot_coolant / (2.0 * 800.0 * 1e-4)
        expected = -factor / 0.01 * 800.0 * velocity**2 / 2.0 * (1.0 + x)
        actual = helper.coolant_friction_rate(x, T_cool=100.0, p_cool=1e6)
        assert actual == pytest.approx(expected, rel=1e-12)

    if roughness_kind == "callable":
        assert evaluated_at == [0.25, 0.75]


@pytest.mark.parametrize("roughness_kind", ["scalar", "callable", "none"])
def test_regen_passes_same_resolved_roughness_to_both_factors(
    monkeypatch, friction_case, roughness_kind,
) -> None:
    helper, circuit = friction_case
    evaluations = []
    calls = []
    x = 0.75

    def roughness_at(station):
        evaluations.append(station)
        return 1e-5 * (1.0 + station)

    circuit.roughness = {
        "scalar": 2e-5, "callable": roughness_at, "none": None,
    }[roughness_kind]
    local_roughness = {
        "scalar": 2e-5, "callable": 1.75e-5, "none": 0.0,
    }[roughness_kind]

    def straight(Re, Dh, roughness=0.0, *, model="native"):
        calls.append("straight")
        assert Re == pytest.approx(4_000.0)
        assert Dh == 0.01
        assert roughness == pytest.approx(local_roughness)
        assert model == "fluids"
        return 0.04

    def multiplier(Re, Dh, R, *, roughness=0.0, model="native"):
        calls.append("curved")
        assert Re == pytest.approx(4_000.0)
        assert (Dh, R) == (0.01, -0.1)
        assert roughness == pytest.approx(local_roughness)
        assert model == "fluids"
        return 0.75

    monkeypatch.setattr(coupled_solver.ducts_friction, "darcy_friction_factor", straight)
    monkeypatch.setattr(coupled_solver.ducts_friction, "curvature_multiplier", multiplier)
    actual = helper.coolant_friction_rate(x, T_cool=100.0, p_cool=1e6)

    assert calls == ["straight", "curved"]
    assert evaluations == ([x] if roughness_kind == "callable" else [])
    assert actual == pytest.approx(-0.04 * 0.75 / 0.01 * 800.0 * 0.5**2 / 2.0 * 1.75)


@pytest.mark.parametrize("model", ["native", "fluids"])
def test_regen_disabled_curvature_does_not_evaluate_radius_or_multiplier(
    monkeypatch, friction_case, model,
) -> None:
    helper, circuit = friction_case
    helper.friction_model = model
    helper.pressure_curvature_correction = False

    def unexpected(*args, **kwargs):
        pytest.fail("disabled pressure curvature correction was evaluated")

    circuit.radius_of_curvature = unexpected
    monkeypatch.setattr(coupled_solver.ducts_friction, "curvature_multiplier", unexpected)
    monkeypatch.setattr(physics, "phi_curv_friction", unexpected)
    straight = coupled_solver.ducts_friction.darcy_friction_factor(
        4_000.0, 0.01, circuit.roughness, model=model,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        actual = helper.coolant_friction_rate(0.75, T_cool=100.0, p_cool=1e6)
    assert actual == pytest.approx(-straight / 0.01 * 800.0 * 0.5**2 / 2.0 * 1.75)


@pytest.mark.parametrize("reynolds", [50.0, 200_000.0])
@pytest.mark.parametrize("radius, enabled", [(np.inf, True), (0.01, False)])
def test_regen_skips_curved_warnings_for_straight_or_disabled_flow(
    friction_case, reynolds, radius, enabled,
) -> None:
    helper, circuit = friction_case
    helper.boundary_conditions.mdot_coolant = reynolds * 2e-5
    helper.pressure_curvature_correction = enabled
    circuit.radius_of_curvature = lambda x: radius
    straight = _straight_factor(reynolds, 0.01, 2e-5)
    velocity = helper.boundary_conditions.mdot_coolant / (2.0 * 800.0 * 1e-4)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        actual = helper.coolant_friction_rate(0.75, T_cool=100.0, p_cool=1e6)
    assert actual == pytest.approx(
        -straight / 0.01 * 800.0 * velocity**2 / 2.0 * 1.75,
    )
