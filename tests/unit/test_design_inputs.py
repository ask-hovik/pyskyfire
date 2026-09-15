"""Contracts for keyword-dispatched design-input functions.

A design input declares the quantities it needs by parameter name. Which names
are legal depends on the slot, because the geometry build order decides what is
already known when a slot is evaluated.
"""

from pathlib import Path

import numpy as np
import pytest

from pyskyfire.regen.design_inputs import (
    DesignContext,
    DesignInputError,
    SLOT_KEYS,
    TIER_A,
    TIER_B,
    TIER_C,
    VOCABULARY,
    contour_providers,
    evaluate,
    evaluate_profile,
    layout_providers,
    make_evaluator,
    section_providers,
)


class _Contour:
    """Cone from r=0.1 at x=-1 to a throat of r=0.05 at x=0, opening to 0.2."""

    xs = np.array([-1.0, 0.0, 1.0])
    x_t = 0.0
    r_t = 0.05
    r_c = 0.1
    r_e = 0.2

    def r(self, x):
        return 0.1 - 0.05 * (x + 1.0) if x < 0.0 else 0.05 + 0.15 * x

    def dr_dx(self, x):
        return 0.05 if x < 0.0 else 0.15

    def normal_angle(self, x):
        return np.arctan(self.dr_dx(x))


def _context(x=0.5, **overrides):
    providers = contour_providers(_Contour(), x)
    providers.update(
        layout_providers(
            n_channel_positions=60,
            n_channels=60,
            n_channels_all_circuits=60.0,
            t_wall=0.001,
            t_wall_total=0.0015,
            t_hot_layer=0.0005,
        )
    )
    providers.update(
        section_providers(
            r=_Contour().r(x),
            r_base=_Contour().r(x) + 0.001,
            h=0.004,
            gamma=0.0,
            theta=2.0 * np.pi / 60.0,
            t_wall=0.001,
        )
    )
    providers.update(overrides)
    return DesignContext(providers)


# --- dispatch ---------------------------------------------------------------


def test_a_constant_needs_no_context() -> None:
    assert evaluate(2.5, _context(), slot="sidewall_half_thickness", label="test") == 2.5


def test_a_function_of_nothing_is_called_with_nothing() -> None:
    assert evaluate(lambda: 1e-3, _context(), slot="sidewall_half_thickness", label="test") == 1e-3


def test_a_function_receives_only_the_quantities_it_names() -> None:
    seen = {}

    def sidewall_half_thickness(pitch_mid, h):
        seen.update(pitch_mid=pitch_mid, h=h)
        return 0.25 * pitch_mid

    result = evaluate(sidewall_half_thickness, _context(), slot="sidewall_half_thickness", label="test")

    assert set(seen) == {"pitch_mid", "h"}
    assert result == pytest.approx(0.25 * seen["pitch_mid"])


def test_parameter_order_does_not_matter() -> None:
    context = _context()
    forwards = evaluate(lambda x, r: r - x, context, slot="sidewall_half_thickness", label="t")
    backwards = evaluate(lambda r, x: r - x, context, slot="sidewall_half_thickness", label="t")

    assert forwards == backwards


def test_a_kwargs_catch_all_receives_everything_available() -> None:
    captured = {}

    def sidewall_half_thickness(**quantities):
        captured.update(quantities)
        return 1e-3

    evaluate(sidewall_half_thickness, _context(), slot="sidewall_half_thickness", label="test")

    assert {"x", "r", "pitch_mid", "n_channels", "t_wall"} <= set(captured)


def test_evaluate_profile_walks_every_station() -> None:
    contexts = [_context(x) for x in (-0.5, 0.0, 0.5)]
    values = evaluate_profile(lambda x: 2.0 * x, contexts, slot="sidewall_half_thickness", label="t")

    np.testing.assert_allclose(values, [-1.0, 0.0, 1.0])


# --- tiers ------------------------------------------------------------------


def test_a_wall_thickness_may_only_ask_about_the_contour() -> None:
    make_evaluator(lambda r, x_frac: 0.001, slot="wall_thickness", label="t")

    with pytest.raises(DesignInputError, match="tier B"):
        make_evaluator(lambda t_wall: t_wall, slot="wall_thickness", label="t")


def test_a_wall_thickness_cannot_depend_on_the_wall_stack_it_defines() -> None:
    with pytest.raises(DesignInputError, match="would be circular"):
        make_evaluator(
            lambda t_wall_total: 0.001, slot="wall_thickness", label="t"
        )


def test_a_helix_angle_cannot_depend_on_the_pitch_its_own_angle_sets() -> None:
    # theta = 2*pi*cos(gamma)/n, so letting gamma depend on theta is a cycle.
    with pytest.raises(DesignInputError, match="tier C"):
        make_evaluator(lambda theta: theta, slot="helix_angle", label="t")

    make_evaluator(
        lambda circumference, n_channels_all_circuits: 0.1,
        slot="helix_angle",
        label="t",
    )


def test_a_sidewall_half_thickness_may_ask_for_anything() -> None:
    for name in VOCABULARY:
        assert name in SLOT_KEYS["sidewall_half_thickness"]


def test_an_unknown_quantity_is_reported_with_the_legal_ones() -> None:
    with pytest.raises(DesignInputError) as excinfo:
        make_evaluator(lambda radius: 1.0, slot="sidewall_half_thickness", label="mine")

    message = str(excinfo.value)
    assert "mine" in message
    assert "'radius' is not a design-input quantity" in message
    assert "pitch_mid" in message


def test_every_vocabulary_entry_belongs_to_exactly_one_tier() -> None:
    tiers = [set(TIER_A), set(TIER_B), set(TIER_C)]
    assert set().union(*tiers) == set(VOCABULARY)
    for first in range(len(tiers)):
        for second in range(first + 1, len(tiers)):
            assert not tiers[first] & tiers[second]


def test_every_slot_only_admits_vocabulary_names() -> None:
    for slot, keys in SLOT_KEYS.items():
        assert keys <= set(VOCABULARY), slot


# --- rejected signatures ----------------------------------------------------


def test_star_args_is_rejected_with_an_explanation() -> None:
    with pytest.raises(DesignInputError, match="declare the quantities"):
        make_evaluator(lambda *args: 0.0, slot="helix_angle", label="t")


def test_positional_only_parameters_are_rejected() -> None:
    def sidewall_half_thickness(x, /):
        return x

    with pytest.raises(DesignInputError, match="positional-only"):
        make_evaluator(sidewall_half_thickness, slot="sidewall_half_thickness", label="t")


def test_returning_an_array_is_rejected() -> None:
    evaluator = make_evaluator(
        lambda: np.zeros(3), slot="sidewall_half_thickness", label="t"
    )
    with pytest.raises(DesignInputError, match="must return a scalar"):
        evaluator(_context())


# --- providers --------------------------------------------------------------


def test_contour_quantities_are_derived_consistently() -> None:
    context = _context(0.5)

    assert context["r"] == pytest.approx(0.125)
    assert context["circumference"] == pytest.approx(2.0 * np.pi * 0.125)
    assert context["area_ratio"] == pytest.approx((0.125 / 0.05) ** 2)
    assert context["x_frac"] == pytest.approx(0.5)
    assert _context(-0.5)["x_frac"] == pytest.approx(-0.5)


def test_pitches_are_measured_at_the_radius_they_name() -> None:
    context = _context(0.5)
    theta = 2.0 * np.pi / 60.0

    # The base point is resolved geometry, not ``r + t_wall``; the context
    # reports whatever radius the placed section actually sits at.
    assert context["r_channel_base"] == pytest.approx(0.126)
    assert context["r_channel_top"] == pytest.approx(0.130)
    assert context["pitch"] == pytest.approx(0.125 * theta)
    assert context["pitch_base"] == pytest.approx(0.126 * theta)
    assert context["pitch_mid"] == pytest.approx(0.128 * theta)


def test_section_quantities_are_omitted_until_they_exist() -> None:
    partial = section_providers(r=0.1, h=0.004, gamma=0.0, theta=None, t_wall=0.001)

    assert set(partial) == {"gamma", "h", "r_channel_base", "r_channel_top"}
    # An unplaced section has no base point of its own yet.
    assert partial["r_channel_base"] == pytest.approx(0.1)


def test_a_missing_contour_is_reported_where_it_is_used() -> None:
    context = DesignContext(contour_providers(None, 0.0))

    assert context["x"] == 0.0
    with pytest.raises(DesignInputError, match="no contour is attached"):
        context["r"]


def test_a_provider_is_evaluated_at_most_once() -> None:
    calls = []

    context = DesignContext({"r": lambda: calls.append(1) or 0.1})

    assert context["r"] == 0.1
    assert context["r"] == 0.1
    assert len(calls) == 1


# --- documentation ----------------------------------------------------------


def test_the_glossary_documents_every_quantity_and_slot() -> None:
    """The explanation page promises an exhaustive glossary; hold it to that."""
    page = (
        Path(__file__).resolve().parents[2]
        / "docs"
        / "explanations"
        / "design-input-functions.md"
    ).read_text(encoding="utf-8")

    for name in VOCABULARY:
        assert f"| `{name}` |" in page, f"{name} is missing from the glossary"
    for slot in SLOT_KEYS:
        argument = "thickness" if slot == "wall_thickness" else slot
        assert f"`{argument}`" in page, f"{slot} is missing from the slot table"
