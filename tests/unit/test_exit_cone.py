"""Geometry contracts for the volute-to-flange exit cone."""

import math

import numpy as np
import pytest

from pyskyfire.ducts.exit_cone import ExitCone, resolve_exit_cone

START = np.array([0.06, 0.0, 0.0])  # volute outlet on the +x axis
OUTLET_DIAMETER = 0.03


def resolve(**kwargs):
    kwargs.setdefault("final_diameter", 0.05)
    kwargs.setdefault("arc_length_parameter", 2.0)
    cone = ExitCone(**kwargs)
    return resolve_exit_cone(cone, start=START, outlet_diameter=OUTLET_DIAMETER,
                             max_section_radius=OUTLET_DIAMETER / 2.0)


def azimuth(point):
    return math.atan2(point[1], point[0])


# -- construction validation -------------------------------------------------


def test_final_diameter_must_be_positive():
    with pytest.raises(ValueError, match="final_diameter"):
        ExitCone(final_diameter=0.0)


def test_first_angle_rejects_an_unknown_name():
    with pytest.raises(ValueError, match="tangential.*radial"):
        ExitCone(final_diameter=0.05, first_angle="sideways")


def test_first_angle_rejects_a_negative_number():
    with pytest.raises(ValueError, match="negative turn"):
        ExitCone(final_diameter=0.05, first_angle=-0.1)


def test_second_angle_is_bounded_to_a_right_angle():
    with pytest.raises(ValueError, match="second_angle"):
        ExitCone(final_diameter=0.05, second_angle=math.radians(100))


def test_arc_length_parameter_has_a_one_diameter_floor():
    with pytest.raises(ValueError, match="at least 1"):
        ExitCone(final_diameter=0.05, arc_length_parameter=0.5)


# -- tangential ---------------------------------------------------------------


def test_tangential_outlet_runs_straight_with_no_turn():
    result = resolve(first_angle="tangential")
    assert result.is_straight
    assert result.turn_angle == pytest.approx(0.0, abs=1e-12)
    assert result.first_angle == pytest.approx(0.0)
    # Tangential at azimuth zero is +y.
    assert result.end_direction == pytest.approx([0.0, 1.0, 0.0], abs=1e-12)
    assert result.end == pytest.approx(START + result.arc_length * result.end_direction)


def test_tangential_is_the_default():
    default = resolve()
    assert default.first_angle == pytest.approx(0.0)
    assert default.is_straight


# -- radial ---------------------------------------------------------------


def test_radial_outlet_points_along_the_true_radial_direction_at_the_flange():
    """The flange has moved round the volute, so 'radial' is evaluated there.

    This is the fixed-point condition the user specified: having turned by
    ``first`` from the fixed tangential direction and travelled to azimuth
    ``phi``, radial means ``first + phi == pi/2`` -- not that the flange
    normal equals the *starting* radial direction.
    """
    result = resolve(first_angle="radial")
    end_azimuth = azimuth(result.end)
    radial_at_flange = np.array([math.cos(end_azimuth), math.sin(end_azimuth), 0.0])
    assert result.end_direction == pytest.approx(radial_at_flange, abs=1e-8)
    assert result.first_angle + end_azimuth == pytest.approx(math.pi / 2.0, abs=1e-8)
    # Having moved forward, radial is closer than a right angle from tangential.
    assert 0.0 < result.first_angle < math.pi / 2.0


def test_radial_first_angle_matches_the_reported_radial_angle():
    result = resolve(first_angle="radial")
    assert result.first_angle == pytest.approx(result.radial_angle)


def test_a_longer_arc_reaches_further_round_so_radial_needs_less_turn():
    """More azimuthal travel is available, so less of it has to come from turn."""
    short = resolve(first_angle="radial", arc_length_parameter=1.2)
    long = resolve(first_angle="radial", arc_length_parameter=4.0)
    assert long.first_angle < short.first_angle


# -- explicit numeric first angle ---------------------------------------------


def test_an_intermediate_first_angle_lies_between_the_named_extremes():
    tangential = resolve(first_angle="tangential")
    radial = resolve(first_angle="radial")
    halfway = radial.radial_angle / 2.0
    result = resolve(first_angle=halfway)
    assert result.first_angle == pytest.approx(halfway)
    assert not result.is_straight
    assert tangential.turn_angle < result.turn_angle < radial.turn_angle


def test_a_first_angle_beyond_radial_is_rejected():
    radial = resolve(first_angle="radial")
    with pytest.raises(ValueError, match="exceeds the radial angle"):
        resolve(first_angle=radial.radial_angle * 1.5)


def test_a_first_angle_exactly_at_the_radial_bound_is_accepted():
    radial = resolve(first_angle="radial")
    result = resolve(first_angle=radial.radial_angle)
    assert result.end_direction == pytest.approx(radial.end_direction, abs=1e-6)


# -- second angle (elevation) --------------------------------------------------


@pytest.mark.parametrize("degrees", [30.0, -45.0, 89.0])
def test_second_angle_elevates_the_outlet_by_exactly_its_sine(degrees):
    result = resolve(first_angle="tangential", second_angle=math.radians(degrees))
    assert result.end_direction[2] == pytest.approx(math.sin(math.radians(degrees)))
    assert np.sign(result.end[2] - START[2]) == np.sign(degrees)


def test_second_angle_leaves_the_horizontal_component_tangential_when_first_is_zero():
    result = resolve(first_angle="tangential", second_angle=math.radians(40))
    horizontal = np.array([result.end_direction[0], result.end_direction[1], 0.0])
    horizontal /= np.linalg.norm(horizontal)
    assert horizontal == pytest.approx([0.0, 1.0, 0.0], abs=1e-8)


def test_second_angle_applies_on_top_of_a_radial_first_angle():
    """A radial-plus-elevation outlet stays in the meridional half-plane."""
    result = resolve(first_angle="radial", second_angle=math.radians(20))
    end_azimuth = azimuth(result.end)
    tangential_at_flange = np.array([-math.sin(end_azimuth), math.cos(end_azimuth), 0.0])
    # The direction has no component along the flange's own tangential axis,
    # i.e. it lies in the meridional plane through the axis and the flange.
    assert np.dot(result.end_direction, tangential_at_flange) == pytest.approx(0.0, abs=1e-6)


# -- arc length and self-intersection -----------------------------------------


def test_arc_length_scales_with_the_outlet_diameter():
    result = resolve(arc_length_parameter=1.8)
    assert result.arc_length == pytest.approx(1.8 * OUTLET_DIAMETER)


def test_a_bend_too_tight_for_its_own_duct_is_rejected():
    with pytest.raises(ValueError, match="fold through itself"):
        resolve(first_angle="radial", arc_length_parameter=1.0)


def test_frames_start_at_the_volute_outlet_and_end_at_the_flange():
    result = resolve(first_angle="radial", arc_length_parameter=3.0)
    points, u, v = result.frames(n_stations=9)
    assert points[0] == pytest.approx(result.start)
    assert points[-1] == pytest.approx(result.end, abs=1e-8)
    # u, v, and the local tangent form a right-handed orthonormal frame throughout.
    tangent = np.gradient(points, axis=0)
    tangent /= np.linalg.norm(tangent, axis=1, keepdims=True)
    for i in range(len(points)):
        assert np.dot(u[i], v[i]) == pytest.approx(0.0, abs=1e-8)
        assert np.linalg.norm(u[i]) == pytest.approx(1.0)
        assert np.linalg.norm(v[i]) == pytest.approx(1.0)


def test_frames_do_not_twist_a_straight_duct():
    result = resolve(first_angle="tangential")
    _, u, v = result.frames(n_stations=5)
    assert np.allclose(u, u[0])
    assert np.allclose(v, v[0])
