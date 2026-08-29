"""Tests for the measurable-height interpretation of rounded channels."""

import numpy as np
import pytest

from pyskyfire.regen.cross_section import CrossSectionRounded, SectionProfiles


def _profiles(*, measured_height: float, radius=0.1, theta=0.2, wall=0.001):
    return SectionProfiles(
        h=np.array([measured_height]),
        theta=np.array([theta]),
        t_wall=np.array([wall]),
        centerline=np.array([[0.0, radius, 0.0]]),
    )


def _measured_height(*, base_height, radius=0.1, theta=0.2, wall=0.001):
    half_theta = theta / 2.0
    inner_arc_radius = radius * np.sin(half_theta) - wall
    return (
        2.0 * inner_arc_radius
        + base_height * (np.sin(half_theta) + np.cos(half_theta))
    )


def test_measurable_height_is_converted_back_to_base_height():
    section = CrossSectionRounded()
    expected_base_height = 0.004
    profiles = _profiles(
        measured_height=_measured_height(base_height=expected_base_height)
    )

    np.testing.assert_allclose(section.base_height(profiles), expected_base_height)


def test_minimum_measurable_height_produces_a_circle():
    section = CrossSectionRounded()
    profiles = _profiles(measured_height=_measured_height(base_height=0.0))
    arc_radius = 0.1 * np.sin(0.2 / 2.0) - 0.001

    np.testing.assert_allclose(section.base_height(profiles), 0.0, atol=1e-15)
    np.testing.assert_allclose(section.A_coolant(profiles), np.pi * arc_radius**2)
    np.testing.assert_allclose(section.Dh_coolant(profiles), 2.0 * arc_radius)


def test_geometry_matches_old_equations_after_height_conversion():
    section = CrossSectionRounded()
    radius = 0.1
    theta = 0.2
    wall = 0.001
    base_height = 0.004
    profiles = _profiles(
        measured_height=_measured_height(
            base_height=base_height,
            radius=radius,
            theta=theta,
            wall=wall,
        ),
        radius=radius,
        theta=theta,
        wall=wall,
    )

    beta = (np.pi - theta) / 2.0
    alpha = np.pi - beta
    q_inner = radius * np.sin(theta / 2.0) - wall
    q_outer = (radius + base_height) * np.sin(theta / 2.0) - wall
    side_length = base_height * np.cos(theta / 2.0)
    expected_area = (
        beta * q_inner**2
        + alpha * q_outer**2
        + (q_inner + q_outer) * side_length
    )
    expected_perimeter = (
        2.0 * beta * q_inner
        + 2.0 * alpha * q_outer
        + 2.0 * side_length
    )

    np.testing.assert_allclose(section.A_coolant(profiles), expected_area)
    np.testing.assert_allclose(
        section.Dh_coolant(profiles), 4.0 * expected_area / expected_perimeter
    )


def test_height_below_rounded_minimum_is_rejected():
    section = CrossSectionRounded()
    minimum = _measured_height(base_height=0.0)
    profiles = _profiles(measured_height=minimum - 1e-4)

    with pytest.raises(ValueError, match="below the geometric minimum"):
        section.A_coolant(profiles)
