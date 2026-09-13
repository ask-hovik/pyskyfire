"""Tests for the measurable-height interpretation of rounded channels."""

import numpy as np
import pytest

from pyskyfire.regen.cross_section import CrossSectionRounded, SectionProfiles


def _profiles(*, measured_height: float, radius=0.1, theta=0.2, wall=0.001):
    return SectionProfiles(
        h=np.array([measured_height]),
        theta=np.array([theta]),
        t_wall=np.array([wall]),
        cross_section_base_point_path=np.array([[0.0, radius, 0.0]]),
    )


def _measured_height(*, base_height, radius=0.1, theta=0.2, wall=0.001):
    # Both lobes are tangent to the sector rays, so raising the outer lobe
    # centre by ``base_height`` grows its radius by ``base_height*sin(theta/2)``
    # and the whole passage by ``base_height*(1 + sin(theta/2))``.
    half_theta = theta / 2.0
    inner_arc_radius = radius * np.sin(half_theta) - wall
    return 2.0 * inner_arc_radius + base_height * (1.0 + np.sin(half_theta))


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


@pytest.mark.parametrize("section_name", ["CrossSectionRounded", "CrossSectionRoundedInternal"])
@pytest.mark.parametrize("base", [0.0, 0.004, 0.02])
def test_base_input_matches_full_input_and_legacy_section_equations(section_name, base):
    from dataclasses import replace
    from pyskyfire.regen import cross_section

    cls = getattr(cross_section, section_name)
    legacy_input = cls(height_input="base")
    default_input = cls()
    full = _measured_height(base_height=base)
    prof_full = _profiles(measured_height=full)
    prof_base = replace(prof_full, h=np.array([base]))
    np.testing.assert_allclose(legacy_input.base_height(prof_base), base)
    np.testing.assert_allclose(legacy_input.full_height(prof_base), full)
    np.testing.assert_allclose(default_input.full_height(prof_full), full)
    for method in ("A_coolant", "Dh_coolant", "P_thermal", "P_coolant"):
        np.testing.assert_allclose(
            getattr(legacy_input, method)(prof_base), getattr(default_input, method)(prof_full),
        )
    # These are the direct prof.h equations from ca6d36a, before the full-
    # height conversion. The option restores that parameterization.
    q_inner = 0.1 * np.sin(0.1) - 0.001
    q_outer = (0.1 + base) * np.sin(0.1) - 0.001
    side = base * np.cos(0.1)
    beta, alpha = (np.pi - 0.2) / 2, (np.pi + 0.2) / 2
    area = beta * q_inner**2 + alpha * q_outer**2 + (q_inner + q_outer) * side
    perimeter = 2 * beta * q_inner + 2 * alpha * q_outer + 2 * side
    np.testing.assert_allclose(legacy_input.A_coolant(prof_base), area)
    np.testing.assert_allclose(legacy_input.Dh_coolant(prof_base), 4 * area / perimeter)
    if section_name == "CrossSectionRounded":
        np.testing.assert_allclose(
            legacy_input.R_coolant_per_len(prof_base, np.array([1000.0]), 350.0),
            default_input.R_coolant_per_len(prof_full, np.array([1000.0]), 350.0),
        )


@pytest.mark.parametrize("height", [-0.001, np.nan, np.inf])
def test_base_input_rejects_invalid_heights(height):
    section = CrossSectionRounded(height_input="base")
    with pytest.raises(ValueError, match="height"):
        section.full_height(_profiles(measured_height=height))


def test_height_input_validation_and_serialization():
    import pickle

    with pytest.raises(ValueError, match="height_input"):
        CrossSectionRounded(height_input="internal")
    for mode in ("full", "base"):
        section = pickle.loads(pickle.dumps(CrossSectionRounded(24, height_input=mode)))
        assert section.height_input == mode
        assert section.n_points == 24
    # Saved objects from before this option keep the full-height behavior.
    old = object.__new__(CrossSectionRounded)
    old.n_points = 16
    assert pickle.loads(pickle.dumps(old)).height_input == "full"
