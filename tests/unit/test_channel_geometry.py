"""Analytic contracts for channel sections and placement."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen.channel_placement import InternalPlacement, SurfacePlacement
from pyskyfire.regen.cross_section import (
    CrossSectionRounded,
    CrossSectionRoundedInternal,
    CrossSectionSquared,
    SectionProfiles,
)


def _profiles(*, radius=0.1, height=0.025, theta=0.2, wall=0.001):
    return SectionProfiles(
        h=np.atleast_1d(float(height)),
        theta=np.atleast_1d(float(theta)),
        t_wall=np.atleast_1d(float(wall)),
        centerline=np.array([[0.0, radius, 0.0]]),
    )


def test_section_profiles_require_consistent_station_shapes() -> None:
    with pytest.raises(ValueError, match="theta length"):
        SectionProfiles(
            h=np.ones(2),
            theta=np.ones(1),
            t_wall=np.ones(2),
            centerline=np.zeros((2, 3)),
        )


def test_squared_section_matches_annular_sector_geometry() -> None:
    profiles = _profiles(height=0.004)
    section = CrossSectionSquared(blockage_ratio=0.25)
    effective_theta = 0.2 * 0.75
    inner = 0.1 + 0.001
    outer = inner + 0.004
    expected_area = 0.5 * effective_theta * (outer**2 - inner**2)
    expected_perimeter = effective_theta * (inner + outer) + 2.0 * 0.004

    assert section.A_coolant(profiles)[0] == pytest.approx(expected_area)
    assert section.Dh_coolant(profiles)[0] == pytest.approx(
        4.0 * expected_area / expected_perimeter
    )
    assert section.P_thermal(profiles)[0] == pytest.approx(0.1 * effective_theta)


def test_more_coolant_side_convection_reduces_section_resistance() -> None:
    profiles = _profiles()
    section = CrossSectionRounded()

    low_h = section.R_coolant_per_len(profiles, np.array([1_000.0]), 20.0)
    high_h = section.R_coolant_per_len(profiles, np.array([2_000.0]), 20.0)

    assert high_h[0] < low_h[0]
    assert np.all(np.isfinite([low_h[0], high_h[0]]))


def test_internal_rounded_section_doubles_contact_perimeters() -> None:
    profiles = _profiles()
    external = CrossSectionRounded()
    internal = CrossSectionRoundedInternal()

    np.testing.assert_allclose(internal.A_coolant(profiles), external.A_coolant(profiles))
    np.testing.assert_allclose(internal.Dh_coolant(profiles), external.Dh_coolant(profiles))
    np.testing.assert_allclose(internal.P_thermal(profiles), 2.0 * external.P_thermal(profiles))
    np.testing.assert_allclose(internal.P_coolant(profiles), 2.0 * external.P_coolant(profiles))


def test_surface_placement_offsets_the_flow_centroid_and_honors_helix_angle() -> None:
    contour = SimpleNamespace(r=lambda x: 0.1, dr_dx=lambda x: 0.2)
    placement = SurfacePlacement(60, helix_angle=np.radians(30.0))

    radius = placement.compute_flow_centerline_radius(
        0.0, contour, wall_thickness=0.002, channel_height=0.004
    )
    dtheta_dx = placement.centerline_dtheta_dx(0.0, contour, radius, 0.2)

    assert radius == pytest.approx(0.104)
    assert radius * dtheta_dx / np.sqrt(1.0 + 0.2**2) == pytest.approx(
        np.tan(np.radians(30.0))
    )


def test_internal_placement_preserves_leaf_and_stack_metadata() -> None:
    placement = InternalPlacement(
        n_channel_positions=12,
        n_channels_per_leaf=3,
        channel_width=lambda x: 0.01,
    )

    assert placement.channel_count() == 12
    assert placement.n_channels_per_leaf == 3
    assert not placement.occludes
