"""Analytic contracts for channel sections and placement."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen.channel_placement import InternalPlacement, SurfacePlacement
from pyskyfire.regen.cross_section import (
    CrossSectionRounded,
    CrossSectionSquared,
    SectionProfiles,
)


def _profiles(*, radius=0.1, height=0.025, theta=0.2, wall=0.001):
    return SectionProfiles(
        h=np.atleast_1d(float(height)),
        theta=np.atleast_1d(float(theta)),
        t_wall=np.atleast_1d(float(wall)),
        cross_section_base_point_path=np.array([[0.0, radius, 0.0]]),
    )


def test_section_profiles_require_consistent_station_shapes() -> None:
    with pytest.raises(ValueError, match="theta length"):
        SectionProfiles(
            h=np.ones(2),
            theta=np.ones(1),
            t_wall=np.ones(2),
            cross_section_base_point_path=np.zeros((2, 3)),
        )


def test_squared_section_removes_a_constant_width_rib_from_the_pitch() -> None:
    profiles = _profiles(height=0.004)
    rib = 0.0015
    section = CrossSectionSquared(sidewall_half_thickness=rib / 2)
    # The base point is the channel floor, so the hot wall lies t_wall inside it.
    inner = 0.1
    outer = inner + 0.004
    # Straight rib faces make these the first-order forms; the section metrics
    # integrate the real boundary, which differs by about one part in a million.
    expected_area = 0.5 * 0.2 * (outer**2 - inner**2) - rib * 0.004
    expected_perimeter = 0.2 * (inner + outer) - 2.0 * rib + 2.0 * 0.004

    assert section.A_coolant(profiles)[0] == pytest.approx(expected_area, rel=1e-5)
    assert section.Dh_coolant(profiles)[0] == pytest.approx(
        4.0 * expected_area / expected_perimeter, rel=1e-5
    )
    assert section.P_thermal(profiles)[0] == pytest.approx((inner - 0.001) * 0.2 - rib)
    assert section.P_coolant(profiles)[0] == pytest.approx(inner * 0.2 - rib, rel=1e-5)


def test_squared_section_rib_keeps_its_width_at_every_radius() -> None:
    """A taller channel adds sector area but no extra rib blockage."""
    rib = 0.0015
    section = CrossSectionSquared(sidewall_half_thickness=rib / 2)
    inner = 0.1

    short = section.A_coolant(_profiles(height=0.004))[0]
    tall = section.A_coolant(_profiles(height=0.008))[0]

    added_sector = 0.5 * 0.2 * ((inner + 0.008) ** 2 - (inner + 0.004) ** 2)
    assert tall - short == pytest.approx(added_sector - rib * 0.004, rel=1e-5)


def test_squared_section_sidewall_half_thickness_may_be_a_design_input_function() -> None:
    profiles = _profiles(height=0.004)
    constant = CrossSectionSquared(sidewall_half_thickness=(0.1 * (0.1 + 0.002) * 0.2) / 2)
    proportional = CrossSectionSquared(sidewall_half_thickness=lambda pitch_mid: (0.1 * pitch_mid) / 2)

    assert proportional.A_coolant(profiles)[0] == pytest.approx(
        constant.A_coolant(profiles)[0]
    )


def test_squared_section_rejects_a_rib_wider_than_the_pitch() -> None:
    section = CrossSectionSquared(sidewall_half_thickness=(1.0) / 2)
    with pytest.raises(ValueError, match="leaves no room"):
        section.A_coolant(_profiles(height=0.004))


def test_squared_equivalent_wedge_encloses_the_same_area() -> None:
    profiles = _profiles(height=0.004)
    section = CrossSectionSquared(sidewall_half_thickness=(0.0015) / 2)

    wedge = section.equivalent_angular_width(profiles)[0]
    inner = 0.1
    outer = inner + 0.004

    assert 0.5 * wedge * (outer**2 - inner**2) == pytest.approx(
        section.A_coolant(profiles)[0]
    )
    assert section.open_fraction(profiles)[0] == pytest.approx(wedge / 0.2)


def test_more_coolant_side_convection_reduces_section_resistance() -> None:
    profiles = _profiles()
    section = CrossSectionRounded()

    low_h = section.R_coolant_per_len(profiles, np.array([1_000.0]), 20.0)
    high_h = section.R_coolant_per_len(profiles, np.array([2_000.0]), 20.0)

    assert high_h[0] < low_h[0]
    assert np.all(np.isfinite([low_h[0], high_h[0]]))


def test_surface_placement_winds_the_reference_path_at_the_requested_helix_angle() -> None:
    contour = SimpleNamespace(r=lambda x: 0.1, dr_dx=lambda x: 0.2)
    placement = SurfacePlacement(60, helix_angle=np.radians(30.0))

    dtheta_dx = placement.dtheta_dx(0.0, contour)

    # The helix angle is measured from the meridian on the ideal contour, so
    # every circuit that asks for it shares one winding regardless of where its
    # own section sits.
    assert 0.1 * dtheta_dx / np.sqrt(1.0 + 0.2**2) == pytest.approx(
        np.tan(np.radians(30.0))
    )


def test_internal_placement_preserves_leaf_and_stack_metadata() -> None:
    placement = InternalPlacement(
        n_channel_positions=12,
        n_channels_per_leaf=3,
        channel_angular_width=lambda: 0.01,
    )

    assert placement.channel_count() == 12
    assert placement.n_channels_per_leaf == 3
    assert placement.channel_angular_width(0.0) == pytest.approx(0.01)
    assert not hasattr(placement, "channel_width")
    assert not placement.occludes
