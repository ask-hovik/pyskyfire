"""The cold-side envelope of a channelled wall.

A manifold wrapped around the chamber sits against this surface, so what it
reports has to be the real outermost metal, and the area-equivalent radius has
to be the honest one for bookkeeping rather than either extreme.
"""

import numpy as np
import pytest

from pyskyfire.regen import OuterWall
from pyskyfire.regen.channel_placement import SurfacePlacement
from pyskyfire.regen.cross_section import (ChannelSection, CrossSectionRounded,
                                           CrossSectionSquared)
from pyskyfire.regen.thrust_chamber import CoolingCircuit, ThrustChamber


class _Cone:
    """Straight-walled contour, so the outward normal has a known tilt."""

    def __init__(self, radius=0.1, slope=0.0, length=1.0):
        self.radius = float(radius)
        self.slope = float(slope)
        self.xs = np.array([0.0, length])
        self.x_t = 0.0

    def r(self, x):
        return self.radius + self.slope*np.asarray(x, dtype=float)

    def dr_dx(self, x):
        return np.full_like(np.asarray(x, dtype=float), self.slope)


class _Wall:
    def __init__(self, thickness=0.0):
        self._thickness = float(thickness)

    def thickness(self, x):
        return np.full_like(np.asarray(x, dtype=float), self._thickness)


def chamber(section, *, slope=0.0, outer_walls=None, n_channel_positions=60,
            height=0.012, t_wall=0.001, radius=0.1):
    contour = _Cone(radius=radius, slope=slope)
    circuit = CoolingCircuit(
        name="jacket", contour=contour, cross_section=section, span=[0.0, 1.0],
        placement=SurfacePlacement(n_channel_positions=n_channel_positions),
        channel_height=lambda x: height, walls=[_Wall(t_wall)],
        coolant_transport=None, roughness=0.0, outer_walls=outer_walls,
    )
    return ThrustChamber(contour=contour, cooling_circuits=[circuit],
                         combustion_transport=None, n_nodes=40,
                         compute_gas=False), circuit


def rounded():
    return CrossSectionRounded()


def squared():
    return CrossSectionSquared(sidewall_half_thickness=0.0006)


def test_a_squared_roof_is_a_cylinder_so_both_radii_agree():
    """Its roof is an arc at constant radius and the rib tops sit on it too.

    There is no scalloping to average over, so the crest and the area-equivalent
    cylinder are the same surface.
    """
    wall = chamber(squared())[0].outer_wall
    assert wall.r(0.5, "envelope") == pytest.approx(wall.r(0.5, "equivalent"),
                                                    rel=1e-12)


def test_a_squared_roof_lands_where_the_stack_puts_it():
    section = squared()
    chamber_, circuit = chamber(section, radius=0.1, height=0.012, t_wall=0.001)
    prof = circuit.section_profiles_at(0.5)
    floor = float(prof.cross_section_base_point_path[0, 1])
    expected = floor + float(prof.h[0]) + float(section.roof_thickness(prof)[0])
    assert chamber_.outer_wall.r(0.5) == pytest.approx(expected, rel=1e-9)


def test_rounded_tubes_scallop_so_the_equivalent_sits_inside_the_crest():
    chamber_, circuit = chamber(rounded())
    wall = chamber_.outer_wall
    crest, equivalent = wall.r(0.5, "envelope"), wall.r(0.5, "equivalent")
    assert crest > equivalent
    # It must also sit outside the trough, where neighbouring tubes touch on
    # the sector ray, or it would not be an average of the real surface.
    prof = circuit.section_profiles_at(0.5)
    theta = float(prof.theta[0])
    centre = crest/(1.0 + np.sin(theta/2.0))
    trough = centre*np.cos(theta/2.0)
    assert trough < equivalent < crest


def test_the_equivalent_radius_encloses_the_area_the_tubes_really_enclose():
    """Cross-checked against a closed form rather than against itself.

    A rounded tube's exterior is tangent to both sector rays, so one sector is
    a triangle out to the two tangent points plus the circular segment bulging
    past their chord.
    """
    chamber_, circuit = chamber(rounded())
    wall = chamber_.outer_wall
    prof = circuit.section_profiles_at(0.5)
    theta = float(prof.theta[0])
    centre = wall.r(0.5, "envelope")/(1.0 + np.sin(theta/2.0))

    half_chord = centre*np.cos(theta/2.0)
    lobe = centre*np.sin(theta/2.0)
    sector = (0.5*half_chord**2*np.sin(theta)
              + 0.5*lobe**2*(np.pi + theta + np.sin(theta)))
    expected = np.sqrt(circuit.n_channels*sector/np.pi)
    assert wall.r(0.5, "equivalent") == pytest.approx(expected, rel=1e-6)


def test_outer_liners_push_the_wall_out_by_their_own_thickness():
    bare = chamber(rounded())[0].outer_wall
    lined = chamber(rounded(), outer_walls=[_Wall(0.0005)])[0].outer_wall
    for kind in ("envelope", "equivalent"):
        assert lined.r(0.5, kind) - bare.r(0.5, kind) == pytest.approx(0.0005,
                                                                       abs=1e-12)


def test_an_outer_liner_does_not_reach_the_channel_walls():
    """It lies outboard of the section, so the channel it closes is unchanged."""
    _, bare = chamber(rounded())
    _, lined = chamber(rounded(), outer_walls=[_Wall(0.0005)])
    for name in ("h", "theta", "t_wall"):
        assert np.allclose(getattr(bare.section_profiles_at(0.5), name),
                           getattr(lined.section_profiles_at(0.5), name))
    assert lined.channel_wall_thickness(0.5) == pytest.approx(
        bare.channel_wall_thickness(0.5))
    assert lined.outer_wall_thickness(0.5) == pytest.approx(0.0005)
    assert bare.outer_wall_thickness(0.5) == pytest.approx(0.0)


@pytest.mark.parametrize("slope", [0.0, 0.05, -0.05])
def test_the_normal_keeps_the_sign_of_the_wall_slope(slope):
    """Which way along the axis the wall tips is what a manifold has to know.

    A cone that widens downstream faces partly upstream and one that narrows
    faces partly downstream, and the two are not interchangeable.
    """
    wall = chamber(squared(), slope=slope)[0].outer_wall
    assert float(wall.dr_dx(0.5)) == pytest.approx(slope, abs=1e-6)
    normal = wall.normal(0.5)
    assert np.linalg.norm(normal) == pytest.approx(1.0)
    assert normal[0] == pytest.approx(-slope/np.hypot(1.0, slope), abs=1e-6)
    assert float(wall.normal_angle(0.5)) == pytest.approx(np.arctan(-slope),
                                                          abs=1e-6)


def test_a_cylindrical_wall_faces_straight_out():
    wall = chamber(squared())[0].outer_wall
    assert wall.normal(0.5) == pytest.approx([0.0, 1.0], abs=1e-9)
    assert float(wall.normal_angle(0.5)) == pytest.approx(0.0, abs=1e-9)


def test_a_section_with_no_exterior_is_refused_by_name():
    class _Analytic(ChannelSection):
        def A_coolant(self, prof):
            return prof.cross_section_base_point_path[:, 1]*prof.theta*prof.h

        def Dh_coolant(self, prof):
            return 2.0*prof.h

        def P_thermal(self, prof):
            return prof.cross_section_base_point_path[:, 1]*prof.theta

        def P_coolant(self, prof):
            return self.P_thermal(prof)

    chamber_, _ = chamber(_Analytic())
    with pytest.raises(ValueError, match="does not define an exterior boundary"):
        _ = chamber_.outer_wall


def test_the_wall_is_measured_once_and_kept():
    chamber_, _ = chamber(squared())
    assert chamber_.outer_wall is chamber_.outer_wall
    # Derived from resolved geometry, so it must not be carried in a pickle.
    assert "_outer_wall" not in chamber_.__getstate__()


def test_an_unknown_kind_names_the_two_that_exist():
    wall = chamber(squared())[0].outer_wall
    with pytest.raises(ValueError, match="kind must be one of"):
        wall.r(0.5, "outermost")


def test_a_hand_built_wall_refuses_an_envelope_inside_its_equivalent():
    with pytest.raises(ValueError, match="cannot fall inside"):
        OuterWall([0.0, 1.0], [0.10, 0.10], [0.11, 0.11])


@pytest.mark.parametrize("slope", [0.0, 0.05, -0.05, 0.2])
def test_a_squared_roof_stays_a_cylinder_however_the_wall_slopes(slope):
    """Pins which frame the sector is measured in.

    A squared roof is an arc at constant radius whatever the wall does, so its
    two radii have to coincide at every slope. On a sloped wall the section
    plane tilts away from the transverse one, and measuring the sector in the
    tilted plane instead splits them by microns -- so this failing is the sign
    that the boundary has been mapped out of the transverse frame it arrives in.
    """
    wall = chamber(squared(), slope=slope)[0].outer_wall
    assert wall.r(0.5, "envelope") == pytest.approx(wall.r(0.5, "equivalent"),
                                                    rel=1e-12)


def test_a_sloped_wall_still_brackets_the_rounded_crest_and_trough():
    # A widening cone has a bigger sector at its far end, so the lobes need
    # more room than the cylindrical fixture gives them.
    chamber_, circuit = chamber(rounded(), slope=0.05, height=0.020)
    wall = chamber_.outer_wall
    crest, equivalent = wall.r(0.5, "envelope"), wall.r(0.5, "equivalent")
    theta = float(circuit.section_profiles_at(0.5).theta[0])
    centre = crest/(1.0 + np.sin(theta/2.0))
    assert centre*np.cos(theta/2.0) < equivalent < crest


def _constant(value):
    """Height callable whose only parameter is ``x``.

    Design inputs are resolved by parameter name, so a default argument named
    after a channel quantity would be read as a request for it.
    """
    def height(x):
        return value
    return height


def interleaved(section_a, section_b, *, height_a=0.012, height_b=0.012,
                radius=0.1, t_wall=0.001):
    """Two circuits sharing one circumference, as a half pass and a return do."""
    contour = _Cone(radius=radius)
    circuits = [
        CoolingCircuit(
            name=name, contour=contour, cross_section=section, span=[0.0, 1.0],
            placement=SurfacePlacement(n_channel_positions=30),
            channel_height=_constant(height), walls=[_Wall(t_wall)],
            coolant_transport=None, roughness=0.0)
        for name, section, height in (("out", section_a, height_a),
                                      ("back", section_b, height_b))]
    return ThrustChamber(contour=contour, cooling_circuits=circuits,
                         combustion_transport=None, n_nodes=40,
                         compute_gas=False), circuits


def test_interleaved_circuits_share_the_circumference_they_sit_on():
    """Two circuits at the same station are side by side, not stacked.

    Each owns half the channels, so neither tiles the circle alone. Scaling one
    circuit's sector by its own channel count would size the wall as though it
    had the whole circumference, which lands it far inside the chamber.
    """
    chamber_, circuits = interleaved(squared(), squared())
    wall = chamber_.outer_wall
    assert sum(c.n_channels*float(c.section_profiles_at(0.5).theta[0])
               for c in circuits) == pytest.approx(2.0*np.pi, rel=1e-9)
    # Identical circuits must give the answer either one gives on its own.
    alone = OuterWall.from_circuit(circuits[0])
    for kind in ("envelope", "equivalent"):
        assert wall.r(0.5, kind) == pytest.approx(alone.r(0.5, kind), rel=1e-9)
    assert wall.r(0.5, "equivalent") > float(chamber_.contour.r(0.5))


def test_a_deeper_circuit_sets_the_crest_but_only_its_share_of_the_mean():
    """The envelope is a clearance; the equivalent is an area average."""
    chamber_, circuits = interleaved(squared(), squared(),
                                     height_a=0.012, height_b=0.018)
    wall = chamber_.outer_wall
    shallow = OuterWall.from_circuit(circuits[0])
    deep = OuterWall.from_circuit(circuits[1])

    assert wall.r(0.5, "envelope") == pytest.approx(deep.r(0.5, "envelope"),
                                                    rel=1e-9)
    assert (shallow.r(0.5, "equivalent") < wall.r(0.5, "equivalent")
            < deep.r(0.5, "equivalent"))


# -- handing the wall to a manifold section ------------------------------


def test_the_profile_is_measured_along_the_wall_not_along_the_axis():
    """A sloped wall offers more surface than its station span.

    The manifold section reserves the slit on that surface, so the profile has
    to hand it arc length rather than axial extent.
    """
    slope = 0.25
    wall = chamber(squared(), slope=slope)[0].outer_wall
    profile = wall.profile(0.5)
    low, high = profile.span
    axial_low, axial_high = wall.span
    stretch = np.hypot(1.0, slope)
    assert high - low == pytest.approx((axial_high - axial_low)*stretch, rel=1e-6)


def test_the_profile_sits_on_the_wall_it_came_from():
    wall = chamber(squared(), slope=0.05)[0].outer_wall
    profile = wall.profile(0.42)
    anchor = profile.at(0.0)
    assert anchor[0] == pytest.approx(0.42, abs=1e-6)
    assert anchor[1] == pytest.approx(float(wall.r(0.42)), rel=1e-9)
    for arc in np.linspace(*profile.span, 9):
        point = profile.at(arc)
        assert point[1] == pytest.approx(float(wall.r(point[0])), rel=1e-6)


def test_the_profile_can_follow_either_radius():
    wall = chamber(rounded())[0].outer_wall
    for kind in ("envelope", "equivalent"):
        assert wall.profile(0.5, kind).at(0.0)[1] == pytest.approx(
            float(wall.r(0.5, kind)), rel=1e-9)


def test_a_limited_reach_still_straddles_the_anchor():
    wall = chamber(squared())[0].outer_wall
    profile = wall.profile(0.5, reach=0.02)
    low, high = profile.span
    assert low < 0.0 < high
    assert max(-low, high) <= 0.02 + 1e-9


def test_an_anchor_outside_the_wall_is_refused():
    wall = chamber(squared())[0].outer_wall
    with pytest.raises(ValueError, match="strictly inside the wall's span"):
        wall.profile(wall.span[1])


def test_a_section_can_be_built_straight_off_a_chamber():
    from pyskyfire.ducts import HalfMoonSection

    chamber_, _ = chamber(squared(), slope=0.05)
    section = HalfMoonSection(chamber_.outer_wall.profile(0.5),
                              slit_height=0.014, fillet_radius=0.003,
                              aspect_ratio=1.2, max_aspect_ratio=5.0)
    areas = [float(section.metrics(size)[0][0]) for size in (0.014, 0.030, 0.050)]
    assert all(np.diff(areas) > 0.0)
    assert float(section.exposed_wall(0.030)[0]) == pytest.approx(0.030, rel=1e-12)
