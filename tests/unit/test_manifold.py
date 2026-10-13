"""Distributing and collecting rings.

A manifold is not a volute with a different law. It has no swirl to conserve,
it may run one way round or split and run both, and what it owes its channels
is an even share rather than a recovered head.
"""

import math

import numpy as np
import pytest
from matplotlib.path import Path

from pyskyfire.ducts import (ConstantAngularMomentum, ConstantVelocity,
                             EllipticalSlottedSection, ExitCone,
                             HalfMoonSection, InnerProfile, Manifold,
                             RadialPort, UniformArea)
from pyskyfire.ducts.volute_section import section_frames

MASS_FLOW, DENSITY = 20.0, 300.0


def wall(radius=0.20, reach=0.12, samples=501):
    s = np.linspace(-reach, reach, samples)
    return InnerProfile(s, np.column_stack((s, np.full_like(s, radius))))


def cone(half_angle, radius=0.20, reach=0.12, samples=501):
    """A sloping wall, so its sections lean and the two angle references differ."""
    s = np.linspace(-reach, reach, samples)
    return InnerProfile(s, np.column_stack((s*math.cos(half_angle),
                                            radius + s*math.sin(half_angle))))


def section(**kwargs):
    kwargs.setdefault("aspect_ratio", 1.2)
    kwargs.setdefault("max_aspect_ratio", 5.0)
    return HalfMoonSection(kwargs.pop("profile", wall()),
                           kwargs.pop("slit_height", 0.014),
                           kwargs.pop("fillet_radius", 0.003), **kwargs)


def tapered(velocity=30.0):
    return ConstantVelocity.for_duty(velocity=velocity, mass_flow=MASS_FLOW,
                                     density=DENSITY)


def manifold(**kwargs):
    kwargs.setdefault("section", section())
    kwargs.setdefault("area_law", tapered())
    kwargs.setdefault("n_stations", 61)
    if kwargs.get("branches", 2) == 1:
        kwargs.setdefault("rotation", "ccw")
    return Manifold(kwargs.pop("section"), kwargs.pop("area_law"), **kwargs)


# -- branches, not a direction ------------------------------------------


def test_a_split_manifold_has_no_rotation_to_name():
    """It is its own mirror image, so a sense would mean nothing."""
    assert manifold(branches=2).rotation is None
    with pytest.raises(ValueError, match="mirror-symmetric"):
        manifold(branches=2, rotation="ccw")


def test_a_single_branch_manifold_must_say_which_way_it_runs():
    with pytest.raises(ValueError, match="needs a rotation"):
        Manifold(section(), tapered(), branches=1, n_stations=21)
    assert manifold(branches=1, rotation="cw").rotation == "clockwise"


def test_more_branches_wait_for_more_ports():
    with pytest.raises(ValueError, match="branches must be 1 or 2"):
        manifold(branches=4)


@pytest.mark.parametrize("branches, degrees, share",
                         [(1, 360.0, 1.0), (2, 180.0, 0.5)])
def test_a_branch_covers_its_share_of_the_ring(branches, degrees, share):
    """Each branch gathers only what its own arc collects."""
    geometry = manifold(branches=branches).geometry
    assert math.degrees(geometry.theta[-1]) == pytest.approx(degrees)
    assert geometry.collected[0] == pytest.approx(0.0)
    assert geometry.collected[-1] == pytest.approx(share)


def test_splitting_the_ring_halves_the_area_the_port_needs():
    single = manifold(branches=1).geometry
    split = manifold(branches=2).geometry
    assert split.area[-1] == pytest.approx(single.area[-1]/2.0, rel=1e-6)


# -- the laws ------------------------------------------------------------


def test_a_tapered_ring_really_holds_its_velocity():
    """Checked as a velocity, not as the area the law was stated in."""
    ring = manifold(area_law=tapered(velocity=30.0))
    live = ~ring.geometry.clipped
    velocity = ring.collected_velocity(MASS_FLOW, DENSITY)
    assert velocity[live] == pytest.approx(30.0, rel=1e-6)


def test_a_uniform_ring_really_holds_its_area():
    ring = manifold(area_law=UniformArea(9.0e-4))
    assert ring.geometry.area == pytest.approx(9.0e-4, rel=1e-9)
    assert not ring.geometry.clipped.any()


def test_a_uniform_ring_speeds_up_as_it_fills():
    """The same passage carrying a growing share has to carry it faster."""
    ring = manifold(area_law=UniformArea(9.0e-4))
    velocity = ring.collected_velocity(MASS_FLOW, DENSITY)
    assert np.all(np.diff(velocity) > 0.0)
    assert velocity[0] == pytest.approx(0.0)


def test_a_duty_and_an_area_are_two_ways_of_saying_the_same_law():
    law = ConstantVelocity.for_duty(velocity=30.0, mass_flow=MASS_FLOW,
                                    density=DENSITY)
    assert law.full_area == pytest.approx(MASS_FLOW/(DENSITY*30.0))
    assert law.velocity(MASS_FLOW, DENSITY) == pytest.approx(30.0)


def test_the_law_runs_out_of_section_at_the_dead_point():
    """It wants zero area there, which no family can build."""
    geometry = manifold(area_law=tapered()).geometry
    assert geometry.clipped[0]
    assert not geometry.clipped[-1]
    assert geometry.dead_clip_angle > 0.0
    assert geometry.size[0] == pytest.approx(manifold().section.minimum_size)


def test_a_volute_law_is_refused_by_name():
    with pytest.raises(TypeError, match="must be a ManifoldLaw"):
        Manifold(section(), ConstantAngularMomentum(), n_stations=21)


@pytest.mark.parametrize("bad", [0.0, -1.0])
def test_the_laws_refuse_an_impossible_size(bad):
    with pytest.raises(ValueError, match="finite and positive"):
        UniformArea(bad)
    with pytest.raises(ValueError, match="finite and positive"):
        ConstantVelocity(bad)


# -- the ring it sweeps --------------------------------------------------


def test_a_split_ring_is_swept_whole_and_is_its_own_mirror():
    ring = manifold(branches=2).sweep(n_points=32)
    stations = len(ring)
    assert stations == 2*61 - 1
    # Station i and its opposite are reflections in the plane through the dead
    # point, which leaves the section shape alone and only moves it round.
    for index in range(1, 61):
        mirrored = ring[60 - index].copy()
        mirrored[:, 1] *= -1.0
        assert ring[60 + index] == pytest.approx(mirrored, abs=1e-12)


def test_a_single_branch_ring_runs_the_way_it_was_told_to():
    forward = manifold(branches=1, rotation="ccw").sweep(n_points=32)
    reverse = manifold(branches=1, rotation="cw").sweep(n_points=32)
    mirrored = reverse.copy()
    mirrored[..., 1] *= -1.0
    assert forward == pytest.approx(mirrored, abs=1e-12)


def test_a_split_ring_meets_itself_cleanly_at_the_port():
    """Both branches arrive there the same size, so the cavity closes."""
    ring = manifold(branches=2).sweep(n_points=32)
    assert ring[0] == pytest.approx(ring[-1], abs=1e-12)


def test_a_single_branch_ring_has_a_tongue_where_it_meets_itself():
    """It comes back round to the dead point at full size, so there is a step.

    That step is the cutwater, and it is hardware rather than a defect: a ring
    running one way round has to part the full passage from the empty one
    somewhere. A split ring has no such place, which is one of the reasons it
    is the usual choice.
    """
    ring = manifold(branches=1)
    swept = ring.sweep(n_points=32)
    assert ring.geometry.size[-1] > 2.0*ring.geometry.size[0]
    assert np.abs(swept[0] - swept[-1]).max() > 1.0e-3


def test_the_dead_point_sits_at_azimuth_zero():
    ring = manifold(branches=2).sweep(n_points=32)
    dead = ring[60].mean(axis=0)
    assert dead[1] == pytest.approx(0.0, abs=1e-12)
    assert dead[0] > 0.0


# -- drawing it more finely than it was resolved -------------------------


def test_the_sweep_says_how_many_stations_it_will_return():
    ring = manifold(branches=2)
    assert ring.sweep_stations() == len(ring.sweep(n_points=4)) == 2*61 - 1
    assert ring.sweep_stations(121) == len(ring.sweep(n_points=4, n_stations=121))


def test_a_finer_drawing_grid_agrees_with_the_law_where_the_law_was_resolved():
    """Refinement is drawing, so it may not move the passage it draws.

    Doubling the grid keeps every resolved station and puts one new station
    between each pair, so the sections the law sized are reproduced exactly and
    only the gaps are filled in.
    """
    ring = manifold(branches=2)
    fine = ring.display_geometry(2*61 - 1)
    assert fine.size[::2] == pytest.approx(ring.geometry.size, abs=1e-15)
    assert fine.area[::2] == pytest.approx(ring.geometry.area, abs=1e-15)
    assert fine.radial[::2] == pytest.approx(ring.geometry.radial, abs=1e-15)
    swept = ring.sweep(n_points=16, n_stations=2*61 - 1)
    assert swept[::2] == pytest.approx(ring.sweep(n_points=16), abs=1e-15)


def test_a_finer_drawing_grid_really_smooths_the_ring():
    """Azimuth is the direction a swept ring facets in, not its section.

    The ring is a whole revolution across and its section a few centimetres, so
    more stations is what a smooth picture needs -- and what it buys is flats
    shorter in proportion.
    """
    ring = manifold(branches=2)
    def facet(**kwargs):
        swept = ring.sweep(n_points=16, **kwargs)
        return np.linalg.norm(np.diff(swept, axis=0), axis=-1).max()
    assert facet(n_stations=181) == pytest.approx(facet()/3.0, rel=0.02)


def test_the_reserved_slit_is_drawn_on_the_same_grid_as_the_ring():
    """It tiles the chamber-facing face, so it has to be refined alongside it."""
    ring = manifold(branches=2)
    stations = ring.sweep_stations(121)
    assert len(ring.sweep_slit(n_points=6, n_stations=121)) == stations
    assert all(len(strip) == stations
               for strip in ring.sweep_margins(n_points=6, n_stations=121))
    assert len(ring.fan_origins(121)) == stations


def test_a_drawing_grid_of_two_stations_is_no_ring_at_all():
    with pytest.raises(ValueError, match="at least three"):
        manifold(branches=2).sweep(n_points=16, n_stations=2)


# -- the port ------------------------------------------------------------


def test_the_port_is_aimed_from_the_ellipse_centre_at_the_port_end():
    ring = manifold(branches=2, port=RadialPort(length=0.05))
    duct = ring.port_geometry
    depth, _, _ = ring.section.dome(ring.geometry.size[-1:],
                                    ring.geometry.aspect[-1:])
    # The stub is aimed from where the port station's section sits, and leaves
    # the dome its own depth out along the normal from there.
    assert duct.centre == pytest.approx([-float(ring.geometry.radial[-1]), 0.0,
                                         float(ring.geometry.axial[-1])],
                                        abs=1e-9)
    assert duct.start == pytest.approx(
        duct.centre + float(depth[0])*duct.start_direction)
    assert duct.arc_length == pytest.approx(0.05 - float(depth[0]))
    assert duct.is_straight


def test_the_stub_is_sized_from_the_dome_it_leaves():
    ring = manifold(branches=2, port=RadialPort(length=0.05))
    _, half_width, arc_limit = ring.section.dome(ring.geometry.size[-1:],
                                                 ring.geometry.aspect[-1:])
    # Not the whole ellipse axis: past the arc limit it is fillet, not passage.
    across = float(half_width[0])*math.sin(float(arc_limit[0]))
    assert ring.port_geometry.final_diameter == pytest.approx(0.95*2.0*across)


def test_the_swept_stub_is_a_cylinder_of_that_diameter():
    """Round, and square across the duct rather than flattened into it.

    A constant radius alone does not say the ring is a circle in the right
    plane: spanning the cross-section with a vector along the duct axis keeps
    every vertex one radius from the centre while collapsing the tube into the
    plane it is supposed to be normal to.
    """
    ring = manifold(branches=2, port=RadialPort(length=0.05))
    rings = ring.sweep_port(n_points=24, n_stations=5)
    assert rings.shape == (5, 24, 3)
    radius = ring.port_geometry.final_diameter/2.0
    axis = ring.port_geometry.start_direction
    # All but the first, which is trimmed onto the dome rather than laid flat.
    for station in rings[1:]:
        offsets = station - station.mean(axis=0)
        assert np.linalg.norm(offsets, axis=-1) == pytest.approx(radius, rel=1e-9)
        assert offsets @ axis == pytest.approx(np.zeros(len(offsets)), abs=1e-12)
    assert rings[-1].mean(axis=0) == pytest.approx(ring.port_geometry.end)


def test_the_drawn_stub_is_sunk_into_the_passage_rather_than_floating_on_it():
    """Every vertex of its first ring is inside the cavity at its own azimuth.

    The trim is taken against the port station's own ellipse, which is not the
    surface: the ring sweeps round a circle and the section shrinks away from
    the port, so an untrimmed tube shows daylight underneath. Erring inward is
    the whole point -- a tube buried a little in metal reads correctly, a gap
    does not.
    """
    ring = manifold(branches=2, port=RadialPort(length=0.05))
    rim = ring.sweep_port(n_points=48, n_stations=4)[0]
    azimuth = np.abs(np.arctan2(rim[:, 1], rim[:, 0]))
    geometry = ring.geometry
    sizes = np.interp(azimuth, geometry.theta, geometry.size)
    aspects = np.interp(azimuth, geometry.theta, geometry.aspect)
    polygon = ring.section.polygon(sizes, n_points=200, aspect=aspects)
    axial, radial, tilt = ring.section.meridional_frame(sizes, aspects)
    outward, u_hat, v_hat = section_frames(azimuth, tilt=tilt, sign=1.0)
    axis = np.zeros_like(outward)
    axis[:, 2] = 1.0
    local = rim - (radial[:, None]*outward + axial[:, None]*axis)
    for index, vertex in enumerate(np.column_stack((
            np.einsum("ij,ij->i", local, u_hat),
            np.einsum("ij,ij->i", local, v_hat)))):
        assert Path(polygon[index]).contains_point(vertex, radius=-1e-12)


def test_a_port_stated_against_the_manifold_plane_ignores_the_sections_lean():
    """On a sloping wall the two references differ by exactly that lean."""
    slope = math.radians(15.0)
    ring = manifold(branches=2, section=section(profile=cone(slope)),
                    port=RadialPort(length=0.05,
                                    angle_reference="manifold_plane"))
    assert ring.geometry.tilt[-1] == pytest.approx(-slope, abs=1e-6)
    duct = ring.port_geometry
    assert duct.angle == pytest.approx(slope, abs=1e-6)
    # Straight out from the machine axis, so it goes nowhere along it.
    assert duct.start_direction[2] == pytest.approx(0.0, abs=1e-9)
    assert duct.end[2] == pytest.approx(duct.centre[2], abs=1e-9)
    # Against the section instead, the same stub leans with the wall.
    leaning = manifold(branches=2, section=section(profile=cone(slope)),
                       port=RadialPort(length=0.05)).port_geometry
    assert leaning.start_direction[2] < 0.0


def test_an_angled_stub_is_narrower_than_one_straight_out():
    """It has turned toward the fillet, which leaves it less dome to sit on."""
    straight = manifold(branches=2, port=RadialPort(length=0.05))
    angled = manifold(branches=2,
                      port=RadialPort(length=0.05, angle=math.radians(30.0)))
    assert (angled.port_geometry.final_diameter
            < straight.port_geometry.final_diameter)
    assert angled.sweep_port(n_points=16, n_stations=4).shape == (4, 16, 3)


def test_a_manifold_may_be_fed_tangentially_instead():
    """The whole point of sharing the port protocol with the volute."""
    ring = manifold(branches=1, rotation="ccw",
                    port=ExitCone(0.05, first_angle="tangential"))
    assert ring.port_geometry.is_straight
    assert ring.port_geometry.final_diameter == pytest.approx(0.05)


def test_a_manifold_without_a_port_says_so():
    ring = manifold(branches=2)
    assert ring.port_geometry is None
    with pytest.raises(ValueError, match="no port"):
        ring.sweep_port()


def test_a_port_that_is_not_one_is_refused():
    with pytest.raises(TypeError, match="port must be a Port"):
        manifold(branches=2, port="a pipe")


# -- where the section sits ----------------------------------------------


def test_a_section_that_places_itself_needs_no_radius():
    ring = manifold(branches=2)
    assert ring.geometry.radial == pytest.approx(0.20, abs=1e-6)
    assert ring.geometry.tilt == pytest.approx(0.0, abs=1e-9)
    with pytest.raises(ValueError, match="places itself"):
        manifold(branches=2, radius=0.2)


def test_a_section_that_does_not_place_itself_needs_a_radius():
    plain = EllipticalSlottedSection(0.012, 0.003, 0.002, max_aspect_ratio=4.0)
    with pytest.raises(ValueError, match="needs a radius"):
        Manifold(plain, UniformArea(4.0e-4), branches=2, n_stations=21)
    ring = Manifold(plain, UniformArea(4.0e-4), branches=2, radius=0.20,
                    n_stations=21)
    assert ring.geometry.radial == pytest.approx(0.20)
    assert ring.sweep(n_points=16).shape[0] == 2*21 - 1


def test_a_family_with_no_dome_cannot_carry_a_radial_stub():
    """It offers a stub nowhere to sit, so the site carries no dome to ask for."""
    plain = EllipticalSlottedSection(0.012, 0.003, 0.002, max_aspect_ratio=4.0)
    with pytest.raises(ValueError, match="'dome'"):
        Manifold(plain, UniformArea(4.0e-4), branches=2, radius=0.20,
                 n_stations=21, port=RadialPort(length=0.04))


def test_a_family_with_no_dome_still_takes_a_tangential_duct():
    """Which is how the chamber headers carry an exit cone."""
    plain = EllipticalSlottedSection(0.012, 0.003, 0.002, max_aspect_ratio=4.0)
    ring = Manifold(plain, UniformArea(4.0e-4), branches=2, radius=0.20,
                    n_stations=21, port=ExitCone(0.05, first_angle="radial"))
    assert ring.port_geometry.final_diameter == pytest.approx(0.05)
    assert ring.port_geometry.turn_angle > 0.0


# -- the reserved slit ---------------------------------------------------


def test_the_slit_is_swept_inside_the_face_it_is_reserved_on():
    ring = manifold(branches=2)
    slit = ring.sweep_slit(n_points=12)
    face = ring.sweep(n_points=32)[:, :32]
    assert slit.shape[0] == face.shape[0]
    # Every slit vertex lies on the wall, so it sits within the face's reach.
    radius = np.linalg.norm(slit, axis=-1)
    assert radius.min() >= np.linalg.norm(face, axis=-1).min() - 1e-9
    assert radius.max() <= np.linalg.norm(face, axis=-1).max() + 1e-9


@pytest.mark.parametrize("placement, edge", [(0.0, "low"), (1.0, "high")])
def test_the_slit_sits_where_the_exposure_puts_it(placement, edge):
    """All the spare wall goes one way, so the slit ends up against the other."""
    sec = section(exposure_placement=placement)
    ring = Manifold(sec, tapered(), branches=2, n_stations=41)
    axial = ring.sweep_slit(n_points=12)[..., 2]
    face = ring.sweep(n_points=32)[:, :32, 2]
    if edge == "low":
        assert axial.min() == pytest.approx(face.min(), abs=1e-9)
        assert axial.max() < face.max() - 1e-3
    else:
        assert axial.max() == pytest.approx(face.max(), abs=1e-9)
        assert axial.min() > face.min() + 1e-3


def test_the_slit_keeps_its_height_however_the_section_is_sized():
    """It is the one thing the area law may not close over."""
    ring = manifold(branches=2, area_law=tapered())
    slit = ring.sweep_slit(n_points=2)
    span = np.linalg.norm(slit[:, 0] - slit[:, -1], axis=-1)
    assert span == pytest.approx(0.014, rel=1e-6)


def test_a_family_that_reserves_nothing_has_no_slit_to_sweep():
    plain = EllipticalSlottedSection(0.012, 0.003, 0.002, max_aspect_ratio=4.0)
    ring = Manifold(plain, UniformArea(4.0e-4), branches=2, radius=0.20,
                    n_stations=21)
    with pytest.raises(ValueError, match="reserves no slit"):
        ring.sweep_slit()


def test_a_section_says_how_big_it_can_get_and_the_law_is_told_when_it_cannot():
    sec = section(exposure_placement=0.5)
    low, high = sec.profile.span
    # Spare wall divides evenly, so whichever side runs out first caps it.
    assert sec.maximum_size == pytest.approx(
        sec.inlet_width + 2.0*min(sec.slit[0] - low, high - sec.slit[1]))
    with pytest.raises(RuntimeError, match="more area than this section can build"):
        Manifold(sec, ConstantVelocity(1.0), branches=2, n_stations=9)


# -- what the scene shows ------------------------------------------------


def component(ring, results=None, **kwargs):
    import pyskyfire as psf
    kwargs.setdefault("n_points", 12)
    return psf.viz.ManifoldComponent(ring, results, **kwargs).build()


def test_only_the_faces_flow_crosses_are_painted():
    """Walls carry fields; openings carry the colour that says they are holes.

    Painting the whole chamber-facing face would claim the manifold is open all
    along the wall, when only the slit is.
    """
    node = component(manifold(branches=2, port=RadialPort(length=0.05)))
    by_name = {child.name: child for child in node.children}
    assert set(by_name) == {"passage", "slit", "outlet"}
    assert by_name["passage"].color != "red"
    assert by_name["slit"].color == "red"
    assert by_name["outlet"].color == "red"
    # The stub's tube is wall, so it is far larger than the flange closing it.
    assert len(by_name["passage"].mesh.points) > 10*len(by_name["outlet"].mesh.points)


def test_the_slit_and_the_wall_around_it_do_not_overlap():
    """They tile the face, so neither is drawn on top of the other."""
    ring = manifold(branches=2)
    above, below = ring.sweep_margins(n_points=6)
    slit = ring.sweep_slit(n_points=6)
    assert above[:, -1] == pytest.approx(slit[:, 0], abs=1e-12)
    assert slit[:, -1] == pytest.approx(below[:, 0], abs=1e-12)
    face = ring.sweep(n_points=32)[:, :32]
    assert above[:, 0] == pytest.approx(face[:, 0], abs=1e-12)
    assert below[:, -1] == pytest.approx(face[:, -1], abs=1e-12)


def test_fields_cover_every_wall_vertex_once():
    ring = manifold(branches=2, port=RadialPort(length=0.05))
    count = len(ring.sweep(n_points=12))
    values = np.arange(count, dtype=float)
    import pyskyfire as psf
    node = component(ring, {"v": psf.viz.FieldData(values, "m/s", "v")})
    passage = node.children[0]
    assert len(passage.fields["v"].values) == len(passage.mesh.points)


def test_fields_stay_on_the_grid_they_were_solved_on_when_the_mesh_is_refined():
    """A finer picture is not a finer solution.

    Values belong to the stations the ring was resolved on, so a refined mesh
    reads between them rather than demanding one value per drawn vertex. A
    lumped result is constant around the ring and stays constant, which is the
    case that matters here: it says the network resolved one state, not many.
    """
    import pyskyfire as psf
    ring = manifold(branches=2, port=RadialPort(length=0.05))
    resolved = ring.sweep_stations()
    lumped = {"v": psf.viz.FieldData(np.full(resolved, 7.0), "m/s", "v")}
    passage = component(ring, lumped, n_stations=121).children[0]
    assert len(passage.fields["v"].values) == len(passage.mesh.points)
    assert passage.fields["v"].values == pytest.approx(7.0)
    # One value per drawn station is the wrong shape, however plausible it looks.
    with pytest.raises(ValueError, match="resolved ring"):
        component(ring, {"v": psf.viz.FieldData(
            np.zeros(ring.sweep_stations(121)), "m/s", "v")}, n_stations=121)


def test_a_section_with_a_real_hole_paints_its_whole_face():
    """A slotted section is open all along its mouth, unlike a walled one."""
    plain = EllipticalSlottedSection(0.012, 0.003, 0.002, max_aspect_ratio=4.0)
    ring = Manifold(plain, UniformArea(4.0e-4), branches=2, radius=0.20,
                    n_stations=21)
    names = {child.name: child.color for child in component(ring).children}
    assert names["opening"] == "red"
    assert "slit" not in names
