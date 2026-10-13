"""The shared port protocol, and the straight radial stub a manifold uses.

Both kinds of outlet resolve to the same thing, so a passage can offer either
without knowing which it was handed. That is what stops a manifold fed
tangentially from costing a second implementation of the same duct.
"""

import math

import numpy as np
import pytest

from pyskyfire.ducts import (DomeSite, ExitCone, Port, PortSite, RadialPort,
                             ResolvedPort, ResolvedRadialPort)
from pyskyfire.ducts.exit_cone import ResolvedExitCone

#: Where the dome sits, and how far round it is passage rather than fillet.
SEMI_NORMAL, SEMI_IN_PLANE, ARC_LIMIT = 0.018, 0.022, math.radians(80.0)


def site(**kwargs):
    kwargs.setdefault("start", np.array([0.070, 0.0, 0.0]))
    kwargs.setdefault("outlet_diameter", 0.030)
    kwargs.setdefault("max_section_radius", 0.020)
    return PortSite(**kwargs)


def dome(**kwargs):
    """A section facing outward and up, as a sloping chamber wall leaves it."""
    kwargs.setdefault("normal", np.array([0.0, 0.6, 0.8]))
    kwargs.setdefault("in_plane", np.array([0.0, -0.8, 0.6]))
    kwargs.setdefault("semi_normal", SEMI_NORMAL)
    kwargs.setdefault("semi_in_plane", SEMI_IN_PLANE)
    kwargs.setdefault("arc_limit", ARC_LIMIT)
    return DomeSite(**kwargs)


def manifold_site(**kwargs):
    kwargs.setdefault("dome", dome())
    return site(**kwargs)


# -- the protocol --------------------------------------------------------


def test_both_outlets_are_ports_and_resolve_to_the_same_shape():
    cone = ExitCone(0.05, first_angle="radial", arc_length_parameter=3.0)
    stub = RadialPort(length=0.06)
    assert isinstance(cone, Port) and isinstance(stub, Port)

    resolved = (cone.resolve(site()), stub.resolve(manifold_site()))
    for duct in resolved:
        assert isinstance(duct, ResolvedPort)
        for name in ("start", "end", "start_direction", "end_direction",
                     "arc_length", "turn_angle", "bend_radius",
                     "outlet_diameter", "final_diameter"):
            assert hasattr(duct, name)
        points, u, v = duct.frames(9)
        assert points.shape == (9, 3) and u.shape == (9, 3) and v.shape == (9, 3)
        assert duct.sweep(n_points=16, n_stations=9).shape == (9, 16, 3)


def test_a_cone_keeps_its_own_angles_on_top_of_the_shared_ones():
    cone = ExitCone(0.05, first_angle="radial", arc_length_parameter=3.0)
    duct = cone.resolve(site())
    assert isinstance(duct, ResolvedExitCone)
    assert duct.first_angle == pytest.approx(duct.radial_angle)
    assert duct.end_azimuth > 0.0


def test_a_site_says_which_field_a_port_wanted():
    with pytest.raises(ValueError, match="'dome'"):
        RadialPort(length=0.05).resolve(site())


# -- the radial stub -----------------------------------------------------


def test_the_stub_leaves_along_the_way_the_section_faces():
    where = manifold_site()
    duct = RadialPort(length=0.06).resolve(where)
    assert isinstance(duct, ResolvedRadialPort)
    assert duct.start_direction == pytest.approx(where.dome.normal)
    assert duct.end_direction == pytest.approx(where.dome.normal)
    assert duct.end == pytest.approx(where.start + 0.06*where.dome.normal)


def test_the_length_is_measured_from_the_ellipse_centre():
    """A defined point, unlike the curved surface, so the flange stays put."""
    where = manifold_site()
    duct = RadialPort(length=0.06).resolve(where)
    assert duct.centre == pytest.approx(where.start)
    # Straight out, the dome is its own depth away from that centre.
    assert duct.start == pytest.approx(
        where.start + SEMI_NORMAL*where.dome.normal)
    assert duct.arc_length == pytest.approx(0.06 - SEMI_NORMAL)


def test_a_stub_too_short_to_clear_the_dome_is_refused():
    with pytest.raises(ValueError, match="does not clear the dome"):
        RadialPort(length=0.5*SEMI_NORMAL).resolve(manifold_site())


def test_a_stub_is_steered_by_its_angle_within_the_section_plane():
    angle = math.radians(25.0)
    where = manifold_site()
    duct = RadialPort(length=0.06, angle=angle).resolve(where)
    expected = (math.cos(angle)*where.dome.normal
                + math.sin(angle)*where.dome.in_plane)
    assert duct.start_direction == pytest.approx(expected)
    assert duct.angle == pytest.approx(angle)
    # Aimed from the same centre, so the flange is still one length out.
    assert np.linalg.norm(duct.end - where.start) == pytest.approx(0.06)


def test_the_plane_reference_turns_the_stub_through_the_sections_own_lean():
    """Stated against the machine, the same number is the same hardware."""
    lean = math.radians(12.0)
    where = manifold_site(dome=dome(tilt=lean))
    stated = math.radians(30.0)
    duct = RadialPort(length=0.06, angle=stated,
                      angle_reference="manifold_plane").resolve(where)
    assert duct.angle == pytest.approx(stated - lean)
    same = RadialPort(length=0.06, angle=stated - lean).resolve(where)
    assert duct.start_direction == pytest.approx(same.start_direction)


def test_a_reference_that_is_not_one_is_refused():
    with pytest.raises(ValueError, match="angle_reference must be"):
        RadialPort(length=0.06, angle_reference="outward")


def test_the_stub_runs_straight_and_is_charged_no_turn():
    duct = RadialPort(length=0.06).resolve(manifold_site())
    assert duct.is_straight
    assert duct.turn_angle == pytest.approx(0.0)
    assert not math.isfinite(duct.bend_radius)


def test_a_straight_stub_carries_its_frame_without_twisting():
    duct = RadialPort(length=0.06).resolve(manifold_site())
    points, u, v = duct.frames(11)
    assert np.allclose(u, u[0]) and np.allclose(v, v[0])
    # The frame spans the plane the stub's own cross-section lives in.
    for vector in (u[0], v[0]):
        assert np.dot(vector, duct.start_direction) == pytest.approx(0.0, abs=1e-12)
        assert np.linalg.norm(vector) == pytest.approx(1.0)
    assert np.dot(u[0], v[0]) == pytest.approx(0.0, abs=1e-12)
    assert points[-1] == pytest.approx(duct.end)


# -- how wide the dome lets it be ----------------------------------------


def test_a_stub_is_sized_by_how_far_the_dome_reaches_across_it():
    """Not the ellipse's own axis: the part of it past the fillets is not passage."""
    across = manifold_site().dome.clearance(0.0)
    assert across == pytest.approx(SEMI_IN_PLANE*math.sin(ARC_LIMIT))
    assert across < SEMI_IN_PLANE
    duct = RadialPort(length=0.06).resolve(manifold_site())
    assert duct.final_diameter == pytest.approx(0.95*2.0*across)


@pytest.mark.parametrize("fraction", [0.4, 0.75, 1.0])
def test_a_fractional_stub_is_sized_against_that_room(fraction):
    across = manifold_site().dome.clearance(0.0)
    duct = RadialPort(length=0.06, diameter_fraction=fraction).resolve(
        manifold_site())
    assert duct.final_diameter == pytest.approx(fraction*2.0*across)


def test_the_room_closes_down_as_the_stub_turns_toward_the_fillets():
    where = manifold_site()
    across = [where.dome.clearance(math.radians(deg)) for deg in (0, 15, 30, 45)]
    assert np.all(np.diff(across) < 0.0)
    # Symmetric, since the dome is.
    assert (where.dome.clearance(math.radians(-30.0))
            == pytest.approx(where.dome.clearance(math.radians(30.0))))


def test_a_stated_diameter_is_taken_as_given():
    duct = RadialPort(length=0.06, diameter=0.030).resolve(manifold_site())
    assert duct.final_diameter == pytest.approx(0.030)


def test_a_stub_wider_than_its_dome_is_refused():
    """It would leave through the fillet rather than the dome."""
    across = manifold_site().dome.clearance(0.0)
    with pytest.raises(ValueError, match="does not fit"):
        RadialPort(length=0.06, diameter=2.2*across).resolve(manifold_site())


def test_a_stub_exactly_as_wide_as_its_dome_is_allowed():
    across = manifold_site().dome.clearance(0.0)
    duct = RadialPort(length=0.06, diameter=2.0*across).resolve(manifold_site())
    assert duct.final_diameter == pytest.approx(2.0*across)


def test_a_stub_aimed_past_the_dome_is_refused():
    """Past the arc limit the ellipse is outboard of the fillet, not passage."""
    limit = manifold_site().dome.angle_limit
    assert limit == pytest.approx(math.atan2(SEMI_IN_PLANE*math.sin(ARC_LIMIT),
                                            SEMI_NORMAL*math.cos(ARC_LIMIT)))
    with pytest.raises(ValueError, match="through the fillet"):
        RadialPort(length=0.06, angle=limit + 1.0e-3).resolve(manifold_site())
    # Aimed at the very end of the dome there is no room left across it.
    with pytest.raises(ValueError, match="no room across"):
        RadialPort(length=0.06, angle=limit).resolve(manifold_site())
    assert RadialPort(length=0.06, angle=0.98*limit).resolve(
        manifold_site()).final_diameter > 0.0


# -- what gets drawn -----------------------------------------------------


def test_the_drawn_stub_is_trimmed_back_onto_the_dome():
    """Its first ring follows the dome rather than floating a chord above it."""
    duct = RadialPort(length=0.06).resolve(manifold_site())
    rings = duct.sweep(n_points=32, n_stations=6)
    assert rings.shape == (6, 32, 3)
    along = (rings[0] - duct.start) @ duct.start_direction
    # Sunk below the outlet, and not by the same amount all the way round.
    assert np.all(along < 0.0)
    assert np.ptp(along) > 0.0
    # The rest of the tube is plain circles, the last of them on the flange.
    for station in rings[1:]:
        offsets = station - station.mean(axis=0)
        assert np.linalg.norm(offsets, axis=-1) == pytest.approx(
            duct.final_diameter/2.0)
        assert offsets @ duct.start_direction == pytest.approx(
            np.zeros(len(offsets)), abs=1e-15)
    assert rings[-1].mean(axis=0) == pytest.approx(duct.end)
    # Nothing in the trim band overtakes the ring that closes it.
    assert np.all(along.max() <= (rings[1] - duct.start) @ duct.start_direction)


def test_the_trim_follows_the_dome_across_the_section_and_not_around_it():
    """The passage carries straight on azimuthally, so only one way is trimmed."""
    duct = RadialPort(length=0.06).resolve(manifold_site())
    across, azimuthal = duct._start_frame()
    rings = duct.sweep(n_points=4, n_stations=4)
    along = (rings[0] - duct.start) @ duct.start_direction
    # Vertices at +/- the same offset across the section are sunk alike; the
    # two out of it are the shallowest, being nearest the crest.
    assert along[0] == pytest.approx(along[2])
    assert along[1] == pytest.approx(along[3])
    assert along[1] > along[0]
    assert np.dot(across, duct.dome.normal) == pytest.approx(0.0, abs=1e-15)
    assert np.dot(azimuthal, duct.dome.in_plane) == pytest.approx(0.0, abs=1e-15)


@pytest.mark.parametrize("bad", [0.0, -0.1])
def test_a_stub_needs_a_real_length(bad):
    with pytest.raises(ValueError, match="length must be finite and positive"):
        RadialPort(length=bad)


@pytest.mark.parametrize("bad", [0.0, -0.5, 1.5])
def test_a_fraction_outside_its_range_is_refused(bad):
    with pytest.raises(ValueError, match="diameter_fraction must lie"):
        RadialPort(length=0.06, diameter_fraction=bad)


def test_a_stated_diameter_must_be_positive():
    with pytest.raises(ValueError, match="diameter must be finite and positive"):
        RadialPort(length=0.06, diameter=-0.01)


# -- what the volute still expects --------------------------------------


def test_a_volute_takes_any_port_and_refuses_anything_else():
    from pyskyfire import ducts

    section = ducts.EllipticalSlottedSection(0.012, 0.003, 0.002)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    law = ducts.ConstantAngularMomentum()
    volute = ducts.Volute(section, inlet, law, rotation="ccw", n_stations=9,
                          exit_cone=ExitCone(0.05, first_angle="tangential"))
    assert isinstance(volute.exit_geometry, ResolvedPort)
    with pytest.raises(TypeError, match="must be a Port"):
        ducts.Volute(section, inlet, law, rotation="ccw", n_stations=9,
                     exit_cone="a pipe")


def test_the_cone_still_reverses_with_the_volute_it_leaves():
    """Tangential is the direction the flow goes, so the sense still matters."""
    cone = ExitCone(0.05, first_angle="radial", arc_length_parameter=3.0)
    forward = cone.resolve(site(sign=1.0))
    reversed_ = cone.resolve(site(sign=-1.0))
    assert reversed_.start_direction == pytest.approx(-forward.start_direction)
    assert reversed_.turn_angle == pytest.approx(forward.turn_angle)
    assert reversed_.end == pytest.approx(forward.end*[1.0, -1.0, 1.0])
