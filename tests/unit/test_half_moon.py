"""The half-moon manifold section.

The inner boundary has to be the chamber wall itself rather than an
approximation of it, the slit has to survive every size the area law asks for,
and the fillets have to actually touch the ellipse they are supposed to blend
into.
"""

import math

import numpy as np
import pytest

from pyskyfire.ducts import HalfMoonSection, InnerProfile


def flat_wall(radius=0.20, reach=0.09, samples=401):
    """A cylindrical wall, where chord and arc coincide."""
    s = np.linspace(-reach, reach, samples)
    return InnerProfile(s, np.column_stack((s, np.full_like(s, radius))))


def curved_wall(curvature_radius, radius=0.20, reach=0.09, samples=401):
    """A wall bent to a known radius of curvature in the meridional plane."""
    s = np.linspace(-reach, reach, samples)
    angle = s/curvature_radius
    return InnerProfile(s, np.column_stack((
        curvature_radius*np.sin(angle),
        radius + curvature_radius*(1.0 - np.cos(angle)))))


def section(profile=None, slit_height=0.014, fillet_radius=0.003, **kwargs):
    kwargs.setdefault("aspect_ratio", 1.2)
    kwargs.setdefault("max_aspect_ratio", 6.0)
    return HalfMoonSection(profile if profile is not None else flat_wall(),
                           slit_height, fillet_radius, **kwargs)


def shoelace(ring):
    return 0.5*abs(np.sum(ring[:-1, 0]*np.diff(ring[:, 1])
                          - ring[:-1, 1]*np.diff(ring[:, 0])))


# -- the shape itself ---------------------------------------------------


def test_the_boundary_closes_and_encloses_the_area_it_reports():
    sec = section()
    ring = sec.polygon(0.040, n_points=6000)[0]
    assert ring[0] == pytest.approx(ring[-1], abs=1e-15)
    assert float(sec.metrics(0.040)[0][0]) == pytest.approx(shoelace(ring), rel=1e-6)


@pytest.mark.parametrize("aspect", [0.6, 1.0, 1.8, 4.0])
def test_every_join_is_continuous_in_position_and_slope(aspect):
    """A flat wall leaves the whole boundary smooth, fillets included."""
    sec = section(aspect_ratio=aspect, max_aspect_ratio=6.0)
    pieces = sec.boundary(0.040)
    step = 1.0e-7
    for index, piece in enumerate(pieces):
        following = pieces[(index + 1) % len(pieces)]
        assert piece.evaluate(np.array([1.0]))[0, 0] == pytest.approx(
            following.evaluate(np.array([0.0]))[0, 0], abs=1e-12)
        before = piece.derivative(np.array([1.0 - step]))[0, 0]
        after = following.derivative(np.array([step]))[0, 0]
        cosine = (before @ after)/np.linalg.norm(before)/np.linalg.norm(after)
        assert math.degrees(math.acos(min(1.0, cosine))) < 1.0e-3


@pytest.mark.parametrize("aspect", [0.6, 1.0, 2.5])
def test_each_fillet_really_touches_the_ellipse(aspect):
    """Checked by distance to the ellipse, not by the solve that placed it."""
    sec = section(aspect_ratio=aspect, max_aspect_ratio=6.0)
    size, radius = 0.040, sec.fillet_radius
    *_, semi_u, semi_v, _, _ = sec._layout(size)
    parameters = np.linspace(-np.pi, np.pi, 400001)
    ellipse = np.column_stack((float(semi_u[0])*np.cos(parameters),
                               float(semi_v[0])*np.sin(parameters)))
    for centre in sec.boundary(size)[1].centers, sec.boundary(size)[3].centers:
        gap = np.hypot(*(ellipse - centre[0]).T).min()
        assert gap == pytest.approx(radius, rel=1e-5)


def test_a_zero_fillet_leaves_the_dome_sitting_on_the_chord():
    """With nothing to blend, the ellipse spans the exposed wall exactly."""
    sec = section(fillet_radius=0.0, aspect_ratio=1.0)
    size = 0.040
    *_, chord, semi_u, semi_v, _, _ = sec._layout(size)
    assert float(semi_v[0]) == pytest.approx(float(chord[0])/2.0, rel=1e-9)
    assert float(semi_u[0]) == pytest.approx(float(chord[0])/2.0, rel=1e-9)
    # Half an ellipse on a flat wall, so half the ellipse's area.
    assert float(sec.metrics(size)[0][0]) == pytest.approx(
        0.5*np.pi*float(semi_u[0])*float(semi_v[0]), rel=1e-9)


@pytest.mark.parametrize("radius", [0.0, 0.002, 0.005, 0.009])
def test_a_round_corner_needs_a_fuller_dome_to_stay_tangent(radius):
    """Checked against the closed form a circular dome admits.

    ``size`` is the wall the section actually touches, so rounding the corner
    without giving up any of it means the dome has to grow to meet the fillet
    tangentially. For a circular dome the radius that does it falls out in
    closed form, which pins the whole tangency solve to something independent.
    """
    sec = section(fillet_radius=radius, aspect_ratio=1.0)
    size = 0.040
    *_, semi_u, semi_v, _, _ = sec._layout(size)
    expected = radius + math.hypot(radius, size/2.0)
    assert float(semi_v[0]) == pytest.approx(expected, rel=1e-9)
    assert float(semi_u[0]) == pytest.approx(expected, rel=1e-9)


def test_a_rounded_corner_never_overhangs_its_wall_by_more_than_the_fillet():
    """The corner rolls off the wall, it does not flare out past it."""
    radius, size = 0.005, 0.040
    sec = section(fillet_radius=radius, aspect_ratio=1.2)
    ring = sec.polygon(size, n_points=4000)[0]
    reach = np.abs(ring[:, 1]).max()
    assert size/2.0 < reach <= size/2.0 + radius


# -- the wall it sits on ------------------------------------------------


def test_the_inner_boundary_is_the_wall_itself():
    """Not a chord or a spline fitted to it: the wall, sampled where asked."""
    profile = curved_wall(0.5)
    sec = section(profile)
    size = 0.040
    lower, upper = sec.window(size)
    *_, origin, axes, _ = sec._frame(np.atleast_1d(size))

    wall = sec.boundary(size)[0]
    fractions = np.linspace(0.0, 1.0, 41)
    drawn = wall.evaluate(fractions)[0]
    arc = float(upper[0]) + fractions*(float(lower[0]) - float(upper[0]))
    expected = np.einsum("ij,kj->ki", axes[0], profile.at(arc) - origin[0])
    assert drawn == pytest.approx(expected, abs=1e-15)


def test_a_curved_wall_meets_its_fillets_within_its_own_turn():
    """The fillets ride the chord, so the kink is the wall's turn, no more.

    The inner boundary still follows the wall exactly; what the chord costs is
    a slope mismatch where the wall hands over to the fillet, and it is bounded
    by how far the wall turns across half a section.
    """
    size = 0.040
    for curvature_radius in (2.0, 0.8, 0.4):
        sec = section(curved_wall(curvature_radius))
        pieces = sec.boundary(size)
        step = 1.0e-7
        for index in (0, 3):
            before = pieces[index].derivative(np.array([1.0 - step]))[0, 0]
            after = pieces[(index + 1) % 4].derivative(np.array([step]))[0, 0]
            cosine = (before @ after)/np.linalg.norm(before)/np.linalg.norm(after)
            kink = math.acos(min(1.0, cosine))
            assert kink == pytest.approx(size/2.0/curvature_radius, rel=1e-3)


def test_width_is_measured_along_the_wall_not_along_the_axis():
    """A sloped wall offers more surface than its station span suggests.

    Reserving the slit axially would hand back the difference, which is the
    room the channel breakout needs.
    """
    slope = 0.25
    s = np.linspace(-0.09, 0.09, 401)
    scale = math.hypot(1.0, slope)
    profile = InnerProfile(s, np.column_stack((s/scale, 0.2 + slope*s/scale)))
    sec = section(profile, slit_height=0.014)
    lower, upper = sec.window(0.014)
    axial = float(profile.at(upper)[0, 0] - profile.at(lower)[0, 0])
    assert float(upper[0] - lower[0]) == pytest.approx(0.014, rel=1e-9)
    assert axial == pytest.approx(0.014/scale, rel=1e-6)


# -- placement and sizing -----------------------------------------------


def test_the_section_closes_onto_the_slit_and_never_over_it():
    sec = section(slit_height=0.014)
    assert sec.minimum_size == pytest.approx(0.014)
    lower, upper = sec.window(sec.minimum_size)
    assert (float(lower[0]), float(upper[0])) == pytest.approx(sec.slit)
    with pytest.raises(ValueError, match="at least the slit height"):
        sec.window(0.013)


@pytest.mark.parametrize("placement", [0.0, 0.25, 0.5, 1.0])
def test_slit_placement_slides_the_slit_about_the_anchor(placement):
    sec = section(slit_height=0.014, slit_placement=placement)
    low, high = sec.slit
    assert high - low == pytest.approx(0.014)
    assert -low == pytest.approx(placement*0.014)


@pytest.mark.parametrize("placement", [0.0, 0.25, 0.75, 1.0])
def test_exposure_placement_divides_the_spare_wall_and_holds_it(placement):
    """The ratio is the same at every size, so it holds at every azimuth."""
    sec = section(slit_height=0.014, exposure_placement=placement)
    low, high = sec.slit
    for size in (0.020, 0.040, 0.070):
        lower, upper = sec.window(size)
        spare = size - 0.014
        assert low - float(lower[0]) == pytest.approx(placement*spare)
        assert float(upper[0]) - high == pytest.approx((1.0 - placement)*spare)


def test_the_slit_stays_put_however_the_section_is_sized():
    """Channels break out where they break out; the cavity moves around them."""
    sec = section(slit_height=0.014, exposure_placement=0.75)
    low, high = sec.slit
    for size in (0.014, 0.025, 0.040, 0.080):
        lower, upper = sec.window(size)
        assert float(lower[0]) <= low + 1e-15
        assert float(upper[0]) >= high - 1e-15


def test_the_exposed_wall_is_the_size_that_was_asked_for():
    sec = section()
    for size in (0.014, 0.030, 0.060):
        assert float(sec.exposed_wall(size)[0]) == pytest.approx(size, rel=1e-12)


def test_area_rises_with_width_then_falls_as_the_dome_flattens():
    """The two knobs the area law works through, in the order it uses them."""
    sec = section(aspect_ratio=1.0, max_aspect_ratio=6.0)
    sizes = np.linspace(0.014, 0.070, 30)
    areas = np.array([float(sec.metrics(s)[0][0]) for s in sizes])
    assert np.all(np.diff(areas) > 0.0)

    aspects = np.linspace(1.0, 6.0, 30)
    flattened = np.array([float(sec.metrics(0.030, a)[0][0]) for a in aspects])
    assert np.all(np.diff(flattened) < 0.0)
    # Flattening has to reach below the smallest the width alone can give, or
    # it would buy the area law nothing.
    assert flattened[-1] < float(sec.metrics(0.030)[0][0])


def test_the_fan_origin_sees_the_whole_boundary():
    """Including a wall that curves into the cavity and would hide itself."""
    for profile in (flat_wall(), curved_wall(0.4), curved_wall(-0.4)):
        sec = section(profile)
        origin = sec.fan_origin(0.040)[0]
        ring = sec.polygon(0.040, n_points=1200)[0]
        angles = np.unwrap(np.arctan2(ring[:, 1] - origin[1], ring[:, 0] - origin[0]))
        assert abs(abs(angles[-1] - angles[0]) - 2.0*np.pi) < 1e-6


# -- refusals ------------------------------------------------------------


def test_a_section_reaching_past_its_wall_says_so():
    sec = section(flat_wall(reach=0.03))
    with pytest.raises(ValueError, match="reaches past the wall"):
        sec.window(0.075)


def test_a_profile_must_straddle_its_anchor():
    s = np.linspace(0.01, 0.09, 40)
    with pytest.raises(ValueError, match="either side of its anchor"):
        InnerProfile(s, np.column_stack((s, np.full_like(s, 0.2))))


@pytest.mark.parametrize("name", ["slit_placement", "exposure_placement"])
def test_placements_outside_the_unit_interval_are_refused(name):
    with pytest.raises(ValueError, match=f"{name} must lie between zero and one"):
        section(**{name: 1.5})


def test_a_negative_fillet_is_refused():
    with pytest.raises(ValueError, match="fillet_radius must be finite"):
        section(fillet_radius=-0.001)


def test_the_profile_has_to_be_one():
    with pytest.raises(TypeError, match="profile must be an InnerProfile"):
        HalfMoonSection(object(), 0.014, 0.003)
