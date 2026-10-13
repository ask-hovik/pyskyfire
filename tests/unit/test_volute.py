"""Geometry, sizing-law and loss-model contracts for volutes and collectors."""

import math
import warnings

import numpy as np
import pytest

from pyskyfire import ducts
from pyskyfire.thermo import FluidState, PHASE_UNKNOWN

#: The shaped volute, which most of these tests exercise.
SHAPED = ducts.ConstantAngularMomentum()


class ConstantTransport:
    def __init__(self, rho=1000.0, mu=1.0e-3):
        self._rho = rho
        self._mu = mu

    def get_rho(self, T, p):
        return self._rho

    def get_mu(self, T, p):
        return self._mu


def inlet_state(p=8.0e6, T=300.0):
    return FluidState(
        T=T,
        p=p,
        h=np.nan,
        quality=np.nan,
        rho=np.nan,
        phase=PHASE_UNKNOWN,
        T_sat=np.nan,
    )


def pump_volute(SP=1.0, area_law=None):
    """A pump-discharge volute at turbopump scale, with its design mass flow."""
    section = ducts.EllipticalSlottedSection(
        inlet_width=0.012, fillet_radius=0.003, throat_length=0.002
    )
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    if area_law is None:
        area_law = ducts.ConstantAngularMomentum(SP=SP)
    volute = ducts.Volute(section, inlet, area_law, rotation="ccw")
    mass_flow = (
        1000.0 * inlet.C_m * 2.0 * math.pi * inlet.radius * section.inlet_width
    )
    return volute, mass_flow


def solve(volute, mass_flow, **kwargs):
    kwargs.setdefault("exit_diameter", 0.05)
    return volute.solve(
        inlet_state=inlet_state(),
        mass_flow=mass_flow,
        transport=ConstantTransport(),
        **kwargs,
    )


# -- cross-section geometry -------------------------------------------------


def test_sharp_cornered_section_matches_the_analytic_circle_minus_segment():
    """With no fillet the section is a bore less the segment its slot cuts."""
    width, bore = 0.02, 0.05
    section = ducts.EllipticalSlottedSection(inlet_width=width, fillet_radius=0.0)
    psi = math.asin(width / (2.0 * bore))
    expected = bore**2 * (math.pi - psi + math.sin(psi) * math.cos(psi))
    assert section.area(bore)[0] == pytest.approx(expected, rel=1e-12)


def test_throat_adds_exactly_its_rectangle_of_area():
    width, bore, throat = 0.02, 0.05, 0.004
    without = ducts.EllipticalSlottedSection(width, 0.0).area(bore)[0]
    with_throat = ducts.EllipticalSlottedSection(width, 0.0, throat).area(bore)[0]
    assert with_throat - without == pytest.approx(width * throat, rel=1e-12)


@pytest.mark.parametrize("fillet_radius", [0.0, 0.001, 0.005, 0.009])
def test_green_theorem_metrics_agree_with_a_dense_polygon(fillet_radius):
    """The exact-arc integral matches a finely sampled shoelace integration."""
    bore = 0.05
    section = ducts.EllipticalSlottedSection(0.02, fillet_radius, throat_length=0.002)
    ring = section.polygon(bore, n_points=20000)[0]
    u, v = ring[:, 0], ring[:, 1]
    cross = u[:-1] * v[1:] - u[1:] * v[:-1]
    area = 0.5 * cross.sum()
    centroid_u = np.sum((u[:-1] + u[1:]) * cross) / (6.0 * area)
    assert section.area(bore)[0] == pytest.approx(area, rel=1e-6)
    assert section.centroid_offset(bore)[0] == pytest.approx(centroid_u, rel=1e-6)


def test_section_is_symmetric_about_the_slot_centreline():
    section = ducts.EllipticalSlottedSection(0.02, 0.004, throat_length=0.001)
    assert section.metrics(0.05)[1][0, 1] == pytest.approx(0.0, abs=1e-15)


def test_bore_narrower_than_its_slot_is_rejected():
    section = ducts.EllipticalSlottedSection(inlet_width=0.02, fillet_radius=0.002)
    assert section.minimum_size == pytest.approx(0.01)
    with pytest.raises(ValueError, match="wider than the slot"):
        section.area(0.009)


def test_area_and_centroid_grow_with_bore_size():
    section = ducts.EllipticalSlottedSection(0.02, 0.004)
    sizes = np.linspace(0.011, 0.08, 40)
    assert np.all(np.diff(section.area(sizes)) > 0.0)
    assert np.all(np.diff(section.centroid_offset(sizes)) > 0.0)


def test_wetted_perimeter_excludes_the_slot_mouth():
    section = ducts.EllipticalSlottedSection(0.02, 0.004)
    full = section.perimeter(0.05, wetted=False)[0]
    wetted = section.perimeter(0.05, wetted=True)[0]
    assert full - wetted == pytest.approx(section.inlet_width, rel=1e-12)


def test_a_throat_pushes_the_bore_clear_of_the_inlet_plane():
    """Overhang is a reported diagnostic, and lengthening the throat cures it."""
    flush = ducts.EllipticalSlottedSection(0.012, 0.003, throat_length=0.0)
    assert flush.inboard_overhang(0.008)[0] > 0.0
    deep = ducts.EllipticalSlottedSection(0.012, 0.003, throat_length=0.01)
    assert deep.inboard_overhang(0.008)[0] == 0.0


# -- sizing law -------------------------------------------------------------


def test_constant_angular_momentum_makes_area_over_radius_linear_in_azimuth():
    """Aungier Eq. 11-5: A / r_c grows at exactly b5 * tan(alpha5) per radian."""
    volute, _ = pump_volute()
    geometry = volute.geometry
    unclipped = ~geometry.clipped
    ratio = geometry.area[unclipped] / geometry.centroid_radius[unclipped]
    slope = np.diff(ratio) / np.diff(geometry.theta[unclipped])
    expected = volute.section.inlet_width * volute.inlet.tan_alpha
    assert np.allclose(slope, expected, rtol=1e-6)


def test_sizing_parameter_scales_the_passage_against_its_own_centroid():
    """Eq. 11-5 scales ``A / r_c``, not ``A``.

    Oversizing pushes the section's centroid outward as well, so the area grows
    slightly faster than SP does. It is the ratio that is exactly linear.
    """
    base, _ = pump_volute(SP=1.0)
    scaled, _ = pump_volute(SP=1.2)
    ratio = lambda v: v.geometry.area[-1] / v.geometry.centroid_radius[-1]  # noqa: E731
    assert ratio(scaled) == pytest.approx(1.2 * ratio(base), rel=1e-6)
    assert scaled.geometry.area[-1] > 1.2 * base.geometry.area[-1]


def test_the_law_is_clipped_near_the_tongue_and_satisfied_thereafter():
    """A shaped volute wants zero area at the tongue; the section has a floor."""
    volute, _ = pump_volute()
    geometry = volute.geometry
    assert geometry.clipped[0]
    assert not geometry.clipped[-1]
    assert 0.0 < geometry.tongue_clip_angle < 2.0 * math.pi
    clipped_sizes = geometry.size[geometry.clipped]
    assert np.allclose(clipped_sizes, volute.section.minimum_size)
    assert np.all(np.diff(geometry.area) >= 0.0)


def test_constant_area_law_holds_the_passage_uniform():
    target = 8.0e-4
    volute, _ = pump_volute(area_law=ducts.ConstantArea(area=target))
    assert np.allclose(volute.geometry.area, target, rtol=1e-9)


# -- loss model -------------------------------------------------------------


def test_designed_sizing_parameter_is_realised_by_the_analysis():
    """Design Eq. 11-5 and loss Eq. 5-106 are consistent, so SP round-trips."""
    for SP in (0.8, 1.0, 1.5):
        volute, mass_flow = pump_volute(SP=SP)
        assert solve(volute, mass_flow).SP == pytest.approx(SP, rel=1e-9)


def test_tangential_loss_vanishes_at_the_matched_sizing_parameter():
    volute, mass_flow = pump_volute(SP=1.0)
    assert solve(volute, mass_flow).omega_tangential == pytest.approx(0.0, abs=1e-12)


def test_undersizing_is_penalised_far_more_steeply_than_oversizing():
    """Eq. 5-107 diverges as SP falls but saturates as SP grows."""
    losses = {}
    for SP in (0.4, 0.5, 1.0, 2.0, 4.0):
        volute, mass_flow = pump_volute(SP=SP)
        losses[SP] = solve(volute, mass_flow).omega_tangential
    assert losses[0.4] > losses[0.5] > losses[1.0]
    assert losses[1.0] < losses[2.0] < losses[4.0]
    assert losses[0.4] > losses[4.0]


def test_meridional_head_is_lost_in_full():
    """Eq. 5-105 charges the whole through-slot velocity head."""
    volute, mass_flow = pump_volute()
    result = solve(volute, mass_flow)
    expected = (volute.inlet.C_m / volute.inlet.C) ** 2
    assert result.omega_meridional == pytest.approx(expected, rel=1e-12)


def test_exit_cone_loss_vanishes_when_the_pipe_matches_the_collection_plane():
    volute, mass_flow = pump_volute()
    area = volute.geometry.area[-1]
    matched = math.sqrt(4.0 * area / math.pi)
    result = solve(volute, mass_flow, exit_diameter=matched)
    assert result.omega_exit_cone == pytest.approx(0.0, abs=1e-12)
    assert result.d_H == pytest.approx(matched, rel=1e-12)


def test_friction_follows_the_equivalent_diameter_and_mean_path_length():
    """Eq. 5-109, with c_f the Fanning factor and L the mean spiral path."""
    volute, mass_flow = pump_volute()
    result = solve(volute, mass_flow, roughness=1.0e-5)
    expected_length = math.pi * (volute.inlet.radius + result.r_6) / 2.0
    assert result.path_length == pytest.approx(expected_length, rel=1e-12)
    expected = (
        4.0
        * result.friction_factor
        * (result.C_6 / result.C_5) ** 2
        * result.path_length
        / result.d_H
    )
    assert result.omega_friction == pytest.approx(expected, rel=1e-12)


def test_total_loss_is_the_sum_of_its_mechanisms_on_the_inlet_head():
    volute, mass_flow = pump_volute()
    result = solve(volute, mass_flow)
    assert result.omega_total == pytest.approx(sum(result.loss_breakdown().values()))
    expected = result.omega_total * 1000.0 * result.C_5**2 / 2.0
    assert result.pressure_loss == pytest.approx(expected, rel=1e-9)
    assert result.total_pressure_loss == result.pressure_loss
    assert result.states[1].p == pytest.approx(
        result.states[0].p - result.pressure_loss, rel=1e-9
    )


def test_a_collector_loses_half_the_tangential_head_even_when_well_sized():
    """Eq. 5-108: a uniform passage diffuses toward rest where collection starts."""
    shaped, mass_flow = pump_volute(SP=1.0)
    area = shaped.geometry.area[-1]
    collector, _ = pump_volute(area_law=ducts.ConstantArea(area=area))
    result = solve(collector, mass_flow)
    assert result.SP == pytest.approx(1.0, rel=1e-9)
    assert result.omega_tangential == pytest.approx(0.5, rel=1e-9)
    assert result.pressure_loss > solve(shaped, mass_flow).pressure_loss


def test_geometry_is_reusable_across_operating_points():
    volute, mass_flow = pump_volute()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ducts.VoluteConsistencyWarning)
        off_design = solve(volute, 0.7 * mass_flow)
    assert off_design.geometry is volute.geometry
    assert off_design.SP > solve(volute, mass_flow).SP


def test_inlet_velocities_inconsistent_with_mass_flow_are_reported():
    volute, mass_flow = pump_volute()
    with pytest.warns(ducts.VoluteConsistencyWarning, match="inlet slot"):
        solve(volute, 1.5 * mass_flow)


def test_a_swirl_free_inlet_is_rejected():
    """Without incoming swirl there is no angular momentum to conserve."""
    with pytest.raises(ValueError, match="C_u"):
        ducts.VoluteInlet(radius=0.06, C_m=6.0, C_u=0.0)


# -- three-dimensional geometry ---------------------------------------------


def test_sweep_places_an_untilted_volute_in_the_plane_of_its_inlet_circle():
    volute, _ = pump_volute()
    # Densely sampled, so the bore's extreme point is actually resolved.
    rings = volute.sweep(n_points=512)
    assert rings.shape[0] == len(volute.geometry.theta)
    radius = np.hypot(rings[..., 0], rings[..., 1])
    # The slot mouth spans the axis only, so nothing reaches inboard of r5.
    assert radius.min() == pytest.approx(volute.inlet.radius, rel=1e-9)
    # Untilted, the passage's half-width about the inlet plane is its bore.
    assert rings[..., 2].max() == pytest.approx(volute.geometry.size[-1], rel=1e-4)
    assert rings[..., 2].min() == pytest.approx(-volute.geometry.size[-1], rel=1e-4)


def test_tilt_rotates_the_section_plane_without_changing_the_analysis():
    """Tilt is placement: it moves the passage but not a single loss scalar."""
    section = ducts.EllipticalSlottedSection(0.012, 0.003, 0.002)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    flat = ducts.Volute(section, inlet, SHAPED, rotation="ccw")
    tilted = ducts.Volute(section, inlet, SHAPED, rotation="ccw", tilt=math.radians(45.0))
    assert np.allclose(flat.geometry.area, tilted.geometry.area)
    mass_flow = 1000.0 * inlet.C_m * 2.0 * math.pi * inlet.radius * section.inlet_width
    assert solve(tilted, mass_flow).pressure_loss == pytest.approx(
        solve(flat, mass_flow).pressure_loss, rel=1e-12
    )
    # Distance from the mouth centre is preserved; only its direction turns.
    origin = inlet.radius * np.array([1.0, 0.0, 0.0])
    for volute in (flat, tilted):
        offset = volute.sweep(n_points=24)[0] - origin
        assert np.isfinite(offset).all()
    assert not np.allclose(flat.sweep(n_points=24), tilted.sweep(n_points=24))


def test_conserved_angular_momentum_shows_up_as_a_slowing_passage():
    """Velocity times centroid radius is what the sizing law holds constant."""
    volute, mass_flow = pump_volute()
    result = solve(volute, mass_flow)
    geometry = result.geometry
    unclipped = ~geometry.clipped
    velocity = result.collected_velocity()
    momentum = velocity[unclipped] * geometry.centroid_radius[unclipped]
    assert np.allclose(momentum, momentum[0], rtol=1e-6)
    assert velocity[-1] == pytest.approx(result.C_6, rel=1e-9)
    assert velocity[unclipped][0] > velocity[-1]


def test_the_bore_centre_is_a_valid_fan_origin_where_the_centroid_is_not():
    """Capping a swept passage assumes star-shapedness; only one point has it."""
    section = ducts.EllipticalSlottedSection(0.012, 0.004, throat_length=0.040)
    sizes = np.array([0.0065, 0.008, 0.012, 0.02])

    def min_normalised_cross(origin):
        ring = section.polygon(sizes, n_points=400)[:, :-1]
        a = ring - origin[:, None, :]
        b = np.roll(ring, -1, axis=1) - origin[:, None, :]
        cross = a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]
        return (cross / section.area(sizes)[:, None]).min()

    assert min_normalised_cross(section.fan_origin(sizes)) > -1e-12
    # A long throat drags the area centroid inboard of the fillets.
    assert min_normalised_cross(section.metrics(sizes)[1]) < -1e-6


# -- exit cone --------------------------------------------------------------


def test_volute_without_an_exit_cone_requires_an_explicit_exit_diameter():
    volute, mass_flow = pump_volute()
    assert volute.exit_geometry is None
    with pytest.raises(ValueError, match="exit_diameter is required"):
        volute.solve(inlet_state=inlet_state(), mass_flow=mass_flow,
                     transport=ConstantTransport())


def test_an_attached_exit_cone_supplies_the_default_exit_diameter():
    section = ducts.EllipticalSlottedSection(0.012, 0.003, 0.002)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    cone = ducts.ExitCone(final_diameter=0.05, first_angle="radial")
    volute = ducts.Volute(section, inlet, SHAPED, rotation="ccw", exit_cone=cone, n_stations=61)
    mass_flow = 1000.0 * inlet.C_m * 2.0 * math.pi * inlet.radius * section.inlet_width

    defaulted = solve(volute, mass_flow, exit_diameter=None)
    explicit = solve(volute, mass_flow, exit_diameter=cone.final_diameter)
    assert defaulted.pressure_loss == pytest.approx(explicit.pressure_loss)

    # An explicit exit_diameter still overrides the attached cone.
    overridden = solve(volute, mass_flow, exit_diameter=0.08)
    assert overridden.C_7 != pytest.approx(defaulted.C_7)


def test_outlet_center_is_the_bore_centre_not_the_area_centroid():
    """The fan origin (bore centre) and the area centroid are different points.

    A slotted section is not point-symmetric about its bore centre -- the slot
    bulges to one side only -- so the two radii genuinely differ; this is not
    an approximation either way.
    """
    volute, _ = pump_volute()
    centre = volute.outlet_center()
    assert centre.shape == (3,)
    assert centre[2] == pytest.approx(0.0, abs=1e-12)

    size = volute.geometry.size[-1:]
    bore_radius = volute.inlet.radius + volute.section.fan_origin(size)[0, 0]
    assert np.hypot(centre[0], centre[1]) == pytest.approx(bore_radius, rel=1e-9)
    assert bore_radius != pytest.approx(volute.geometry.centroid_radius[-1], rel=1e-3)


def test_sweep_exit_cone_starts_at_the_volute_outlet_and_ends_at_the_flange():
    """The first ring is bit-identical to the volute's own last one.

    ``n_points`` means the same thing here as in :meth:`Volute.sweep` -- points
    per boundary piece, six pieces total -- specifically so the two calls at
    matching ``n_points`` produce rings that stitch into one seamless duct
    rather than merely tracing the same curve at different points along it.
    """
    section = ducts.EllipticalSlottedSection(0.012, 0.003, 0.002)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    cone = ducts.ExitCone(final_diameter=0.05, first_angle="radial",
                          arc_length_parameter=3.0)
    volute = ducts.Volute(section, inlet, SHAPED, rotation="ccw", exit_cone=cone, n_stations=61)

    n_points = 48
    rings = volute.sweep_exit_cone(n_points=n_points, n_stations=13)
    assert rings.shape == (13, 6 * n_points, 3)

    volute_last_ring = volute.sweep(n_points=n_points)[-1, :-1]
    assert rings[0] == pytest.approx(volute_last_ring, abs=1e-12)

    # The last ring is a true circle at the flange diameter.
    last_radii = np.linalg.norm(rings[-1] - volute.exit_geometry.end, axis=-1)
    assert last_radii == pytest.approx(cone.final_diameter/2.0, rel=1e-9)


def test_sweep_exit_cone_requires_an_attached_exit_cone():
    volute, _ = pump_volute()
    with pytest.raises(ValueError, match="no exit cone"):
        volute.sweep_exit_cone()


def test_a_tangential_exit_cone_produces_a_straight_swept_duct():
    """Every ring lies exactly in the plane through its centreline point.

    The frame vectors are always perpendicular to the local travel direction
    (by construction, and trivially so for a straight duct where the frame
    never rotates), so ``(ring - centreline_point) . axis`` is exactly zero at
    every vertex -- unlike the mean of the ring's vertices, which need not sit
    on the axis for a non-point-symmetric cross-section.
    """
    section = ducts.EllipticalSlottedSection(0.012, 0.003, 0.002)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    cone = ducts.ExitCone(final_diameter=0.05, first_angle="tangential")
    volute = ducts.Volute(section, inlet, SHAPED, rotation="ccw", exit_cone=cone, n_stations=61)
    assert volute.exit_geometry.is_straight

    n_stations = 9
    rings = volute.sweep_exit_cone(n_points=32, n_stations=n_stations)
    points, _, _ = volute.exit_geometry.frames(n_stations)
    axis = volute.exit_geometry.start_direction
    offsets = rings - points[:, None, :]
    assert offsets @ axis == pytest.approx(0.0, abs=1e-9)


def test_sweep_exit_cone_matches_a_tilted_volutes_own_section_frame():
    """A tilted volute's cutwater must continue smoothly, not snap to meet it.

    The exit cone's own frame math is tilt-agnostic (correctly, since
    tangential/radial/axial are machine-frame directions), but the *rendered*
    starting frame must be the volute's own tilted section basis, not the
    untilted default, or the first ring would not match the volute's own last
    ring for any nonzero tilt.
    """
    section = ducts.EllipticalSlottedSection(0.012, 0.003, 0.002)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    cone = ducts.ExitCone(final_diameter=0.05, first_angle="radial",
                          second_angle=math.radians(15), arc_length_parameter=2.5)
    volute = ducts.Volute(section, inlet, SHAPED, rotation="ccw", exit_cone=cone,
                          tilt=math.radians(45.0), n_stations=61)

    n_points = 24
    volute_last_ring = volute.sweep(n_points=n_points)[-1, :-1]
    exit_first_ring = volute.sweep_exit_cone(n_points=n_points, n_stations=9)[0]
    assert exit_first_ring == pytest.approx(volute_last_ring, abs=1e-12)


# -- elliptical sections and tongue flattening ------------------------------


def test_unit_aspect_reproduces_the_circular_closed_form():
    """The ellipse is solved numerically; at aspect one it must match algebra.

    The bore centre of a circular section has a closed form, so it pins the
    bisected tangency solve independently of the rest of the construction.
    """
    width, fillet, throat = 0.012, 0.003, 0.002
    section = ducts.EllipticalSlottedSection(width, fillet, throat)
    sizes = np.linspace(width / 2.0, 0.15, 60)
    h = width / 2.0 + fillet
    expected = throat + np.sqrt(np.maximum((sizes + fillet) ** 2 - h**2, 0.0))
    assert section.fan_origin(sizes)[:, 0] == pytest.approx(expected, abs=1e-9)


def test_aspect_stretches_the_bore_across_the_slot():
    section = ducts.EllipticalSlottedSection(0.012, 0.003, max_aspect_ratio=4.0)
    size = 0.03
    # Densely sampled, so the bore's extreme point is actually resolved.
    for aspect in (0.5, 1.0, 2.0, 4.0):
        ring = section.polygon(size, n_points=2000, aspect=aspect)[0]
        across = ring[:, 1].max() - ring[:, 1].min()
        # The semi-axis across the slot is the size itself, whatever the aspect.
        assert across == pytest.approx(2.0 * size, rel=1e-6)
    # Reaching outward is what the aspect ratio trades away.
    outward = [section.polygon(size, n_points=2000, aspect=a)[0][:, 0].max()
               for a in (0.5, 1.0, 2.0, 4.0)]
    assert outward[0] > outward[1] > outward[2] > outward[3]


def test_flattening_shrinks_the_minimum_area_inversely_with_aspect():
    """With no throat rectangle to floor it, area scales as 1/aspect."""
    section = ducts.EllipticalSlottedSection(0.012, 0.003, throat_length=0.0,
                                             max_aspect_ratio=8.0)
    minimum = section.minimum_size
    circle = section.area(minimum, aspect=1.0)[0]
    for aspect in (2.0, 4.0, 8.0):
        assert section.area(minimum, aspect=aspect)[0] == pytest.approx(
            circle / aspect, rel=1e-6
        )


def test_the_volute_flattens_toward_the_tongue_instead_of_clipping():
    """Aspect is the second lever once size has bottomed out on the slot."""
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    circle = ducts.Volute(
        ducts.EllipticalSlottedSection(0.012, 0.003, 0.002, max_aspect_ratio=1.0),
        inlet, SHAPED, rotation="ccw", n_stations=181)
    ellipse = ducts.Volute(
        ducts.EllipticalSlottedSection(0.012, 0.003, 0.002, max_aspect_ratio=4.0),
        inlet, SHAPED, rotation="ccw", n_stations=181)

    # Flattening buys a materially shorter clipped run at the tongue.
    assert ellipse.geometry.tongue_clip_angle < circle.geometry.tongue_clip_angle
    assert ellipse.geometry.flattened.any()
    assert not circle.geometry.flattened.any()

    # Where it flattens, size sits on its floor and only aspect moves.
    flattened = ellipse.geometry.flattened
    assert np.allclose(ellipse.geometry.size[flattened],
                       ellipse.section.minimum_size)
    assert np.all(ellipse.geometry.aspect[flattened] > 1.0)
    assert np.all(ellipse.geometry.aspect[flattened] <= 4.0 + 1e-9)

    # Beyond the flattened run it is back to the design shape.
    assert ellipse.geometry.aspect[-1] == pytest.approx(1.0)
    # And the passage tracks the area law more closely the whole way in.
    assert np.all(ellipse.geometry.area <= circle.geometry.area + 1e-15)


def test_flattening_is_monotone_and_continuous_into_the_design_aspect():
    """Aspect falls smoothly to the design value as the passage opens out.

    The flattened aspect goes roughly as ``1/theta``, so it is genuinely steep
    near the tongue and the station-to-station step there is large at any fixed
    resolution. Continuity is the statement that the step shrinks with station
    spacing, not that it is small at one particular spacing.
    """
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)

    def aspect_at(n_stations):
        section = ducts.EllipticalSlottedSection(0.012, 0.003, 0.002,
                                                 max_aspect_ratio=4.0)
        return ducts.Volute(section, inlet, SHAPED, rotation="ccw", n_stations=n_stations).geometry.aspect

    coarse, fine = aspect_at(181), aspect_at(721)
    # Never flattens further as the passage opens out, and lands on design.
    assert np.all(np.diff(coarse) <= 1e-9)
    assert coarse.min() == pytest.approx(1.0)
    # Refining the sweep shrinks the largest step, so there is no real jump.
    assert np.max(np.abs(np.diff(fine))) < 0.5 * np.max(np.abs(np.diff(coarse)))


def test_an_oblong_design_aspect_still_flattens_from_its_own_shape():
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    volute = ducts.Volute(
        ducts.EllipticalSlottedSection(0.012, 0.003, 0.002,
                                       aspect_ratio=2.0, max_aspect_ratio=6.0),
        inlet, SHAPED, rotation="ccw", n_stations=121)
    assert volute.geometry.aspect.min() == pytest.approx(2.0)
    assert volute.geometry.aspect.max() == pytest.approx(6.0)
    assert volute.geometry.flattened.any()


def test_max_aspect_ratio_below_the_design_aspect_is_rejected():
    with pytest.raises(ValueError, match="at least aspect_ratio"):
        ducts.EllipticalSlottedSection(0.012, 0.003, aspect_ratio=3.0,
                                       max_aspect_ratio=2.0)


def test_a_family_that_cannot_flatten_clips_exactly_as_before():
    """max_aspect_ratio == aspect_ratio restores the fixed-shape behaviour."""
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    volute = ducts.Volute(
        ducts.EllipticalSlottedSection(0.012, 0.003, 0.002, max_aspect_ratio=1.0),
        inlet, SHAPED, rotation="ccw", n_stations=121)
    assert not volute.geometry.flattened.any()
    assert np.allclose(volute.geometry.aspect, 1.0)
    clipped = volute.geometry.clipped
    assert clipped.any()
    assert np.allclose(volute.geometry.size[clipped], volute.section.minimum_size)


def test_the_sizing_law_still_holds_wherever_it_was_not_clipped():
    """Flattened stations satisfy the area law as exactly as sized ones do."""
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    section = ducts.EllipticalSlottedSection(0.012, 0.003, 0.002,
                                             max_aspect_ratio=4.0)
    geometry = ducts.Volute(section, inlet, SHAPED, rotation="ccw", n_stations=181).geometry
    live = ~geometry.clipped
    target = (geometry.theta[live] * geometry.centroid_radius[live]
              * section.inlet_width * inlet.tan_alpha)
    assert geometry.area[live] == pytest.approx(target, rel=1e-6)
    # Which includes the flattened run, not just the sized one.
    assert (geometry.flattened & live).any()


# -- exit-cone duct losses and pressure --------------------------------------


def with_cone(cone, **kwargs):
    section = ducts.EllipticalSlottedSection(0.012, 0.003, 0.002)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    volute = ducts.Volute(section, inlet, SHAPED, rotation="ccw", exit_cone=cone, n_stations=61, **kwargs)
    mass_flow = 1000.0 * inlet.C_m * 2.0 * math.pi * inlet.radius * section.inlet_width
    return volute, solve(volute, mass_flow, exit_diameter=None)


def test_a_volute_without_an_exit_cone_is_charged_no_duct_loss():
    volute, mass_flow = pump_volute()
    result = solve(volute, mass_flow)
    assert result.omega_exit_friction == 0.0
    assert result.omega_exit_bend == 0.0


def test_a_straight_duct_is_charged_friction_but_no_bend():
    _, result = with_cone(ducts.ExitCone(0.05, "tangential", arc_length_parameter=2.0))
    assert result.omega_exit_friction > 0.0
    assert result.omega_exit_bend == 0.0


def test_duct_friction_grows_with_its_length():
    _, short = with_cone(ducts.ExitCone(0.05, "tangential", arc_length_parameter=2.0))
    _, long = with_cone(ducts.ExitCone(0.05, "tangential", arc_length_parameter=6.0))
    # Three times the duct, near enough three times the friction. Only near
    # enough: lengthening a cone also shallows it, and the correlation charges
    # friction over the sloped wall rather than a cylinder of the mean bore.
    assert long.omega_exit_friction == pytest.approx(
        3.0 * short.omega_exit_friction, rel=0.05
    )


def test_a_turning_duct_is_charged_a_bend_loss():
    _, straight = with_cone(ducts.ExitCone(0.05, "tangential", arc_length_parameter=2.5))
    _, turning = with_cone(ducts.ExitCone(0.05, "radial", arc_length_parameter=2.5))
    assert turning.omega_exit_bend > 0.0
    assert turning.pressure_loss > straight.pressure_loss


def test_a_gentler_bend_costs_less_to_turn_and_more_to_rub():
    """Spreading the same turn over more duct trades bend loss for friction."""
    _, tight = with_cone(ducts.ExitCone(0.05, "radial", arc_length_parameter=2.5))
    _, gentle = with_cone(ducts.ExitCone(0.05, "radial", arc_length_parameter=5.0))
    assert gentle.omega_exit_bend < tight.omega_exit_bend
    assert gentle.omega_exit_friction > tight.omega_exit_friction


def test_a_gentler_cone_diffuses_at_lower_loss():
    """The angle the cone opens out at is most of what decides its loss.

    The whole point of spending duct length on a diffuser, and the reason the
    diffusion term cannot be a function of the area ratio alone.
    """
    _, abrupt = with_cone(ducts.ExitCone(0.05, "tangential", arc_length_parameter=1.0))
    _, gentle = with_cone(ducts.ExitCone(0.05, "tangential", arc_length_parameter=6.0))
    assert gentle.omega_exit_cone < abrupt.omega_exit_cone
    # Even an abrupt cone beats dumping the passage straight into the flange.
    assert abrupt.omega_exit_cone < ((abrupt.C_6 - abrupt.C_7) / abrupt.C_5) ** 2


def test_a_narrowing_cone_is_charged_as_a_contraction():
    """A flange narrower than the passage speeds the flow up, which is cheap.

    Borda-Carnot is not merely angle-blind here, it is the wrong shape of
    formula: it is derived for a flow that separates on expanding, and nothing
    separates on the way into a contraction.
    """
    _, result = with_cone(ducts.ExitCone(0.025, "tangential", arc_length_parameter=2.0))
    assert result.C_7 > result.C_6
    borda = ((result.C_6 - result.C_7) / result.C_5) ** 2
    assert result.omega_exit_cone < 0.1 * borda


def test_a_gentler_contraction_costs_less_to_squeeze_through():
    _, abrupt = with_cone(ducts.ExitCone(0.025, "tangential", arc_length_parameter=2.0))
    _, gentle = with_cone(ducts.ExitCone(0.025, "tangential", arc_length_parameter=6.0))
    assert gentle.omega_exit_cone < abrupt.omega_exit_cone


def test_without_a_cone_the_flange_is_charged_aungiers_dump_loss():
    """No cone to resolve, so Eq. 5-110's sudden expansion is the right model."""
    volute, mass_flow = pump_volute()
    result = solve(volute, mass_flow)
    assert result.omega_exit_cone == pytest.approx(
        ((result.C_6 - result.C_7) / result.C_5) ** 2
    )


def test_the_loss_breakdown_sums_to_the_total():
    _, result = with_cone(ducts.ExitCone(0.05, "radial", arc_length_parameter=2.5))
    assert result.omega_total == pytest.approx(sum(result.loss_breakdown().values()))


def test_stagnation_pressure_spans_exactly_the_collection_loss():
    volute, mass_flow = pump_volute()
    result = solve(volute, mass_flow)
    stagnation = result.stagnation_pressure()
    assert stagnation[0] == pytest.approx(result.states[0].p)
    assert np.all(np.diff(stagnation) <= 1e-9)
    # The passage carries every loss except the exit cone's.
    assert result.states[0].p - stagnation[-1] == pytest.approx(result.collection_loss())


def test_static_pressure_is_highest_at_the_still_tongue():
    volute, mass_flow = pump_volute()
    result = solve(volute, mass_flow)
    static = result.static_pressure()
    velocity = result.collected_velocity()
    # Nothing collected yet, so nothing is moving and none of the head is lost.
    assert static[0] == pytest.approx(result.states[0].p)
    assert np.argmax(static) == 0
    # Lowest where the passage runs fastest, which is just past the tongue.
    assert np.argmin(static) == np.argmax(velocity)


def test_exit_cone_pressure_lands_on_the_solved_outlet():
    _, result = with_cone(ducts.ExitCone(0.05, "radial", arc_length_parameter=2.5))
    stagnation, static = result.exit_cone_pressures(7)
    assert stagnation[0] == pytest.approx(result.stagnation_pressure()[-1])
    assert stagnation[-1] == pytest.approx(result.states[1].p)
    assert np.all(np.diff(stagnation) <= 1e-9)
    # A diffusing duct trades stagnation pressure for static pressure.
    assert result.C_7 < result.C_6
    assert static[-1] > static[0]


# -- rotation sense ------------------------------------------------------


def rotation_pair(**kwargs):
    """The same volute wound each way."""
    section = ducts.EllipticalSlottedSection(
        inlet_width=0.012, fillet_radius=0.003, throat_length=0.002
    )
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    return tuple(
        ducts.Volute(section, inlet, SHAPED, rotation=sense, n_stations=41, **kwargs)
        for sense in ("ccw", "cw")
    )


def test_rotation_mirrors_the_passage_and_leaves_the_sizing_untouched():
    ccw, cw = rotation_pair()
    # The area law is solved in the azimuth alone, so the passage it sizes is
    # the same hardware either way round.
    for field in ("theta", "size", "aspect", "area", "centroid_radius"):
        assert np.allclose(getattr(ccw.geometry, field),
                           getattr(cw.geometry, field))
    mirrored = cw.sweep(n_points=24).copy()
    mirrored[..., 1] *= -1.0
    assert ccw.sweep(n_points=24) == pytest.approx(mirrored, abs=1e-15)


def test_rotation_costs_no_pressure_the_way_tilt_does_not():
    ccw, cw = rotation_pair()
    mass_flow = 1000.0 * 6.0 * 2.0 * math.pi * 0.060 * 0.012
    assert solve(cw, mass_flow).pressure_loss == pytest.approx(
        solve(ccw, mass_flow).pressure_loss, rel=1e-12
    )


def test_a_clockwise_exit_cone_turns_the_same_amount_the_other_way():
    cone = ducts.ExitCone(final_diameter=0.05, first_angle="radial",
                          arc_length_parameter=3.0)
    ccw, cw = rotation_pair(exit_cone=cone)
    # Every scalar the loss model reads is a turn or a length, and mirroring
    # changes neither.
    for field in ("first_angle", "radial_angle", "turn_angle", "bend_radius",
                  "arc_length", "end_azimuth"):
        assert getattr(cw.exit_geometry, field) == pytest.approx(
            getattr(ccw.exit_geometry, field), rel=1e-12
        )
    for field in ("end", "end_direction", "start_direction"):
        expected = np.array(getattr(ccw.exit_geometry, field)) * [1.0, -1.0, 1.0]
        assert getattr(cw.exit_geometry, field) == pytest.approx(expected, abs=1e-15)


def test_a_clockwise_duct_still_stitches_onto_its_own_volute():
    cone = ducts.ExitCone(final_diameter=0.05, first_angle="radial",
                          arc_length_parameter=3.0)
    for volute in rotation_pair(exit_cone=cone, tilt=math.radians(30.0)):
        last = volute.sweep(n_points=24)[-1, :-1]
        first = volute.sweep_exit_cone(n_points=24, n_stations=9)[0]
        assert first == pytest.approx(last, abs=1e-12)


@pytest.mark.parametrize("sense", ["cw", "CW", "Clockwise", " clockwise "])
def test_clockwise_is_spelled_however_the_designer_spells_it(sense):
    section = ducts.EllipticalSlottedSection(inlet_width=0.012, fillet_radius=0.003)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    volute = ducts.Volute(section, inlet, SHAPED, rotation=sense, n_stations=9)
    assert volute.rotation == "clockwise"


def test_an_unrecognised_rotation_is_refused():
    section = ducts.EllipticalSlottedSection(inlet_width=0.012, fillet_radius=0.003)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    with pytest.raises(ValueError, match="rotation must be one of"):
        ducts.Volute(section, inlet, SHAPED, rotation="widdershins", n_stations=9)


def test_the_area_law_and_the_rotation_are_both_stated_not_defaulted():
    section = ducts.EllipticalSlottedSection(inlet_width=0.012, fillet_radius=0.003)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    # A shaped volute and a uniform collector are different hardware, and which
    # way a volute winds has to match the machine feeding it, so neither is a
    # sensible default to fall back on.
    with pytest.raises(TypeError):
        ducts.Volute(section, inlet, rotation="ccw")
    with pytest.raises(TypeError):
        ducts.Volute(section, inlet, SHAPED)
