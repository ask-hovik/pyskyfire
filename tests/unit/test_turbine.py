"""Axial turbine meanline, hardware definition and built geometry.

The meanline is checked against what it has to satisfy whatever the inputs:
Euler's work equals the efficiency times the isentropic head, continuity
holds at every station on the drawn annulus, the last station lands on the
duty's outlet pressure, and a 50 % reaction repeating stage has the symmetric
triangles every textbook draws. The built geometry is checked against a thin
cascade's throat, ``pitch*cos(exit angle)``.
"""

import math
import pickle

import numpy as np
import pytest

from pyskyfire.common import MeridionalPassage
from pyskyfire.thermo import CoolPropEOS, Incompressible
from pyskyfire.turbine import (
    AxialTurbine,
    BladeRow,
    FreeVortexTwist,
    RowSettings,
    SectionProfile,
    Stage,
    ThicknessDistribution,
    TurbineDefinition,
    TurbineDuty,
    metal_angles,
)


@pytest.fixture(scope="module")
def duty():
    return TurbineDuty.from_states(fluid="Nitrogen", mdot=20., p0_in=40e5, T0_in=700.,
                                   p_out=28e5, n=30000.)


def design(duty, **kwargs):
    kwargs.setdefault("stages", [Stage(0.6, reaction=0.5)])
    kwargs.setdefault("mean_diameter", 0.2)
    return AxialTurbine(duty, kwargs.pop("efficiency", 0.85), kwargs.pop("stages"), **kwargs)


# -- equations of state used by the meanline --------------------------------------

def test_state_from_enthalpy_and_entropy_round_trips_even_past_the_eos_range():
    eos = CoolPropEOS("Methane")
    # 633 K is past the methane equation's 625 K, where CoolProp's own (h, s)
    # flash refuses and the isenthalp search takes over.
    for p, T in ((50e5, 300.), (160e5, 633.)):
        state = eos.state_pT(p, T)
        back = eos.state_hs(state.h, eos.entropy(state))
        assert back.p == pytest.approx(p, rel=1e-8) and back.T == pytest.approx(T, rel=1e-8)
    assert eos.speed_of_sound(eos.state_pT(1e5, 300.)) == pytest.approx(450., rel=0.02)


def test_incompressible_entropy_needs_cp():
    with pytest.raises(ValueError, match="cp"):
        Incompressible(1000.).entropy(Incompressible(1000.).state_pT(1e5, 300.))
    eos = Incompressible(1000., cp=4000.)
    state = eos.state_pT(3e5, 290.)
    back = eos.state_hs(state.h, eos.entropy(state))
    assert back.p == pytest.approx(3e5) and back.T == pytest.approx(290.)
    assert eos.speed_of_sound(state) == math.inf


# -- meanline -------------------------------------------------------------------------

def test_euler_work_is_efficiency_times_isentropic_head(duty):
    t = design(duty, efficiency=0.8)
    assert t.stages[0].work == pytest.approx(0.8*duty.isentropic_head, rel=1e-12)
    assert t.stages[0].reaction == pytest.approx(0.5, abs=1e-9)
    assert t.stages[0].flow_coefficient == pytest.approx(0.6)


def test_last_station_lands_on_the_outlet_pressure(duty):
    assert design(duty).stations[-1].p0 == pytest.approx(duty.outlet.p, rel=1e-7)
    static = TurbineDuty.from_states(fluid="Nitrogen", mdot=20., p0_in=40e5, T0_in=700.,
                                     p_out=28e5, n=30000., outlet_basis="static")
    assert design(static, efficiency=0.7).stations[-1].static.p == pytest.approx(28e5, rel=1e-7)
    # Total-to-static, the loss has to cover the leaving kinetic energy at least.
    with pytest.raises(ValueError, match="exit velocity"):
        design(static, efficiency=0.95)


def test_continuity_holds_on_the_drawn_annulus(duty):
    t = design(duty, stages=[Stage(0.5, reaction=0.3), Stage(0.7, reaction=0.5)])
    rows = t.definition.rows
    edges = [(rows[0].r_hub_in, rows[0].r_tip_in)] + [(r.r_hub_out, r.r_tip_out) for r in rows]
    for station, (hub, tip) in zip(t.stations, edges):
        mdot = station.static.rho*station.c_x*np.pi*(tip**2 - hub**2)
        assert mdot == pytest.approx(duty.mdot, rel=1e-12)


def test_fifty_percent_repeating_stage_has_symmetric_triangles(duty):
    again = design(duty)
    for _ in range(60):
        again = design(duty, inlet_swirl=again.exit_swirl)
    s = again.stages[0]
    assert s.inlet.c_u == pytest.approx(s.exit.c_u, rel=1e-6)
    assert s.middle.alpha == pytest.approx(-s.exit.beta, rel=1e-6)
    assert s.inlet.alpha == pytest.approx(-s.middle.beta, rel=1e-6)
    # The textbook swirls for R = 1/2 with c1 = c3.
    u, psi = s.middle.u, s.loading
    assert s.middle.c_u == pytest.approx(u*(0.5 + psi/2), rel=1e-6)


def test_reaction_and_exit_swirl_routes_agree(duty):
    axial = design(duty, stages=[Stage(0.6, exit_swirl=0.0)])
    assert axial.exit_swirl == 0.0 and axial.exit_flow_angle == 0.0
    same = design(duty, stages=[Stage(0.6, reaction=axial.stages[0].reaction)])
    assert same.exit_swirl == pytest.approx(0.0, abs=1e-6)
    assert same.stages[0].middle.c_u == pytest.approx(axial.stages[0].middle.c_u)


def test_blade_jet_ratio_sets_the_mean_diameter(duty):
    t = design(duty, mean_diameter=None, blade_jet_ratio=0.5)
    assert t.blade_speed/duty.spouting_velocity == pytest.approx(0.5)


def test_work_and_loss_splits(duty):
    t = design(duty, stages=[Stage(0.6, reaction=0.5)]*2, work_split=[2, 1], loss_split=[1, 1, 1, 3])
    assert t.stages[0].work == pytest.approx(2*t.stages[1].work)
    s = [t.duty.eos.entropy(st.static) for st in t.stations]
    rises = np.diff(s)
    np.testing.assert_allclose(rises/rises.sum(), [1/6, 1/6, 1/6, 1/2], rtol=1e-6)


@pytest.mark.parametrize("kwargs, message", [
    (dict(stages=[]), "non-empty"),
    (dict(efficiency=1.2), "efficiency"),
    (dict(mean_diameter=None), "one of the two"),
    (dict(work_split=[1, 2]), "work_split"),
])
def test_invalid_meanline_inputs(duty, kwargs, message):
    with pytest.raises(ValueError, match=message):
        design(duty, **kwargs)


def test_stage_and_row_settings_take_one_of_each_pair():
    with pytest.raises(ValueError, match="one of the two"):
        Stage(0.6)
    with pytest.raises(ValueError, match="one of the two"):
        Stage(0.6, reaction=0.5, exit_swirl=0.)
    with pytest.raises(ValueError, match="not both"):
        RowSettings(pitch_chord=0.8, zweifel=0.8)
    assert RowSettings().zweifel == 0.8 and RowSettings(pitch_chord=.8).zweifel is None


# -- hardware ---------------------------------------------------------------------------

def test_rows_follow_the_settings(duty):
    t = design(duty, rows=RowSettings(aspect_ratio=2., pitch_chord=0.75, axial_gap=0.5))
    stator, rotor = t.definition.rows
    assert stator.aspect_ratio == pytest.approx(2.)
    assert rotor.x_le - stator.x_te == pytest.approx(0.5*stator.axial_chord)
    for row in (stator, rotor):
        assert row.pitch/row.chord == pytest.approx(0.75, rel=0.05)
    s = t.stages[0]
    assert stator.outlet_metal_angle == pytest.approx(s.middle.alpha)
    assert rotor.inlet_metal_angle == pytest.approx(s.middle.beta)
    assert rotor.outlet_metal_angle == pytest.approx(s.exit.beta)


def test_stagger_fraction_places_the_chord_between_the_metal_angles(duty):
    for fraction in (0.0, 0.7, 1.0):
        profile = SectionProfile(stagger_fraction=fraction)
        for row in design(duty, rows=RowSettings(stator=profile, rotor=profile)).definition.rows:
            expected = row.inlet_metal_angle + fraction*(row.outlet_metal_angle - row.inlet_metal_angle)
            assert row.stagger == pytest.approx(expected)
    with pytest.raises(ValueError, match="stagger_fraction"):
        SectionProfile(stagger_fraction=1.2)
    fixed = design(duty, rows=RowSettings(rotor=SectionProfile(stagger_deg=-30.))).definition.rows[1]
    assert fixed.stagger == pytest.approx(math.radians(-30.))


def test_stator_and_rotor_take_their_own_section_profiles(duty):
    rows = RowSettings(stator=SectionProfile(max_thickness=0.25),
                       rotor=SectionProfile.impulse(max_thickness=0.12))
    stator, rotor = design(duty, rows=rows).definition.rows
    assert stator.max_thickness/stator.chord == pytest.approx(0.25)
    assert rotor.max_thickness/rotor.chord == pytest.approx(0.12)
    assert rotor.leading_edge_diameter/rotor.chord == pytest.approx(0.03)
    assert stator.thickness_distribution == ThicknessDistribution()
    assert rotor.thickness_distribution == ThicknessDistribution.symmetric()


@pytest.mark.parametrize("kwargs, message", [
    (dict(suction=()), "suction"),
    (dict(pressure=(0.1, float("nan"))), "pressure"),
    (dict(suction_extent=0.), "suction_extent"),
    (dict(spacing=0.), "spacing"),
])
def test_invalid_thickness_distributions(kwargs, message):
    with pytest.raises(ValueError, match=message):
        ThicknessDistribution(**kwargs)


def test_zweifel_sets_the_pitch(duty):
    t = design(duty, rows=RowSettings(zweifel=0.9))
    s = t.stages[0]
    stator = t.definition.rows[0]
    a1, a2 = s.inlet.alpha, s.middle.alpha
    zweifel = 2*(stator.pitch/stator.axial_chord)*math.cos(a2)**2*abs(math.tan(a2) - math.tan(a1))
    assert zweifel == pytest.approx(0.9, rel=0.05)


def test_incidence_and_deviation_move_the_metal():
    assert metal_angles(0.1, 0.9, 0.02, 0.05) == pytest.approx((0.12, 0.95))
    assert metal_angles(0.1, -0.9, 0.0, 0.05) == pytest.approx((0.1, -0.95))


def test_free_vortex_twist():
    twist = FreeVortexTwist(r_mean=.1, c_x_in=100., c_x_out=100., c_u_in=0., c_u_out=150., u_mean=0.)
    assert twist(0.5, .1) == pytest.approx((0., math.atan(1.5)))
    assert math.tan(twist(1., .12)[1]) == pytest.approx(1.5*.1/.12)
    rotor = FreeVortexTwist(.1, 100., 100., 150., 0., u_mean=200.)
    assert math.tan(rotor(0., .08)[0]) == pytest.approx((150*.1/.08 - 200*.8)/100)


def test_passage_records_what_turns_and_every_gap(duty):
    open_rotor = design(duty, rows=RowSettings(tip_clearance=2e-4, hub_clearance=1e-4)).definition
    p = open_rotor.passage
    stator, rotor = open_rotor.rows
    f = lambda x: float(np.interp(x, p.hub[:, 0], p.stations))
    assert p.rotating("hub", f(0.5*(rotor.x_le + rotor.x_te)))
    assert p.rotating("hub", f(0.5*(stator.x_le + stator.x_te)))  # the drum under a cantilevered vane
    assert not p.rotating("shroud").any()
    assert {(c.blade, c.wall, c.gap) for c in p.clearances} == {("R1", "shroud", 2e-4), ("S1", "hub", 1e-4)}
    shrouded = design(duty).definition.passage
    assert shrouded.clearances == ()
    assert shrouded.rotating("shroud", f(0.5*(rotor.x_le + rotor.x_te)))


def test_exit_annulus_matches_the_last_row(duty):
    d = design(duty).definition
    end = d.exit_annulus
    last = d.rows[-1]
    assert (end.r_hub, end.r_shroud) == pytest.approx((last.r_hub_out, last.r_tip_out))
    assert end.area == pytest.approx(np.pi*(last.r_tip_out**2 - last.r_hub_out**2))
    assert end.x == pytest.approx(last.x_te + d.outlet_length)
    np.testing.assert_allclose(end.hub_tangent, [1, 0], atol=1e-9)


def row(**changes):
    base = dict(name="S1", kind="stator", blade_count=40, x_le=0., axial_chord=.01,
                r_hub_in=.09, r_tip_in=.11, r_hub_out=.09, r_tip_out=.11,
                inlet_metal_angle=0., outlet_metal_angle=1.2, stagger=.6,
                leading_edge_diameter=.001, trailing_edge_thickness=.0003, max_thickness=.002)
    base.update(changes)
    return BladeRow(**base)


@pytest.mark.parametrize("changes, message", [
    (dict(kind="nozzle"), "kind"),
    (dict(blade_count=1), "blade_count"),
    (dict(r_tip_out=.08), "above hub"),
    (dict(outlet_metal_angle=1.6), "85"),
    (dict(tip_clearance=1e-4), "casing"),
    (dict(kind="rotor", hub_clearance=1e-4), "hub"),
    (dict(hub_clearance=.03), "smaller than the blade height"),
    (dict(admission=0.), "admission"),
])
def test_invalid_rows(changes, message):
    with pytest.raises(ValueError, match=message):
        row(**changes)


def test_definition_rejects_overlapping_or_repeated_rows():
    with pytest.raises(ValueError, match="overlap"):
        TurbineDefinition((row(), row(name="R1", kind="rotor", x_le=.005)))
    with pytest.raises(ValueError, match="unique"):
        TurbineDefinition((row(), row(x_le=.02)))
    d = TurbineDefinition((row(), row(name="R1", kind="rotor", x_le=.015,
                                      inlet_metal_angle=.2, outlet_metal_angle=-1.)))
    assert isinstance(d.passage, MeridionalPassage)
    assert d.axial_gaps == pytest.approx((.005,))


def test_definition_pickles(duty):
    d = design(duty, rows=RowSettings(tip_clearance=2e-4)).definition
    back = pickle.loads(pickle.dumps(d))
    assert back.rows == d.rows and back.passage.clearances == d.passage.clearances


# -- built geometry --------------------------------------------------------------------

pyturbo = pytest.importorskip("pyturbo")


@pytest.fixture(scope="module")
def built(duty):
    t = design(duty, rows=RowSettings(tip_clearance=2e-4))
    return t, t.definition.build_geometry(span_points=5, chord_points=65)


def test_geometry_has_one_closed_blade_per_row_inside_the_passage(built):
    t, g = built
    assert [r.row.name for r in g.rows] == ["S1", "R1"]
    for rg in g.rows:
        assert rg.blade.suction.shape == (5, 65, 3)
        x = rg.blade.suction[..., 0]
        # The trailing-edge circle stands proud of the camber line's end.
        reach = rg.row.trailing_edge_thickness
        assert x.min() >= rg.row.x_le - 1e-9 and x.max() <= rg.row.x_te + reach
    rotor = g.row("R1")
    tip = rotor.blade.suction[-1, :, 1]
    casing = np.interp(rotor.blade.suction[-1, :, 0], g.passage.shroud[:, 0], g.passage.shroud[:, 1])
    np.testing.assert_allclose(casing - tip, 2e-4, atol=1e-9)


def test_throat_is_measured_and_near_a_thin_cascades(built):
    t, g = built
    for rg in g.rows:
        assert np.all(rg.throat < rg.pitch)
        assert np.all(rg.narrowest <= rg.throat + 1e-12)
        design_exit = rg.at_mid_span(rg.outlet_metal_angle)
        # Thick trailing edges and a curved suction side narrow the throat a
        # little below pitch*cos(exit angle); a few degrees of gauging angle.
        assert rg.at_mid_span(rg.gauging_angle) == pytest.approx(design_exit, abs=math.radians(5))
        np.testing.assert_allclose(rg.max_thickness/rg.row.chord, 0.2, rtol=0.1)
    # A reaction row's passage narrows to its exit, so its narrowest gap is
    # the throat, at the trailing edge.
    assert np.all(g.row("S1").narrowest_at > 0.7)


def test_an_impulse_rotor_drawn_symmetric_does_not_pinch():
    duty = TurbineDuty.from_states(fluid="Oxygen", mdot=26., p0_in=220e5, T0_in=700.,
                                   p_out=140e5, n=30000.)
    rows = RowSettings(aspect_ratio=1.0, stator=SectionProfile(stagger_fraction=0.7),
                       rotor=SectionProfile.impulse())
    t = AxialTurbine(duty, 0.8, [Stage(0.4, reaction=0.0)], blade_jet_ratio=0.47, rows=rows)
    rotor = t.definition.build_geometry(span_points=3, chord_points=65).row("R1")
    np.testing.assert_allclose(rotor.inlet_width, rotor.pitch*np.cos(rotor.inlet_metal_angle))
    # Its steep inlet is as narrow as its exit; drawn with thin edges and its
    # thickness at mid-chord, nothing between them is narrower than either.
    assert np.all(rotor.pinch > 0.95)
    assert np.all((0 <= rotor.narrowest_at) & (rotor.narrowest_at <= 1))
    low, high = sorted((rotor.row.height_in, rotor.row.height_out))
    assert np.all((low <= rotor.narrowest_height) & (rotor.narrowest_height <= high))


def test_free_vortex_twists_the_rotor_toward_the_tip(built):
    _, g = built
    rotor = g.row("R1")
    assert abs(rotor.outlet_metal_angle[-1]) > abs(rotor.outlet_metal_angle[0])


def test_clockwise_turbine_is_the_mirror_image(duty):
    ccw = design(duty).definition.build_geometry(span_points=3, chord_points=33)
    cw = design(duty, rotation="clockwise").definition.build_geometry(span_points=3, chord_points=33)
    for a, b in zip(ccw.rows, cw.rows):
        np.testing.assert_allclose(a.blade.suction[..., 2], -b.blade.suction[..., 2], atol=1e-12)
        np.testing.assert_allclose(a.throat, b.throat)


def test_turbine_scene_repeats_every_blade(built):
    viz = pytest.importorskip("pyskyfire.viz")
    _, g = built
    node = viz.TurbineComponent(g, n_theta=16).build()
    assert len(node.children) == 2 + sum(r.row.blade_count for r in g.rows)
    assert set(node.anchors) == {"inlet", "outlet", "shaft"}
