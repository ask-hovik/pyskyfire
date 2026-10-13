"""Geometry, hydraulic and diagnostic contracts for injectors."""

import math

import numpy as np
import pytest

from pyskyfire.ducts import Orifice
from pyskyfire.injector import (ElementSet, FeedCondition, ImpingingElement, Injector, Jet,
                                JetFlow)
from pyskyfire.thermo import FluidState, Incompressible, PHASE_UNKNOWN

WATER = Incompressible(1000.0)
OIL = Incompressible(800.0)
THIRTY = math.radians(30)


def state(p, T=300.0):
    return FluidState(T=T, p=p, h=math.nan, quality=math.nan, rho=math.nan,
                      phase=PHASE_UNKNOWN, T_sat=math.nan)


def doublet(fuel="fuel", ox="ox", d_fu=1e-3, d_ox=1.2e-3, inclination=THIRTY):
    return ImpingingElement([
        Jet(fuel, Orifice.circular(d_fu, cd=0.7, length=5e-3), inclination, azimuth=0.0),
        Jet(ox, Orifice.circular(d_ox, cd=0.8, length=5e-3), inclination, azimuth=math.pi),
    ], impingement_height=5e-3)


def two_feed_injector():
    film = ImpingingElement([
        Jet("film", Orifice.circular(0.5e-3, cd=0.7), math.radians(20), 0.0),
        Jet("film", Orifice.circular(0.5e-3, cd=0.7), math.radians(20), math.pi),
    ], impingement_height=3e-3)
    return Injector([
        ElementSet(doublet(), [0.01, 0.02, 0.02], [0.0, 0.0, math.pi], name="main"),
        ElementSet(film, np.full(8, 0.03), np.linspace(0, 2*np.pi, 8, endpoint=False),
                   orientation=math.pi/2, name="film"),
    ])


# --------------------------------------------------------------------------
# Elements
# --------------------------------------------------------------------------

def test_jets_are_aimed_at_the_impingement_point():
    element = doublet()
    holes = element.holes()
    offset = 5e-3*math.tan(math.radians(30))

    np.testing.assert_allclose(holes[0].exit, [0, offset, 0], atol=1e-15)
    np.testing.assert_allclose(holes[1].exit, [0, -offset, 0], atol=1e-15)
    for hole in holes:
        np.testing.assert_allclose(hole.target, [5e-3, 0, 0])
        np.testing.assert_allclose(np.linalg.norm(hole.direction), 1.0)
        # The axis passes through the target.
        along = hole.target - hole.exit
        np.testing.assert_allclose(np.cross(along, hole.direction), 0.0, atol=1e-15)
        assert math.degrees(math.acos(hole.direction[0])) == pytest.approx(30.0)
    assert element.free_jet_length(element.jets[0]) == pytest.approx(5e-3/math.cos(math.radians(30)))


def test_orifices_are_filtered_by_feed_with_one_entry_per_hole():
    like = ImpingingElement([Jet("f", Orifice.circular(1e-3, cd=0.7), 0.3),
                             Jet("f", Orifice.circular(1e-3, cd=0.7), 0.3, math.pi)],
                            impingement_height=3e-3)
    assert like.feeds == ("f",)
    assert len(like.orifices("f")) == 2
    assert like.orifices("other") == []
    assert doublet().feeds == ("fuel", "ox")


def test_feeds_can_be_hardware_objects_compared_by_identity():
    class Plenum:
        pass

    a, b = Plenum(), Plenum()
    element = doublet(fuel=a, ox=b)
    assert element.orifices(a)[0].hydraulic_diameter == 1e-3
    assert element.orifices(Plenum()) == []


@pytest.mark.parametrize("kwargs", [
    dict(inclination=math.pi/2),
    dict(inclination=-0.1),
    dict(inclination=0.1, azimuth=math.inf),
])
def test_invalid_jets_are_rejected(kwargs):
    with pytest.raises(ValueError):
        Jet("f", Orifice.circular(1e-3, cd=0.7), **kwargs)


def test_unhashable_feed_is_rejected():
    with pytest.raises(TypeError):
        Jet(["f"], Orifice.circular(1e-3, cd=0.7), 0.1)


def test_element_needs_jets_and_a_positive_height():
    with pytest.raises(ValueError):
        ImpingingElement([], impingement_height=1e-3)
    with pytest.raises(ValueError):
        ImpingingElement([Jet("f", Orifice.circular(1e-3, cd=0.7), 0.1)], impingement_height=0.0)


def test_diagnostics_of_an_unlike_doublet():
    element = doublet(d_fu=1e-3, d_ox=2e-3)
    flows = (JetFlow(mass_flow=0.01, velocity=20.0, density=800.0),
             JetFlow(mass_flow=0.03, velocity=10.0, density=1000.0))
    result = element.diagnostics(flows, frozenset({"ox"}))

    assert result.mixture_ratio == pytest.approx(3.0)
    assert result.momentum_flux_ratio == pytest.approx(1000*10**2/(800*20**2))
    assert result.rupe == pytest.approx(result.momentum_flux_ratio*2.0)
    # 0.2 N sideways one way, 0.3 N the other, both 30° off the axis.
    sideways = (0.3 - 0.2)*math.sin(math.radians(30))
    axial = (0.3 + 0.2)*math.cos(math.radians(30))
    assert result.resultant_angle == pytest.approx(math.atan2(sideways, axial))
    assert result.momentum == pytest.approx(math.hypot(sideways, axial))


def test_like_doublet_has_no_intra_element_mixing():
    like = ImpingingElement([Jet("f", Orifice.circular(1e-3, cd=0.7), 0.3),
                             Jet("f", Orifice.circular(1e-3, cd=0.7), 0.3, math.pi)],
                            impingement_height=3e-3)
    flow = JetFlow(0.01, 20.0, 800.0)
    result = like.diagnostics((flow, flow), frozenset({"ox"}))
    assert result.resultant_angle == pytest.approx(0.0, abs=1e-12)
    assert math.isnan(result.rupe)
    assert result.mixture_ratio == 0.0
    assert math.isnan(like.diagnostics((flow, flow), frozenset()).mixture_ratio)


# --------------------------------------------------------------------------
# Placement
# --------------------------------------------------------------------------

def test_element_sets_place_and_orient_their_elements():
    s = ElementSet(doublet(), [0.02], [math.pi/2], orientation=math.pi/2, name="a")
    holes = Injector([s]).holes()
    offset = 5e-3*math.tan(math.radians(30))

    # The centre is at +z. Orientation pi/2 turns the element's +y from the
    # local radial (+z) to the local tangential direction (-y).
    np.testing.assert_allclose(holes.exit[0], [0, -offset, 0.02], atol=1e-15)
    np.testing.assert_allclose(holes.exit[1], [0, offset, 0.02], atol=1e-15)
    np.testing.assert_allclose(holes.target[0], [5e-3, 0, 0.02], atol=1e-15)
    assert holes.set_name == ("a", "a")
    assert holes.feeds == ("fuel", "ox")
    assert holes.footprint[0] == pytest.approx(0.5e-3/math.cos(math.radians(30)))


@pytest.mark.parametrize("kwargs", [
    dict(r=[0.01, 0.02], theta=[0.0]),
    dict(r=[], theta=[]),
    dict(r=[-0.01], theta=[0.0]),
    dict(r=[0.01], theta=[0.0], name=""),
])
def test_invalid_element_sets_are_rejected(kwargs):
    kwargs.setdefault("name", "a")
    with pytest.raises(ValueError):
        ElementSet(doublet(), **kwargs)


def test_set_names_must_be_unique():
    s = ElementSet(doublet(), [0.01], [0.0], name="a")
    with pytest.raises(ValueError):
        Injector([s, s])


# --------------------------------------------------------------------------
# Hydraulics
# --------------------------------------------------------------------------

def test_equivalent_orifice_sums_every_hole_on_a_feed():
    injector = two_feed_injector()
    assert injector.feeds == ("fuel", "ox", "film")

    fuel = injector.equivalent_orifice("fuel")
    one = Orifice.circular(1e-3, cd=0.7)
    assert fuel.area == pytest.approx(3*one.area)
    assert fuel.effective_area == pytest.approx(3*one.effective_area)
    assert injector.hole_count("film") == 16
    assert injector.equivalent_orifice("film").effective_area == pytest.approx(
        16*Orifice.circular(0.5e-3, cd=0.7).effective_area)
    with pytest.raises(KeyError):
        injector.equivalent_orifice("nothing")


def test_each_feed_is_solved_on_its_own_pressure():
    injector = two_feed_injector()
    result = injector.solve({
        "fuel": FeedCondition(OIL, 300.0, pressure=11e5),
        "ox": FeedCondition(WATER, 300.0, pressure=12e5),
        "film": FeedCondition(OIL, 300.0, pressure=14e5),
    }, p_face=10e5, oxidizer=["ox"])

    for feed, fluid, p in (("fuel", OIL, 11e5), ("ox", WATER, 12e5), ("film", OIL, 14e5)):
        expected = injector.equivalent_orifice(feed).effective_area*math.sqrt(
            2*fluid.rho*(p - 10e5))
        assert result.feeds[feed].mass_flow == pytest.approx(expected, rel=1e-10)
        assert result.feeds[feed].pressure_drop == pytest.approx(p - 10e5)

    assert result.mixture_ratio == pytest.approx(
        result.feeds["ox"].mass_flow/(result.feeds["fuel"].mass_flow + result.feeds["film"].mass_flow))
    assert result.sets["film"].mass_flows == {"film": pytest.approx(result.feeds["film"].mass_flow)}
    assert result.sets["main"].count == 3


def test_flow_divides_between_holes_by_effective_area():
    big = ImpingingElement([Jet("f", Orifice.circular(2e-3, cd=0.8), 0.2)], impingement_height=4e-3)
    small = ImpingingElement([Jet("f", Orifice.circular(1e-3, cd=0.6), 0.2)], impingement_height=4e-3)
    injector = Injector([ElementSet(big, [0.01], [0.0], name="big"),
                         ElementSet(small, [0.02, 0.02], [0.0, math.pi], name="small")])
    result = injector.solve({"f": FeedCondition(WATER, 300.0, mass_flow=1.0)}, p_face=10e5)

    share = result.feeds["f"].set_mass_flows
    ratio = (0.8*4)/(2*0.6*1)
    assert share["big"]/share["small"] == pytest.approx(ratio)
    assert share["big"] + share["small"] == pytest.approx(1.0)
    assert math.isnan(result.mixture_ratio)


def test_pressure_and_flow_directions_are_inverse():
    injector = two_feed_injector()
    flows = {"fuel": 0.05, "ox": 0.08, "film": 0.02}
    fluids = {"fuel": OIL, "ox": WATER, "film": OIL}
    needed = injector.solve({f: FeedCondition(fluids[f], 300.0, mass_flow=m)
                             for f, m in flows.items()}, p_face=10e5)
    rated = injector.solve({f: FeedCondition(fluids[f], 300.0,
                                             pressure=needed.feeds[f].plenum_pressure)
                            for f in flows}, p_face=10e5)
    for feed, m in flows.items():
        assert rated.feeds[feed].mass_flow == pytest.approx(m, rel=1e-8)


def test_solve_feed_is_the_network_direction():
    injector = two_feed_injector()
    result = injector.solve_feed("ox", inlet_state=state(12e5), mass_flow=0.08, fluid=WATER)
    expected = (0.08/injector.equivalent_orifice("ox").effective_area)**2/2000.0
    assert result.pressure_drop == pytest.approx(expected, rel=1e-10)
    assert result.set_mass_flows == {"main": pytest.approx(0.08)}
    assert result.hole_count == 3


def test_real_fluid_design_round_trip():
    injector = Injector([ElementSet(doublet(), [0.01, 0.02], [0.0, 1.0], name="main")])
    conditions = {"fuel": FeedCondition("Methane", 120.0, mass_flow=0.05),
                  "ox": FeedCondition("Oxygen", 95.0, mass_flow=0.15)}
    needed = injector.solve(conditions, p_face=40e5, oxidizer=["ox"])
    rated = injector.solve({f: FeedCondition(c.fluid, c.temperature,
                                             pressure=needed.feeds[f].plenum_pressure)
                            for f, c in conditions.items()}, p_face=40e5, oxidizer=["ox"])
    assert rated.mixture_ratio == pytest.approx(3.0, rel=1e-8)


@pytest.mark.parametrize("conditions, oxidizer, match", [
    ({"fuel": FeedCondition(OIL, 300.0, pressure=11e5)}, (), "missing"),
    ({"fuel": FeedCondition(OIL, 300.0, pressure=11e5),
      "ox": FeedCondition(WATER, 300.0, pressure=11e5),
      "film": FeedCondition(OIL, 300.0, pressure=11e5),
      "extra": FeedCondition(OIL, 300.0, pressure=11e5)}, (), "unknown"),
    ({"fuel": FeedCondition(OIL, 300.0, pressure=11e5),
      "ox": FeedCondition(WATER, 300.0, pressure=11e5),
      "film": FeedCondition(OIL, 300.0, pressure=11e5)}, ("lox",), "oxidizer"),
])
def test_solve_validates_its_conditions(conditions, oxidizer, match):
    with pytest.raises(ValueError, match=match):
        two_feed_injector().solve(conditions, p_face=10e5, oxidizer=oxidizer)


def test_feed_condition_needs_exactly_one_of_pressure_and_flow():
    with pytest.raises(ValueError):
        FeedCondition(WATER, 300.0)
    with pytest.raises(ValueError):
        FeedCondition(WATER, 300.0, pressure=1e6, mass_flow=1.0)


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------

def test_checks_report_webs_reach_and_hole_proportions():
    injector = Injector([ElementSet(doublet(), [0.01, 0.02], [0.0, 0.0], name="main")])
    checks = injector.checks(chamber_radius=0.03)
    holes = injector.holes()

    face = holes.exit[:, 1:]
    webs = [np.linalg.norm(face[i] - face[j]) - holes.footprint[i] - holes.footprint[j]
            for i in range(4) for j in range(i + 1, 4)]
    assert checks.min_web == pytest.approx(min(webs))
    assert checks.overlaps == 0
    reach = np.linalg.norm(face, axis=1) + holes.footprint
    assert checks.outer_radius == pytest.approx(reach.max())
    assert checks.target_radius == pytest.approx(0.02)
    assert checks.wall_clearance == pytest.approx(0.03 - reach.max())
    assert checks.target_clearance == pytest.approx(0.01)
    assert checks.min_length_to_diameter == pytest.approx(5e-3/1.2e-3)
    assert math.isnan(injector.checks().wall_clearance)
    assert math.isnan(injector.checks().target_clearance)


def test_a_single_jet_aimed_at_the_wall():
    """A film jet is an impinging element of one jet whose target is on the wall."""
    wall, height, lean = 0.03, 20e-3, math.radians(15)
    jet = ImpingingElement([Jet("film", Orifice.circular(0.5e-3, cd=0.7), lean, azimuth=math.pi)],
                           impingement_height=height)
    injector = Injector([ElementSet(jet, np.full(12, wall), np.linspace(0, 2*np.pi, 12, endpoint=False),
                                    name="film")])
    holes = injector.holes()

    # Each hole sits inboard of the wall and leans outward onto it.
    radius = np.linalg.norm(holes.exit[:, 1:], axis=1)
    assert radius == pytest.approx(np.full(12, wall - height*math.tan(lean)))
    outward = np.einsum("ij,ij->i", holes.direction[:, 1:], holes.exit[:, 1:])/radius
    assert outward == pytest.approx(np.full(12, math.sin(lean)))

    checks = injector.checks(chamber_radius=wall)
    assert checks.target_clearance == pytest.approx(0.0, abs=1e-15)
    assert checks.wall_clearance > 0.0


def test_checks_count_overlapping_holes():
    crowded = Injector([ElementSet(doublet(), [0.010, 0.0105], [0.0, 0.0], name="main")])
    assert crowded.checks().overlaps > 0
    assert crowded.checks().min_web < 0.0
