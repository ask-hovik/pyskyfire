"""Vaned diffuser: Gülich's sizing, the channel-first vanes, losses and checks."""
import numpy as np
import pytest

from pyskyfire.pump import Impeller, VanedDiffuser, minimum_gap_ratio


class Water:
    fluid = "Water"

    def get_rho(self, T, p):
        return 1000.

    def get_mu(self, T, p):
        return 1e-3


def state(p=50e5):
    from pyskyfire.thermo.state import FluidState, PHASE_UNKNOWN
    nan = float("nan")
    return FluidState(T=300., p=p, h=nan, quality=nan, rho=nan, phase=PHASE_UNKNOWN, T_sat=nan)


@pytest.fixture(scope="module")
def impeller():
    return Impeller(Q=.03, H=400., n=6000., shaft_diameter=.03, inlet_hub_diameter=.04,
                    blade_count=6, blade_thickness=.003, hydraulic_efficiency=.95, rho=1000.)


def design(impeller, **changes):
    return VanedDiffuser.from_impeller(impeller, **{"vane_count": 8, **changes})


def test_inlet_follows_table_7_7_and_the_gap_rule(impeller):
    d = design(impeller)
    dims, mean = impeller.dimensions, impeller.velocity_triangles.outlet[1]
    least = minimum_gap_ratio(nq=dims.nq, rho=1000., H=impeller.H)
    assert least == pytest.approx(1.015+0.08*(1000*400/1e6-0.1)**0.8)
    assert d.d3 == pytest.approx(least*dims.d2)
    assert d.b3 == pytest.approx(1.15*dims.b2)
    # Free vortex across the gap (T7.7.3), continuity at d3 (T7.7.4).
    assert d.c3u == pytest.approx(mean.cu*dims.d2/d.d3)
    assert np.tan(d.flow_angle) == pytest.approx(dims.internal_flow/(np.pi*d.d3*d.b3)/d.c3u)
    # Default throat decelerates to c3q/c2 = 0.8 (T7.7.11).
    assert d.checks.deceleration == pytest.approx(.8)
    X = np.exp(dims.internal_flow/(d.b3*mean.cu*dims.d2/2*8))
    assert d.throat_factor == pytest.approx(d.a3/(d.d3/2*(X-1)))
    # Or from f_a3 (T7.7.10).
    assert design(impeller, throat_factor=1.2).throat_factor == pytest.approx(1.2)


def test_channel_is_an_optimum_diffuser_for_its_length(impeller):
    d = design(impeller)
    R = np.sqrt(d.a3*d.b3/np.pi)
    assert d.channel_length == pytest.approx(4*d.a3)
    assert d.area_ratio == pytest.approx(1.05+0.184*d.channel_length/R)
    assert d.a4*d.b4 == pytest.approx(d.area_ratio*d.a3*d.b3)
    assert d.checks.opening_angle_deg < d.checks.allowed_opening_angle_deg
    # A flared wall opens the channel as well, and a fixed outer diameter
    # sets the length instead.
    flared = design(impeller, outlet_width_ratio=1.2)
    assert flared.b4 == pytest.approx(1.2*flared.b3)
    fixed = design(impeller, outer_diameter_ratio=1.45)
    assert fixed.d4 == pytest.approx(1.45*impeller.dimensions.d2)


def test_vanes_are_thin_and_draw_the_design_channel(impeller):
    d = design(impeller)
    assert d.vane_thickness == pytest.approx(.012*impeller.dimensions.d2)
    # The drawn throat and exit are the design ones, and the leading spiral
    # sits on the approach flow within Gülich's +-3 degrees.
    assert d.checks.drawn_throat_ratio == pytest.approx(1, abs=1e-3)
    assert d.checks.drawn_exit_ratio == pytest.approx(1, abs=1e-3)
    assert abs(d.checks.incidence_deg) < 3
    assert 0 < d.checks.deflection_deg < 25
    # Recovery is the optimum for the length, corrected for curvature (Eq. 1.46).
    R = np.sqrt(d.a3*d.b3/np.pi)
    assert d.pressure_recovery == pytest.approx(
        .36*(d.channel_length/R)**.26*(1-.13*d.checks.deflection_deg/30))
    # Every face point is within half a vane of the vane's own camber, so
    # nothing thickens into a wedge.
    inner, outer = d._faces
    assert np.max(np.linalg.norm(inner-outer, axis=1)) <= d.vane_thickness*1.01
    g = d.build_geometry(chord_points=41)
    assert g.vane.suction.shape == (2, 41, 3)
    # Walls symmetric about the impeller outlet mid-plane; hub at +x.
    np.testing.assert_allclose(g.hub[:, 0], -g.shroud[:, 0])
    assert g.hub[-1, 1] == pytest.approx(d.d4/2)
    assert np.all(g.vane.suction[..., 1] >= d.d3/2-d.vane_thickness/2)
    assert np.all(g.vane.pressure[..., 1] <= d.d4/2+d.vane_thickness/2)
    from pyskyfire.viz.components import _impeller_blade_mesh
    faces = _impeller_blade_mesh(g.vane).triangles
    edges = np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    _, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    assert np.all(counts == 2)


def test_losses_follow_table_3_8(impeller):
    d = design(impeller)
    mdot = 1000.*impeller.dimensions.internal_flow
    r = d.solve(inlet_state=state(), mass_flow=mdot, transport=Water())
    assert r.zeta_total == pytest.approx(sum(r.loss_breakdown().values()))
    assert r.pressure_loss == pytest.approx(r.zeta_total*1000.*d.u2**2/2)
    assert r.states[-1].p == pytest.approx(50e5-r.pressure_loss)
    head = (r.c3q/d.u2)**2
    assert r.zeta_shock == pytest.approx(head*.3*(r.c2/r.c3q-1)**2)
    assert r.zeta_channel == pytest.approx(head*(1-d.pressure_recovery-1/d.area_ratio**2))
    # The outlet handed to the volute conserves mass over the whole ring.
    assert d.outlet.cm*2*np.pi*d.outlet.radius*d.outlet.width == pytest.approx(impeller.dimensions.internal_flow)
    assert 0 < d.outlet.cu < r.c4


def test_checks_flag_what_gulich_warns_against(impeller):
    from pyskyfire.pump import vane_count_conflicts
    assert design(impeller).checks.vane_count_conflicts == vane_count_conflicts(8, 6) == ()
    # Six vanes against six blades coincide; seven differ by one; nine put
    # the third blade harmonic on the second vane harmonic.
    assert (1, 1, 0) in vane_count_conflicts(6, 6)
    assert (1, 1, 1) in vane_count_conflicts(7, 6)
    assert (3, 2, 0) in vane_count_conflicts(9, 6)


@pytest.mark.parametrize("changes, message", [
    (dict(throat_factor=1.2, deceleration=.8), "throat"),
    (dict(outer_diameter_ratio=1.4, length_ratio=4), "outer diameter"),
    (dict(gap_ratio=1.), "exceed 1"),
    (dict(outer_diameter_ratio=1.1, gap_ratio=1.1), "no room"),
    (dict(vane_count=2), "vane_count"),
    (dict(vane_thickness_ratio=.2), "thick"),
])
def test_bad_inputs_are_refused(impeller, changes, message):
    with pytest.raises(ValueError, match=message):
        design(impeller, **changes)


def test_network_block_emits_the_loss(impeller):
    from pyskyfire.network import DiffuserBlock, Station
    d = design(impeller)
    block = DiffuserBlock("diffuser", "in", "out", d, Water())
    mdot = 1000.*impeller.dimensions.internal_flow
    out, signals = block.compute({"in": Station(50e5, 300., mdot)}, {})
    loss = d.solve(inlet_state=state(), mass_flow=mdot, transport=Water()).pressure_loss
    assert signals["dp_diffuser"] == pytest.approx(loss)
    assert out["out"].p == pytest.approx(50e5-loss)
    assert out["out"].T == 300.
    with pytest.raises(TypeError):
        DiffuserBlock("x", "in", "out", object(), Water())
