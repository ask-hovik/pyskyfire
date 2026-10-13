"""Geometry and pressure-loss contracts for annular ducts."""

import math

import fluids
import numpy as np
import pytest

from pyskyfire import ducts, viz
from pyskyfire.common.friction import annulus_laminar_factor, darcy_friction_factor
from pyskyfire.thermo import FluidState, PHASE_UNKNOWN


class ConstantTransport:
    def __init__(self, rho=1000.0, mu=1.0e-3):
        self._rho = rho
        self._mu = mu

    def get_rho(self, T, p):
        return self._rho

    def get_mu(self, T, p):
        return self._mu

    def get_cp(self, T, p):
        return 4180.0


def inlet(p=2.0e6, T=300.0):
    return FluidState(T=T, p=p, h=np.nan, quality=np.nan, rho=np.nan,
                      phase=PHASE_UNKNOWN, T_sat=np.nan)


def axial_duct(r=0.05, height=0.01, length=0.5, stations=21):
    x = np.linspace(0.0, length, stations)
    return ducts.AnnularDuct(np.column_stack((x, np.full_like(x, r))), height)


def crossover(R=0.015, height=0.01, r_in=0.06, r_top=0.10):
    """Radially out, a 180 degree meridional turn, radially back in."""
    turn = np.linspace(0.0, np.pi, 25)
    out = np.column_stack((np.zeros(10), np.linspace(r_in, r_top, 10)))
    arc = np.column_stack((R - R*np.cos(turn), r_top + R*np.sin(turn)))
    back = np.column_stack((np.full(10, 2*R), np.linspace(r_top, r_in, 10)))
    return ducts.AnnularDuct(np.vstack((out[:-1], arc, back[1:])), height)


# -- friction ----------------------------------------------------------------

def test_annulus_laminar_factor_spans_pipe_to_parallel_plates():
    assert annulus_laminar_factor(0.0) == 1.0
    assert annulus_laminar_factor(1.0) == 1.5
    # Shah & London tabulate a Fanning f*Re of 23.81 at a radius ratio of 0.5.
    assert 64*annulus_laminar_factor(0.5) == pytest.approx(4*23.81, rel=1e-3)
    values = [annulus_laminar_factor(k) for k in np.linspace(0, 1, 101)]
    assert np.all(np.diff(values) >= -1e-9)
    with pytest.raises(ValueError):
        annulus_laminar_factor(1.2)


@pytest.mark.parametrize("model", ["native", "fluids"])
def test_laminar_factor_scales_only_the_laminar_branch(model):
    assert darcy_friction_factor(500, 0.01, model=model, laminar_factor=1.5) == pytest.approx(
        1.5*darcy_friction_factor(500, 0.01, model=model))
    assert darcy_friction_factor(1e5, 0.01, model=model, laminar_factor=1.5) == pytest.approx(
        darcy_friction_factor(1e5, 0.01, model=model))


# -- symmetric about the shaft ----------------------------------------------

def test_straight_axial_duct_matches_darcy_weisbach():
    duct = axial_duct()
    area = math.pi*(0.055**2 - 0.045**2)
    assert duct.area == pytest.approx(area)
    assert duct.hydraulic_diameter == pytest.approx(0.02)
    result = duct.solve(inlet_state=inlet(), mass_flow=5.0, transport=ConstantTransport())
    velocity = 5.0/(1000.0*area)
    reynolds = 1000.0*velocity*0.02/1e-3
    expected = fluids.friction.Clamond(reynolds, 0.0)*0.5/0.02*1000.0*velocity**2/2
    assert result.total_pressure_loss == pytest.approx(expected, rel=1e-9)
    np.testing.assert_allclose(result.p_static,
                               result.p_stagnation - 0.5*result.density*result.velocity**2)


def test_span_lines_are_normal_to_the_mean_line():
    duct = crossover()
    span = duct.shroud - duct.hub
    np.testing.assert_allclose(np.einsum("ij,ij->i", span, duct._tangents), 0.0, atol=1e-12)
    np.testing.assert_allclose(duct.height, 0.01)
    # Radial legs sweep a disc-like band, 2*pi*r*b; the top of the turn an annulus.
    assert duct.area[0] == pytest.approx(2*np.pi*0.06*0.01)


def test_swirl_lengthens_the_friction_path_and_decays():
    duct = axial_duct()
    plain = duct.solve(inlet_state=inlet(), mass_flow=5.0, transport=ConstantTransport())
    velocity = plain.meridional_velocity[0]
    swirling = duct.solve(inlet_state=inlet(), mass_flow=5.0, transport=ConstantTransport(),
                          inlet_swirl=velocity)
    assert swirling.total_pressure_loss > 2.0*plain.total_pressure_loss
    assert swirling.swirl_velocity[0] == pytest.approx(velocity)
    assert np.all(np.diff(swirling.swirl_velocity) < 0.0)
    assert swirling.swirl_angle[0] == pytest.approx(np.pi/4)


def test_swirl_keeps_its_angular_momentum_through_a_change_of_radius():
    r = np.linspace(0.10, 0.05, 11)
    duct = ducts.AnnularDuct(np.column_stack((np.zeros_like(r), r)), 0.01)
    # Enough through-flow that the swirl travels a short helix and friction
    # takes only a little of its angular momentum.
    result = duct.solve(inlet_state=inlet(), mass_flow=50.0,
                        transport=ConstantTransport(), inlet_swirl=10.0)
    angular_momentum = result.swirl_velocity*duct.r_mean
    assert np.all(np.diff(angular_momentum) < 0.0)
    assert angular_momentum[-1] == pytest.approx(angular_momentum[0], rel=0.03)
    assert result.swirl_velocity[-1] == pytest.approx(20.0, rel=0.03)


def test_from_walls_recovers_normal_pairing_from_skewed_walls():
    hub = np.column_stack((np.linspace(0, 0.2, 7), np.linspace(0.04, 0.06, 7)))
    shroud = np.column_stack((np.linspace(0, 0.2, 13), np.linspace(0.07, 0.09, 13)))
    duct = ducts.AnnularDuct.from_walls(hub, shroud, stations=11)
    np.testing.assert_allclose(duct.height[1:-1], 0.03*np.cos(np.arctan(0.1)), rtol=1e-9)
    # The ends are the joining planes and stay as given.
    np.testing.assert_allclose(duct.hub[[0, -1]], hub[[0, -1]])
    np.testing.assert_allclose(duct.shroud[[0, -1]], shroud[[0, -1]])


def test_from_walls_reproduces_a_crossover_turn():
    original = crossover()
    duct = ducts.AnnularDuct.from_walls(original.hub, original.shroud, stations=41)
    np.testing.assert_allclose(duct.height, 0.01, atol=1e-5)
    assert duct.s[-1] == pytest.approx(original.s[-1], rel=1e-3)


def test_sharp_meridional_turn_warns():
    with pytest.warns(ducts.PipeCorrelationWarning, match="tighter than one hydraulic"):
        crossover().solve(inlet_state=inlet(), mass_flow=5.0, transport=ConstantTransport())


def test_rejects_backward_and_folding_ducts():
    x = np.linspace(0.0, -0.2, 5)
    with pytest.raises(ValueError, match="backward"):
        ducts.AnnularDuct(np.column_stack((x, np.full_like(x, 0.05))), 0.01)
    with pytest.raises(ValueError, match="folds"):
        crossover(R=0.004, height=0.01)
    with pytest.raises(ValueError, match="height"):
        axial_duct(height=-0.01)


def test_end_is_the_passage_annulus_end():
    end = axial_duct().end("outlet")
    assert end.r_hub == pytest.approx(0.045)
    assert end.r_shroud == pytest.approx(0.055)


def test_widening_annulus_loses_the_diffuser_term_of_its_equivalent_cone():
    # A radial annulus of constant height widens its area as 2*pi*r*b.
    r = np.linspace(0.05, 0.10, 21)
    duct = ducts.AnnularDuct(np.column_stack((np.zeros_like(r), r)), 0.01)
    area = duct.area
    D = np.sqrt(4*area/np.pi)
    mid = duct._mid_area()
    K = float(np.sum(duct.area_change_loss_coefficient*(area[0]/mid)**2))
    angles = np.degrees(duct.cone_angle)
    weight = area[:-1]**-2 - area[1:]**-2
    expected = sum(w*fluids.fittings.diffuser_conical(D[0], D[-1], angle=a, fd=0.0)
                   for w, a in zip(weight/weight.sum(), angles))

    assert np.all(duct.cone_angle > 0)
    assert K == pytest.approx(expected, rel=1e-12)
    result = duct.solve(inlet_state=inlet(), mass_flow=5.0, transport=ConstantTransport())
    np.testing.assert_allclose(result.loss_coefficient - result.darcy_factor*duct.lengths/(
        0.5*(duct.hydraulic_diameter[:-1] + duct.hydraulic_diameter[1:])),
        duct.area_change_loss_coefficient, rtol=1e-9, atol=1e-15)


def test_constant_area_annulus_has_no_area_change_loss():
    assert not np.any(axial_duct().area_change_loss_coefficient)


# -- display ----------------------------------------------------------------

def test_component_draws_both_walls_with_aligned_fields():
    duct = axial_duct(stations=11)
    result = duct.solve(inlet_state=inlet(), mass_flow=2.0, transport=ConstantTransport())
    node = viz.AnnularDuctComponent(duct, result, n_theta=12).build()
    assert [c.name for c in node.children] == ["hub", "shroud"]
    hub = node.children[0]
    assert len(hub.fields["p_static"].values) == len(hub.mesh.points)
    assert len(hub.fields["pressure_loss"].values) == len(hub.mesh.triangles)
    assert "swirl_velocity" in hub.fields
    assert set(node.anchors) == {"inlet", "outlet"}
    with pytest.raises(TypeError, match="AnnularDuct"):
        viz.AnnularDuctComponent(object()).build()
