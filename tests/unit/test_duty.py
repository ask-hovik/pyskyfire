"""Equations of state and the pump and turbine duties built on them.

The duties are cross-checked against the network blocks that already carry
the same thermodynamics: a pump block's ideal work is ``Δp/rho`` at inlet
density, which is exactly an incompressible duty's head, and a turbine
block's exit pressure is where the isentropic drop equals ``w/eta``, which a
duty built across the same two pressures has to reproduce.
"""

import math
import pickle

import CoolProp.CoolProp as CP
import numpy as np
import pytest

from pyskyfire.network.blocks import TurbineBlock
from pyskyfire.network.engine_network import Station
from pyskyfire.pump import Impeller, PumpDuty
from pyskyfire.pump.impeller import G0, specific_speed
from pyskyfire.thermo import (
    CoolPropEOS,
    Fluid,
    FluidTransport,
    Incompressible,
    MachineDuty,
    PHASE_LIQUID,
    as_equation_of_state,
)
from pyskyfire.turbine import TurbineDuty

IMPELLER = dict(shaft_diameter=.02, inlet_hub_diameter=.03, blade_count=6,
                beta2_deg=25, blade_thickness=.0015, hydraulic_efficiency=.85)


# -- equations of state -------------------------------------------------------

def test_incompressible_isentropic_head_is_pressure_over_density():
    eos = Incompressible(800.)
    inlet = eos.state_pT(2e5, 300.)
    outlet = eos.isentropic(inlet, 50e5)
    assert outlet.h - inlet.h == pytest.approx(48e5/800.)
    assert outlet.T == inlet.T and outlet.rho == 800.


def test_incompressible_temperature_needs_cp_to_come_back_from_enthalpy():
    assert math.isnan(Incompressible(800.).state_ph(1e5, 1e3).T)
    eos = Incompressible(800., cp=2000.)
    state = eos.state_pT(3e5, 310.)
    assert eos.state_ph(state.p, state.h).T == pytest.approx(310.)


@pytest.mark.parametrize("rho, cp", [(0., None), (np.nan, None), (1., -1.)])
def test_incompressible_rejects_bad_properties(rho, cp):
    with pytest.raises(ValueError):
        Incompressible(rho, cp)


def test_coolprop_states_match_coolprop():
    eos = CoolPropEOS("Nitrogen")
    state = eos.state_pT(30e5, 400.)
    assert state.rho == pytest.approx(CP.PropsSI("D", "P", 30e5, "T", 400., "Nitrogen"))
    assert state.h == pytest.approx(CP.PropsSI("H", "P", 30e5, "T", 400., "Nitrogen"))
    assert eos.state_ph(state.p, state.h).T == pytest.approx(400.)


def test_coolprop_isentropic_keeps_entropy_and_approaches_ideal_gas():
    eos = CoolPropEOS("Nitrogen")
    inlet = eos.state_pT(2e5, 300.)
    outlet = eos.isentropic(inlet, 1e5)
    assert eos.entropy(outlet) == pytest.approx(eos.entropy(inlet), rel=1e-9)
    gamma = 1.4
    assert outlet.T/inlet.T == pytest.approx(0.5**((gamma - 1)/gamma), rel=2e-3)


def test_coolprop_reads_a_saturation_pair_as_liquid():
    T_sat = CP.PropsSI("T", "P", 5e5, "Q", 0, "Methane")
    state = CoolPropEOS("Methane").state_pT(5e5, T_sat)
    assert state.phase == PHASE_LIQUID
    assert state.rho == pytest.approx(CP.PropsSI("D", "P", 5e5, "Q", 0, "Methane"), rel=1e-4)


def test_coolprop_accepts_every_way_of_naming_a_fluid_and_pickles():
    names = ["Methane", Fluid("fuel", ["Methane"], [1.0]),
             FluidTransport(Fluid("fuel", ["Methane"], [1.0]))]
    eoses = [as_equation_of_state(name) for name in names]
    assert all(isinstance(e, CoolPropEOS) for e in eoses)
    assert len({e.state_pT(50e5, 300.).rho for e in eoses}) == 1
    eos = eoses[0]
    eos.state_pT(50e5, 300.)
    assert pickle.loads(pickle.dumps(eos)).state_pT(50e5, 300.).rho == eos.state_pT(50e5, 300.).rho
    assert as_equation_of_state(eos) is eos
    with pytest.raises(TypeError):
        CoolPropEOS(3)


# -- pump duty -----------------------------------------------------------------

def test_pump_head_duty_without_density_has_no_state():
    duty = PumpDuty.from_head(.01, 100., 10000.)
    assert duty.H == pytest.approx(100.) and duty.Q == .01
    assert duty.mdot is None and duty.inlet is None and duty.rho_in is None
    assert duty.specific_speed == pytest.approx(specific_speed(10000., .01, 100.))


def test_pump_pressure_rise_matches_impeller_and_incompressible_states():
    direct = PumpDuty.from_pressure_rise(mdot=10., rho=1000., dp=980665., n=10000.)
    assert direct.Q == pytest.approx(.01) and direct.H == pytest.approx(100.)
    states = PumpDuty.from_states(fluid=Incompressible(1000.), mdot=10., p_in=2e5, T_in=300.,
                                  p_out=2e5 + 980665., n=10000.)
    assert states.H == pytest.approx(direct.H) and states.Q_out == pytest.approx(direct.Q)
    impeller = Impeller.from_pressure_rise(mdot=10., rho=1000., dp=980665., n=10000., **IMPELLER)
    assert Impeller.from_duty(direct, **IMPELLER).dimensions == impeller.dimensions
    assert Impeller.from_duty(direct, **IMPELLER).rho == 1000.


def test_real_liquid_pump_head_differs_from_incompressible_by_its_compression():
    T_sat = CP.PropsSI("T", "P", 5e5, "Q", 0, "Methane")
    duty = PumpDuty.from_states(fluid="Methane", mdot=20., p_in=5e5, T_in=T_sat,
                                p_out=230e5, n=40000.)
    incompressible = (230e5 - 5e5)/duty.rho_in
    # Methane stiffens by a few percent over 225 bar, so the ideal work is a
    # little less than Δp over inlet density, and the flow leaves smaller.
    assert 0.95 < duty.isentropic_head/incompressible < 1.0
    assert duty.Q_out < duty.Q_in
    assert duty.H == pytest.approx(duty.isentropic_head/G0)


def test_pump_station_duty_reads_inlet_flow_and_ignores_outlet_temperature():
    eos = Incompressible(1000.)
    a = PumpDuty.from_stations(Station(2e5, 300., 10.), Station(50e5, 300., 10.), fluid=eos, n=1e4)
    b = PumpDuty.from_stations(Station(2e5, 300., 10.), Station(50e5, 999., 3.), fluid=eos, n=1e4)
    assert a.mdot == b.mdot == 10.
    assert (a.isentropic_head, a.Q_in, a.Q_out) == (b.isentropic_head, b.Q_in, b.Q_out)


def test_pump_must_raise_pressure():
    with pytest.raises(ValueError, match="above its inlet"):
        PumpDuty.from_states(fluid=Incompressible(1000.), mdot=1., p_in=5e5, T_in=300., p_out=5e5, n=1e4)


def test_duty_fields_must_agree():
    eos = Incompressible(1000.)
    good = PumpDuty.from_states(fluid=eos, mdot=10., p_in=2e5, T_in=300., p_out=5e5, n=1e4)
    with pytest.raises(ValueError, match="Q_in"):
        MachineDuty(n=1e4, isentropic_head=good.isentropic_head, Q_in=.5, Q_out=good.Q_out,
                    mdot=10., inlet=good.inlet, outlet=good.outlet)
    with pytest.raises(ValueError, match="together"):
        MachineDuty(n=1e4, isentropic_head=1., Q_in=1., Q_out=1., mdot=1., inlet=good.inlet)
    with pytest.raises(ValueError, match="n must be"):
        PumpDuty.from_head(.01, 100., 0.)


# -- turbine duty --------------------------------------------------------------

def test_turbine_duty_reproduces_the_turbine_blocks_expansion():
    eta, P, mdot = .7353, 2.0e6, 5.
    block = TurbineBlock("turbine", "in", "out", "P_required", eta, "hydrogen")
    inlet = Station(80e5, 250., mdot)
    out, _ = block.compute({"in": inlet}, {"P_required": P})
    duty = TurbineDuty.from_stations(inlet, out["out"], fluid="hydrogen", n=60000.)
    # The block set the exit pressure so that the isentropic drop is w/eta.
    assert duty.isentropic_head == pytest.approx(P/(mdot*eta), rel=1e-6)
    assert duty.outlet_basis == "total"
    assert duty.Q_out > duty.Q_in


def test_turbine_similarity_numbers():
    duty = TurbineDuty.from_states(fluid="Nitrogen", mdot=4., p0_in=40e5, T0_in=700.,
                                   p_out=20e5, n=30000.)
    c0 = math.sqrt(2*duty.isentropic_head)
    omega = 30000*2*math.pi/60
    assert duty.spouting_velocity == pytest.approx(c0)
    assert duty.pressure_ratio == pytest.approx(2.)
    assert duty.blade_jet_ratio(.2) == pytest.approx(omega*.1/c0)
    assert duty.specific_speed == pytest.approx(omega*math.sqrt(duty.Q_out)/duty.isentropic_head**.75)
    assert duty.specific_diameter(.2) == pytest.approx(.2*duty.isentropic_head**.25/math.sqrt(duty.Q_out))
    s0 = CP.PropsSI("S", "P", 40e5, "T", 700., "Nitrogen")
    h0 = CP.PropsSI("H", "P", 40e5, "T", 700., "Nitrogen")
    assert duty.isentropic_head == pytest.approx(h0 - CP.PropsSI("H", "P", 20e5, "S", s0, "Nitrogen"))


def test_turbine_duty_needs_states_expansion_and_a_known_basis():
    with pytest.raises(ValueError, match="needs its states"):
        TurbineDuty(n=1e4, isentropic_head=1e5, Q_in=1., Q_out=1.)
    with pytest.raises(ValueError, match="below its inlet"):
        TurbineDuty.from_states(fluid="Nitrogen", mdot=1., p0_in=10e5, T0_in=500., p_out=11e5, n=1e4)
    with pytest.raises(ValueError, match="outlet_basis"):
        TurbineDuty.from_states(fluid="Nitrogen", mdot=1., p0_in=10e5, T0_in=500., p_out=5e5,
                                n=1e4, outlet_basis="dynamic")
