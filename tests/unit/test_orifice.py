"""Flow contracts for orifices and orifice plates."""

import math
import warnings

import CoolProp.CoolProp as CP
import fluids
import pytest

from pyskyfire.ducts import (Orifice, OrificeCapacityError, OrificeCorrelationWarning,
                             OrificePlate, required_area)
from pyskyfire.thermo import FluidState, Incompressible, PHASE_UNKNOWN


def inlet(p, T):
    return FluidState(T=T, p=p, h=math.nan, quality=math.nan, rho=math.nan,
                      phase=PHASE_UNKNOWN, T_sat=math.nan)


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

def test_shapes_set_area_and_hydraulic_diameter():
    circle = Orifice.circular(2e-3, cd=0.7, length=6e-3)
    assert circle.area == pytest.approx(math.pi * 1e-6)
    assert circle.hydraulic_diameter == 2e-3
    assert circle.length_to_diameter == pytest.approx(3.0)
    assert circle.effective_area == pytest.approx(0.7 * circle.area)

    annulus = Orifice.annular(3e-3, 4e-3, cd=0.8)
    assert annulus.area == pytest.approx(math.pi * (16e-6 - 9e-6) / 4)
    assert annulus.hydraulic_diameter == pytest.approx(1e-3)

    slot = Orifice.slot(4e-3, 1e-3, cd=0.8)
    assert slot.area == pytest.approx(4e-6)
    assert slot.hydraulic_diameter == pytest.approx(1.6e-3)

    assert math.isnan(Orifice(1e-6, cd=1.0).length_to_diameter)


@pytest.mark.parametrize("kwargs", [
    dict(area=0.0, cd=0.7),
    dict(area=1e-6, cd=0.0),
    dict(area=1e-6, cd=1.2),
    dict(area=1e-6, cd=0.7, length=-1.0),
])
def test_invalid_orifice_is_rejected(kwargs):
    with pytest.raises(ValueError):
        Orifice(**kwargs)


def test_annulus_needs_outer_above_inner():
    with pytest.raises(ValueError):
        Orifice.annular(4e-3, 3e-3, cd=0.8)


# --------------------------------------------------------------------------
# Orifice flow
# --------------------------------------------------------------------------

def test_constant_density_gives_the_orifice_equation():
    orifice = Orifice.circular(2e-3, cd=0.7)
    result = orifice.solve(inlet_state=inlet(50e5, 300.0), mass_flow=0.05,
                           fluid=Incompressible(1000.0))

    expected = (0.05 / (0.7 * orifice.area))**2 / 2000.0
    assert result.pressure_loss == pytest.approx(expected, rel=1e-10)
    assert result.velocity == pytest.approx(math.sqrt(2 * expected / 1000.0))
    assert result.outlet_state.T == 300.0
    assert not result.choked


def test_water_is_close_to_incompressible_and_conserves_enthalpy():
    orifice = Orifice.circular(2e-3, cd=0.7)
    result = orifice.solve(inlet_state=inlet(50e5, 300.0), mass_flow=0.05, fluid="Water")

    rho = CP.PropsSI("D", "P", 50e5, "T", 300.0, "Water")
    assert result.pressure_loss == pytest.approx((0.05 / (0.7 * orifice.area))**2 / (2 * rho),
                                                 rel=1e-3)
    assert result.outlet_state.h == pytest.approx(result.inlet_state.h)


def test_gas_chokes_at_the_ideal_gas_critical_flux():
    orifice = Orifice.circular(2e-3, cd=0.7)
    p0, T0, gamma = 10e5, 300.0, 1.4
    R = 8314.46 / 28.0134
    ideal = (0.7 * orifice.area * p0 * math.sqrt(gamma / (R * T0))
             * (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1))))

    choked = orifice.mass_flow(inlet_state=inlet(p0, T0), p_out=1e5, fluid="Nitrogen")
    assert choked.choked
    assert choked.mass_flow == pytest.approx(ideal, rel=1e-2)
    assert choked.throat_state.p / p0 == pytest.approx((2 / (gamma + 1))**(gamma / (gamma - 1)),
                                                        rel=1e-2)
    assert choked.outlet_state.p == 1e5

    # Lower still, no more flow.
    lower = orifice.mass_flow(inlet_state=inlet(p0, T0), p_out=0.5e5, fluid="Nitrogen")
    assert lower.mass_flow == pytest.approx(choked.mass_flow, rel=1e-6)


def test_mass_flow_and_solve_are_inverse():
    orifice = Orifice.circular(2e-3, cd=0.7)
    forward = orifice.mass_flow(inlet_state=inlet(10e5, 300.0), p_out=8e5, fluid="Nitrogen")
    assert not forward.choked

    back = orifice.solve(inlet_state=inlet(10e5, 300.0), mass_flow=forward.mass_flow,
                         fluid="Nitrogen")
    assert back.outlet_state.p == pytest.approx(8e5, rel=1e-8)


def test_flow_beyond_choking_raises_with_capacity():
    orifice = Orifice.circular(2e-3, cd=0.7)
    choked = orifice.mass_flow(inlet_state=inlet(10e5, 300.0), p_out=1e5, fluid="Nitrogen")

    with pytest.raises(OrificeCapacityError) as info:
        orifice.solve(inlet_state=inlet(10e5, 300.0), mass_flow=1.01 * choked.mass_flow,
                      fluid="Nitrogen")
    assert info.value.max_mass_flow * 0.7 * orifice.area == pytest.approx(choked.mass_flow,
                                                                          rel=1e-6)


def test_liquid_cavitates_near_its_vapour_pressure():
    orifice = Orifice.circular(2e-3, cd=0.7)
    p1, T = 5e5, 90.0
    p_sat = CP.PropsSI("P", "T", T, "Q", 0, "Oxygen")
    rho = CP.PropsSI("D", "P", p1, "T", T, "Oxygen")

    above = orifice.mass_flow(inlet_state=inlet(p1, T), p_out=2 * p_sat, fluid="Oxygen")
    assert not above.choked

    flows = [orifice.mass_flow(inlet_state=inlet(p1, T), p_out=p, fluid="Oxygen")
             for p in (0.5 * p_sat, 0.2 * p_sat)]
    assert all(result.choked for result in flows)
    assert flows[0].mass_flow == pytest.approx(flows[1].mass_flow, rel=1e-6)
    assert flows[0].mass_flow == pytest.approx(
        0.7 * orifice.area * math.sqrt(2 * rho * (p1 - p_sat)), rel=5e-3)


def test_sized_orifice_gives_back_its_drop():
    state = inlet(50e5, 300.0)
    orifice = Orifice.sized(mass_flow=0.05, inlet_state=state, dp=5e5, fluid="Water", cd=0.7)
    result = orifice.solve(inlet_state=state, mass_flow=0.05, fluid="Water")
    assert result.pressure_loss == pytest.approx(5e5, rel=1e-8)


def test_required_area_refuses_a_choking_drop():
    with pytest.raises(ValueError, match="chokes"):
        required_area(mass_flow=0.01, inlet_state=inlet(10e5, 300.0), dp=8e5,
                      fluid="Nitrogen", cd=0.7)


def test_mass_flow_needs_a_lower_discharge_pressure():
    orifice = Orifice.circular(2e-3, cd=0.7)
    with pytest.raises(ValueError):
        orifice.mass_flow(inlet_state=inlet(10e5, 300.0), p_out=10e5, fluid="Nitrogen")


# --------------------------------------------------------------------------
# Orifice plate
# --------------------------------------------------------------------------

def test_plate_matches_iso_5167_for_a_liquid():
    D, d = 0.1, 0.05
    plate = OrificePlate(D, d, taps="flange")
    rho = CP.PropsSI("D", "P", 10e5, "T", 300.0, "Water")
    mu = CP.PropsSI("V", "P", 10e5, "T", 300.0, "Water")

    result = plate.solve(inlet_state=inlet(10e5, 300.0), mass_flow=20.0, fluid="Water")

    p2 = fluids.flow_meter.differential_pressure_meter_solver(
        D=D, D2=d, P1=10e5, m=20.0, rho=rho, mu=mu, meter_type="ISO 5167 orifice",
        taps="flange", epsilon_specified=1.0)
    assert result.differential_pressure == pytest.approx(10e5 - p2, rel=1e-4)
    assert result.pressure_loss == pytest.approx(
        fluids.flow_meter.dP_orifice(D, d, 10e5, p2, result.discharge_coefficient), rel=1e-4)
    assert result.pressure_loss < result.differential_pressure
    assert result.beta == 0.5
    assert result.reynolds == pytest.approx(4 * 20.0 / (math.pi * D * mu))


@pytest.mark.parametrize("mass_flow", [0.3, 0.9, 2.0])
def test_plate_matches_iso_5167_expansibility_for_a_gas(mass_flow):
    D, d = 0.1, 0.05
    plate = OrificePlate(D, d, taps="flange")
    props = {name: CP.PropsSI(key, "P", 10e5, "T", 300.0, "Nitrogen")
             for name, key in (("rho", "D"), ("mu", "V"),
                               ("k", "isentropic_expansion_coefficient"))}

    result = plate.solve(inlet_state=inlet(10e5, 300.0), mass_flow=mass_flow, fluid="Nitrogen")

    p2 = fluids.flow_meter.differential_pressure_meter_solver(
        D=D, D2=d, P1=10e5, m=mass_flow, meter_type="ISO 5167 orifice", taps="flange", **props)
    assert result.differential_pressure == pytest.approx(10e5 - p2, rel=1e-4)


def test_plate_uses_a_supplied_coefficient():
    plate = OrificePlate(0.1, 0.05, cd=0.62)
    result = plate.solve(inlet_state=inlet(10e5, 300.0), mass_flow=20.0,
                         fluid=Incompressible(1000.0))

    b4 = 0.5**4
    expected = (20.0 / (0.62 * plate.area))**2 * (1 - b4) / 2000.0
    assert result.discharge_coefficient == 0.62
    assert result.differential_pressure == pytest.approx(expected, rel=1e-10)
    assert math.isnan(result.reynolds)


def test_correlation_needs_a_viscosity_without_coolprop():
    plate = OrificePlate(0.1, 0.05)
    with pytest.raises(ValueError, match="viscosity"):
        plate.solve(inlet_state=inlet(10e5, 300.0), mass_flow=20.0,
                    fluid=Incompressible(1000.0))

    result = plate.solve(inlet_state=inlet(10e5, 300.0), mass_flow=20.0,
                         fluid=Incompressible(1000.0), viscosity=1e-3)
    assert 0.59 < result.discharge_coefficient < 0.62


def test_small_plate_warns_that_the_correlation_is_extrapolated():
    plate = OrificePlate(0.02, 0.01)
    with pytest.warns(OrificeCorrelationWarning, match="pipe diameter"):
        plate.solve(inlet_state=inlet(10e5, 300.0), mass_flow=1.0, fluid="Water")


def test_plate_in_range_does_not_warn():
    plate = OrificePlate(0.1, 0.05)
    with warnings.catch_warnings():
        warnings.simplefilter("error", OrificeCorrelationWarning)
        plate.solve(inlet_state=inlet(10e5, 300.0), mass_flow=20.0, fluid="Water")


@pytest.mark.parametrize("args, kwargs", [
    ((0.05, 0.05), {}),
    ((0.1, 0.05), dict(taps="pipe")),
    ((0.1, 0.05), dict(cd=1.5)),
])
def test_invalid_plate_is_rejected(args, kwargs):
    with pytest.raises(ValueError):
        OrificePlate(*args, **kwargs)
