"""Expansion thermodynamics of :class:`TurbineBlock`.

The block has to keep two enthalpy drops apart: the isentropic one,
``w / eta``, which sets the exit pressure, and the actual one, ``w``, which
sets the exit temperature. Collapsing them onto a single temperature drop
discards the loss that ``eta`` represents. The exit state is also taken from
CoolProp directly rather than from an ideal-gas ``T``-``p`` relation, so the
tests below pin both the energy bookkeeping and the real-gas path.
"""

import CoolProp.CoolProp as CP
import pytest

from pyskyfire.common.blocks import TurbineBlock
from pyskyfire.common.engine_network import Station


FLUID = "hydrogen"


def make_block(eta: float = 0.7353) -> TurbineBlock:
    return TurbineBlock(
        name="turbine",
        st_in="in",
        st_out="out",
        P_req_key="P_required",
        eta=eta,
        medium=FLUID,
    )


def run(block: TurbineBlock, station: Station, power: float):
    stations, signals = block.compute({"in": station}, {"P_required": power})
    return stations["out"], signals["dp_turbine"]


# RL10 fuel turbine, at the point the cycle balance converges to.
RL10_INLET = Station(5.194567e6, 250.0445, 2.746920)
RL10_POWER = 533711.6


def test_actual_enthalpy_drop_equals_the_extracted_work() -> None:
    inlet = RL10_INLET
    outlet, _ = run(make_block(), inlet, RL10_POWER)

    h_in = CP.PropsSI("Hmass", "T", inlet.T, "P", inlet.p, FLUID)
    h_out = CP.PropsSI("Hmass", "T", outlet.T, "P", outlet.p, FLUID)

    assert inlet.mdot * (h_in - h_out) == pytest.approx(RL10_POWER, rel=1e-9)


def test_exit_pressure_sits_on_the_isentropic_drop() -> None:
    inlet = RL10_INLET
    eta = 0.7353
    outlet, _ = run(make_block(eta), inlet, RL10_POWER)

    h_in = CP.PropsSI("Hmass", "T", inlet.T, "P", inlet.p, FLUID)
    s_in = CP.PropsSI("Smass", "T", inlet.T, "P", inlet.p, FLUID)
    h_out_s = CP.PropsSI("Hmass", "Smass", s_in, "P", outlet.p, FLUID)

    assert (h_in - h_out_s) == pytest.approx(RL10_POWER / inlet.mdot / eta, rel=1e-6)


def test_losses_reappear_as_exit_temperature() -> None:
    # A lossy turbine leaves the fluid hotter than a perfect one taking the
    # same shaft power, because the same pressure ratio must yield less work.
    inlet = RL10_INLET

    lossy, _ = run(make_block(0.70), inlet, RL10_POWER)
    ideal, _ = run(make_block(1.00), inlet, RL10_POWER)

    assert lossy.T > ideal.T
    # The ideal machine is isentropic; the real one is not.
    s_in = CP.PropsSI("Smass", "T", inlet.T, "P", inlet.p, FLUID)
    s_ideal = CP.PropsSI("Smass", "T", ideal.T, "P", ideal.p, FLUID)
    s_lossy = CP.PropsSI("Smass", "T", lossy.T, "P", lossy.p, FLUID)
    assert s_ideal == pytest.approx(s_in, rel=1e-6)
    assert s_lossy > s_in


def test_expansion_cools_and_drops_pressure() -> None:
    inlet = RL10_INLET
    outlet, dp = run(make_block(), inlet, RL10_POWER)

    assert outlet.T < inlet.T
    assert outlet.p < inlet.p
    assert dp == pytest.approx(inlet.p - outlet.p)
    assert outlet.mdot == inlet.mdot


def test_zero_power_is_a_pass_through() -> None:
    inlet = RL10_INLET
    outlet, dp = run(make_block(), inlet, 0.0)

    assert outlet.T == pytest.approx(inlet.T, rel=1e-9)
    assert outlet.p == pytest.approx(inlet.p, rel=1e-9)
    assert dp == pytest.approx(0.0, abs=1e-3)


def test_dense_hydrogen_departs_from_the_ideal_gas_pressure_ratio() -> None:
    # Cold, dense hydrogen is where gamma = cp/cv stops being the isentropic
    # exponent. The real-gas exit pressure must not agree with the ideal-gas
    # relation the block used to apply.
    inlet = Station(20e6, 150.0, 1.0)
    power = 300e3
    eta = 0.75

    outlet, _ = run(make_block(eta), inlet, power)

    c_p = CP.PropsSI("Cpmass", "T", inlet.T, "P", inlet.p, FLUID)
    c_v = CP.PropsSI("Cvmass", "T", inlet.T, "P", inlet.p, FLUID)
    gamma = c_p / c_v
    T_ideal = inlet.T - power / inlet.mdot / (eta * c_p)
    p_ideal = inlet.p * (T_ideal / inlet.T) ** (gamma / (gamma - 1.0))

    assert p_ideal / outlet.p > 1.05


def test_impossible_power_demand_is_rejected() -> None:
    with pytest.raises(ValueError, match="no physical expansion"):
        run(make_block(), RL10_INLET, 1e9)


def test_non_positive_mass_flow_is_rejected() -> None:
    with pytest.raises(ValueError, match="mdot must be positive"):
        run(make_block(), Station(5.19e6, 250.0, 0.0), RL10_POWER)
