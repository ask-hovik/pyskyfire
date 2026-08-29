"""Persistence tests for objects backed by native extension handles."""

import cloudpickle
import numpy as np
import pytest

import pyskyfire as psf
from pyskyfire.skycea.aerothermodynamics import GasState


def test_aerothermodynamics_rebuilds_cea_solvers_after_round_trip() -> None:
    fuel = psf.common.Fluid("fuel", ["CH4"], [1.0])
    oxidizer = psf.common.Fluid("oxidizer", ["O2"], [1.0])
    transport = psf.skycea.Aerothermodynamics.from_F_eps_Lstar(
        fu=fuel,
        ox=oxidizer,
        MR=2.8,
        p_c=1.0e7,
        F=1.0e5,
        eps=10.0,
        L_star=1.1,
        T_fu_in=300.0,
        T_ox_in=300.0,
    )

    expected_state = transport.equilibrium_state("tp", 3000.0, 2.0e6)
    transport._cache[(0.0, None, None, transport.transport)] = GasState(T=1.0)
    restored = cloudpickle.loads(cloudpickle.dumps(transport))
    assert restored.cache_size == 0
    restored_state = restored.equilibrium_state("tp", 3000.0, 2.0e6)

    assert restored.MR == transport.MR
    assert restored.p_c == transport.p_c
    assert restored.minimum_cea_temperature == pytest.approx(200.0)
    assert restored_state.T == pytest.approx(expected_state.T)
    assert restored_state.p == pytest.approx(expected_state.p)
    assert restored_state.h == pytest.approx(expected_state.h)


def test_conditioned_gas_states_are_not_cached() -> None:
    fuel = psf.common.Fluid("fuel", ["CH4"], [1.0])
    oxidizer = psf.common.Fluid("oxidizer", ["O2"], [1.0])
    transport = psf.skycea.Aerothermodynamics.from_F_eps_Lstar(
        fu=fuel,
        ox=oxidizer,
        MR=2.8,
        p_c=1.0e7,
        F=1.0e5,
        eps=10.0,
        L_star=1.1,
    )
    base_key = (0.0, None, None, transport.transport)
    transport._cache[base_key] = GasState(p=2.0e6)
    calls = 0

    def equilibrium_state(problem, value, pressure):
        nonlocal calls
        calls += 1
        return GasState(T=float(value), p=float(pressure))

    transport.equilibrium_state = equilibrium_state

    transport.get_state(0.0, T=1000.0)
    transport.get_state(0.0, T=1000.0)

    assert calls == 2
    assert transport.cache_size == 1


def test_temperature_conditioned_gas_state_uses_c1_tangent_extrapolation() -> None:
    fuel = psf.common.Fluid("fuel", ["CH4"], [1.0])
    oxidizer = psf.common.Fluid("oxidizer", ["O2"], [1.0])
    transport = psf.skycea.Aerothermodynamics.from_F_eps_Lstar(
        fu=fuel,
        ox=oxidizer,
        MR=2.8,
        p_c=1.0e7,
        F=1.0e5,
        eps=10.0,
        L_star=1.1,
    )
    base_key = (0.0, None, None, transport.transport)
    transport._cache[base_key] = GasState(p=2.0e6)
    calls = []

    def equilibrium_state(problem, value, pressure):
        calls.append((problem, value, pressure))
        value = float(value)
        return GasState(
            T=value,
            p=float(pressure),
            rho=np.exp(0.001 * value),
            cp=np.exp(0.002 * value),
            gamma=np.exp(0.0001 * value),
            h=2.0 * value + 5.0,
            a=np.exp(0.003 * value),
            mu=np.exp(-0.001 * value),
            k=np.exp(0.01 * value),
            Pr=np.exp(0.0005 * value),
            mw=np.exp(-0.0002 * value),
            X={"test": 1.0},
        )

    transport.equilibrium_state = equilibrium_state

    state = transport.get_state(0.0, T=150.0)

    assert state.T == pytest.approx(150.0)
    assert state.h == pytest.approx(305.0)
    assert state.k == pytest.approx(np.exp(1.5), rel=1e-6)
    assert state.X == {"test": 1.0}
    assert calls == [
        ("tp", 200.0, 2.0e6),
        ("tp", 200.2, 2.0e6),
        ("tp", 200.4, 2.0e6),
    ]


def test_temperature_extrapolation_boundary_is_dynamic() -> None:
    fuel = psf.common.Fluid("fuel", ["CH4"], [1.0])
    oxidizer = psf.common.Fluid("oxidizer", ["O2"], [1.0])
    transport = psf.skycea.Aerothermodynamics.from_F_eps_Lstar(
        fu=fuel,
        ox=oxidizer,
        MR=2.8,
        p_c=1.0e7,
        F=1.0e5,
        eps=10.0,
        L_star=1.1,
        minimum_cea_temperature=400.0,
    )
    base_key = (0.0, None, None, transport.transport)
    transport._cache[base_key] = GasState(p=2.0e6)
    calls = []

    def equilibrium_state(problem, value, pressure):
        calls.append((problem, value, pressure))
        value = float(value)
        return GasState(
            T=value,
            p=float(pressure),
            rho=1.0,
            cp=2.0,
            gamma=1.2,
            h=3.0 * value,
            a=4.0,
            mu=5.0,
            k=np.exp(0.005 * value),
            Pr=6.0,
            mw=7.0,
        )

    transport.equilibrium_state = equilibrium_state

    state = transport.get_state(0.0, T=350.0)

    assert transport.minimum_cea_temperature == pytest.approx(400.0)
    assert state.T == pytest.approx(350.0)
    assert state.h == pytest.approx(1050.0)
    assert state.k == pytest.approx(np.exp(0.005 * 350.0), rel=1e-6)
    assert calls == [
        ("tp", 400.0, 2.0e6),
        ("tp", 400.4, 2.0e6),
        ("tp", 400.8, 2.0e6),
    ]

    transport.minimum_cea_temperature = 300.0
    assert transport.minimum_cea_temperature == pytest.approx(300.0)
    assert transport._temperature_tangent_cache == {}

    state = transport.get_state(0.0, T=250.0)
    assert state.T == pytest.approx(250.0)
    assert state.h == pytest.approx(750.0)
    assert calls[-3:] == [
        ("tp", 300.0, 2.0e6),
        ("tp", 300.3, 2.0e6),
        ("tp", 300.6, 2.0e6),
    ]


@pytest.mark.parametrize("value", [0.0, -1.0, np.inf, np.nan])
def test_temperature_extrapolation_boundary_must_be_positive(value) -> None:
    transport = psf.skycea.Aerothermodynamics.__new__(
        psf.skycea.Aerothermodynamics
    )
    transport._temperature_tangent_cache = {}

    with pytest.raises(ValueError, match="finite and positive"):
        transport.minimum_cea_temperature = value
