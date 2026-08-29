"""Tests for bounded combustion-gas cache behavior."""

from types import SimpleNamespace

import pytest

from pyskyfire.regen.coupled_solver import _clear_gas_cache_after
from pyskyfire.skycea.aerothermodynamics import Aerothermodynamics


class _CombustionTransport:
    def __init__(self):
        self.cache = {}

    def clear_cache(self):
        self.cache.clear()


class _ThrustChamber:
    def __init__(self):
        self.combustion_transport = _CombustionTransport()


@_clear_gas_cache_after
def _simulation(thrust_chamber, *, fail=False):
    thrust_chamber.combustion_transport.cache["gas-state"] = object()
    if fail:
        raise RuntimeError("simulation failed")
    return "result"


def test_gas_cache_is_cleared_after_simulation() -> None:
    thrust_chamber = _ThrustChamber()

    assert _simulation(thrust_chamber) == "result"
    assert thrust_chamber.combustion_transport.cache == {}


def test_gas_cache_is_cleared_when_simulation_fails() -> None:
    thrust_chamber = _ThrustChamber()

    with pytest.raises(RuntimeError, match="simulation failed"):
        _simulation(thrust_chamber, fail=True)

    assert thrust_chamber.combustion_transport.cache == {}


class _RocketSolver:
    """Minimal stand-in for the native CEA rocket handle."""

    def __init__(self):
        self.calls = 0

    def solve(self, solution, weights, pressure, **station):
        self.calls += 1


def _station_transport():
    """Build a transport object whose only live part is the rocket solve."""
    transport = Aerothermodynamics.__new__(Aerothermodynamics)
    transport.transport = "equilibrium"
    transport._cache_enabled = True
    transport._cache = {}
    transport._temperature_tangent_cache = {}
    transport.p_c = 3.0e6
    transport.reactant_weights = [1.0]
    transport._reactant_enthalpy = 0.0
    transport.contour = SimpleNamespace(A=lambda x: 2.0, A_t=1.0)
    transport._rocket_solver = _RocketSolver()
    transport._rocket_solution = SimpleNamespace(
        converged=True, T=2500.0, P=20.0, Mach=2.0, density=0.5,
        cp_eq=2.0, gamma_s=1.2, enthalpy=-1.0, sonic_velocity=1000.0,
        viscosity=1.0, conductivity_eq=1.0, Pr_eq=0.6, MW=20.0,
        mole_fractions={"H2O": 0.5},
    )
    return transport


def test_station_states_are_solved_once_per_axial_coordinate() -> None:
    transport = _station_transport()

    first = transport.get_state(0.4)
    second = transport.get_state(0.4)

    assert second is first
    assert transport._rocket_solver.calls == 1

    transport.get_state(0.5)
    assert transport._rocket_solver.calls == 2


def test_clearing_the_cache_forces_a_new_station_solve() -> None:
    transport = _station_transport()

    transport.get_state(0.4)
    transport.clear_cache()
    transport.get_state(0.4)

    assert transport._rocket_solver.calls == 2


def test_reset_solver_state_recreates_native_handles_and_clears_cache() -> None:
    transport = Aerothermodynamics.__new__(Aerothermodynamics)
    transport._cache = {"stale": object()}
    calls = []
    transport._make_solvers = lambda: calls.append("reset")

    transport.reset_solver_state()

    assert calls == ["reset"]
    assert transport._cache == {}
