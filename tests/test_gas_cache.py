"""Tests for bounded combustion-gas cache behavior."""

import pytest

from pyskyfire.regen.coupled_solver import _clear_gas_cache_after


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
