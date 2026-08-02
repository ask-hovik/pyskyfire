"""Persistence tests for objects backed by native extension handles."""

import cloudpickle
import pytest

import pyskyfire as psf


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
    restored = cloudpickle.loads(cloudpickle.dumps(transport))
    restored_state = restored.equilibrium_state("tp", 3000.0, 2.0e6)

    assert restored.MR == transport.MR
    assert restored.p_c == transport.p_c
    assert restored_state.T == pytest.approx(expected_state.T)
    assert restored_state.p == pytest.approx(expected_state.p)
    assert restored_state.h == pytest.approx(expected_state.h)
