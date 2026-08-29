"""Analytic limits for regenerative-cooling correlations."""

import numpy as np
import pytest

from pyskyfire.regen import physics


def test_velocity_and_reynolds_number_follow_their_definitions() -> None:
    velocity = physics.u_coolant(rho=20.0, mdot_c_single_channel=0.4, A_cool=0.002)

    assert velocity == pytest.approx(10.0)
    assert physics.reynolds(20.0, velocity, 0.01, 1.0e-3) == pytest.approx(2_000.0)


def test_colburn_and_bartz_obey_expected_scaling() -> None:
    colburn = physics.h_coolant_colburn(0.2, 0.01, 4_000.0, 1e-3, 0.1, 1e-4)
    doubled_flow = physics.h_coolant_colburn(0.2, 0.01, 4_000.0, 1e-3, 0.2, 1e-4)
    assert doubled_flow / colburn == pytest.approx(2.0**0.8)

    bartz = physics.h_gas_bartz(0.1, 1e-4, 2_000.0, 0.7, 3e6, 1_500.0, 0.01, 0.02, 1.0)
    throat = physics.h_gas_bartz(0.1, 1e-4, 2_000.0, 0.7, 3e6, 1_500.0, 0.01, 0.01, 1.0)
    assert throat / bartz == pytest.approx(2.0**0.9)


def test_darcy_factor_has_laminar_limit_and_continuous_transition() -> None:
    assert physics.f_darcy(1_000.0, 0.01, 0.0, None) == pytest.approx(64.0 / 1_000.0)

    below = physics.f_darcy(physics.ReDh_laminar - 1e-6, 0.01, 0.0, None)
    at = physics.f_darcy(physics.ReDh_laminar, 0.01, 0.0, None)
    assert at == pytest.approx(below, rel=1e-6)
    assert physics.f_darcy(100_000.0, 0.01, 0.0, lambda x: 1e-5) > 0.0


def test_recovery_temperature_has_static_and_stagnation_limits() -> None:
    assert physics.T_aw(1.4, 0.0, 300.0, 0.7) == pytest.approx(300.0)

    recovered = physics.T_aw(1.4, 2.0, 300.0, 1.0)
    ideal_stagnation = 300.0 * (1.0 + 0.5 * (1.4 - 1.0) * 2.0**2)
    assert recovered == pytest.approx(ideal_stagnation)
