"""Acceleration term of the coolant momentum equation.

The static-pressure march has to span two limits that a friction-only
stagnation march followed by ``p = p_0 - rho u^2 / 2`` does not: reversible
area-driven acceleration, and irreversible heating-driven acceleration in a
duct of fixed area. Both are checked against their closed forms here.
"""

import numpy as np
import pytest

from pyskyfire.regen.coupled_solver import _acceleration_pressure_drop


def test_constant_density_reduces_to_bernoulli() -> None:
    rho = 71.0
    u_up, u_down = 12.0, 260.0

    drop = _acceleration_pressure_drop(rho, u_up, rho, u_down)

    assert drop == pytest.approx(0.5 * rho * (u_down**2 - u_up**2))


def test_constant_area_reduces_to_mass_flux_form() -> None:
    # Constant area means constant mass flux, so the segment is characterised
    # by G alone and the drop must be G^2 (1/rho_down - 1/rho_up).
    flux = 3120.0
    rho_up, rho_down = 23.7, 9.06
    u_up, u_down = flux / rho_up, flux / rho_down

    drop = _acceleration_pressure_drop(rho_up, u_up, rho_down, u_down)

    assert drop == pytest.approx(flux**2 * (1.0 / rho_down - 1.0 / rho_up))


def test_heated_constant_area_drop_is_twice_the_dynamic_head_change() -> None:
    # This is the discrepancy the fix addresses: reconstructing static pressure
    # from a friction-only stagnation march counts exactly half of this.
    flux = 3120.0
    rho_up, rho_down = 23.7, 9.06
    u_up, u_down = flux / rho_up, flux / rho_down

    drop = _acceleration_pressure_drop(rho_up, u_up, rho_down, u_down)
    dynamic_head_change = 0.5 * rho_down * u_down**2 - 0.5 * rho_up * u_up**2

    assert drop == pytest.approx(2.0 * dynamic_head_change)


def test_deceleration_recovers_pressure() -> None:
    rho = 9.0
    drop = _acceleration_pressure_drop(rho, 340.0, rho, 150.0)

    assert drop < 0.0


def test_no_velocity_change_costs_nothing() -> None:
    assert _acceleration_pressure_drop(9.0, 150.0, 9.0, 150.0) == pytest.approx(0.0)


def test_subdividing_a_segment_is_consistent_at_constant_area() -> None:
    # At constant area the term telescopes exactly, so a refined coolant grid
    # cannot change the integrated acceleration loss.
    flux = 3120.0
    densities = np.array([23.7, 18.0, 13.5, 11.0, 9.06])
    velocities = flux / densities

    subdivided = sum(
        _acceleration_pressure_drop(
            densities[i], velocities[i], densities[i + 1], velocities[i + 1]
        )
        for i in range(len(densities) - 1)
    )
    whole = _acceleration_pressure_drop(
        densities[0], velocities[0], densities[-1], velocities[-1]
    )

    assert subdivided == pytest.approx(whole)
