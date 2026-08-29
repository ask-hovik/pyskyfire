"""Limiting cases and conservation behavior in the film-cooling models."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen.coupled_solver import FilmHeatFluxModel
from pyskyfire.regen.film_solver import (
    GasEmittanceCalculator,
    GaseousFilmResults,
    LiquidFilmResults,
    LiquidFilmSolver,
)


def test_gas_emittance_is_zero_without_radiating_species_and_bounded() -> None:
    transparent = SimpleNamespace(mole_fraction_H2O=0.0, mole_fraction_CO2=0.0)
    radiating = SimpleNamespace(mole_fraction_H2O=0.3, mole_fraction_CO2=0.1)

    assert GasEmittanceCalculator(transparent).gas_emittance(2_000.0, 0.2, 3e6) == 0.0
    epsilon = GasEmittanceCalculator(radiating).gas_emittance(2_000.0, 0.2, 3e6)
    assert 0.0 < epsilon <= 1.0


def test_liquid_film_marches_mass_per_perimeter_until_dryout() -> None:
    solver = object.__new__(LiquidFilmSolver)
    solver.film = SimpleNamespace(x=0.0)
    solver._Gamma_0 = 0.5
    solver._T_sat = 100.0
    solver._evaporation_rate = lambda x, Gamma, x_inj: (0.2, 10.0, 20.0, 0.0)

    result = solver.solve(np.arange(0.0, 4.0, 1.0))

    np.testing.assert_allclose(result.Gamma, [0.5, 0.3, 0.1])
    assert result.x_dryout == pytest.approx(2.5)
    assert result.Gamma_dryout == 0.0


def test_liquid_film_without_heat_load_does_not_dry_out() -> None:
    solver = object.__new__(LiquidFilmSolver)
    solver.film = SimpleNamespace(x=1.0)
    solver._Gamma_0 = 0.5
    solver._T_sat = 100.0
    solver._evaporation_rate = lambda x, Gamma, x_inj: (0.0, 0.0, 0.0, 0.0)

    result = solver.solve([0.0, 1.0, 2.0, 3.0])

    assert result.x == [1.0, 2.0, 3.0]
    np.testing.assert_allclose(result.Gamma, 0.5)
    assert result.x_dryout is None
    assert result.Gamma_dryout == pytest.approx(0.5)


def test_film_boundary_model_switches_from_liquid_to_vapour_to_bare_wall() -> None:
    liquid = LiquidFilmResults(
        x=[1.0, 2.0],
        T_film=[100.0, 110.0],
        x_dryout=2.5,
    )
    gaseous = GaseousFilmResults(
        x=[2.5, 3.5],
        T_aw=[200.0, 300.0],
        Q_rad=[1_000.0, 2_000.0],
    )
    model = FilmHeatFluxModel(liquid, gaseous, x_injection=1.0, h_liquid_wall=5_000.0)

    assert model.state(0.5) is None
    assert model.state(1.5) == {
        "regime": "liquid_film",
        "T_drive": pytest.approx(105.0),
        "q_rad": 0.0,
        "h_wall": 5_000.0,
    }
    assert model.state(3.0) == {
        "regime": "gaseous_film",
        "T_drive": pytest.approx(250.0),
        "q_rad": pytest.approx(1_500.0),
        "h_wall": None,
    }
    assert model.state(4.0) is None
