"""Material-property definitions used by the thermal solver."""

from pathlib import Path

import pytest

from pyskyfire.common import solids


def test_pure_silver_uses_temperature_dependent_conductivity() -> None:
    silver = solids.PureSilver

    assert silver.name == "Pure silver"
    assert silver.get_k(300.0) == pytest.approx(429.0)
    assert silver.get_k(1000.0) == pytest.approx(379.0)
    assert silver.get_k(300.0) > silver.get_k(1000.0)
    assert silver.get_rho(293.0) == pytest.approx(10492.0)


def test_stainless_steel_347_matches_nerva_conductivity_table() -> None:
    """NERVA DRM M-2 page 4 values are converted from US units to SI."""
    steel = solids.StainlessSteel347
    btu_per_hr_ft_f_to_w_per_m_k = 1.730734666371391

    assert steel.name == "AISI 347"
    assert steel.get_k(((-420.0 - 32.0) * 5.0 / 9.0) + 273.15) == pytest.approx(
        1.39 * btu_per_hr_ft_f_to_w_per_m_k
    )
    assert steel.get_k(((540.0 - 32.0) * 5.0 / 9.0) + 273.15) == pytest.approx(
        10.78 * btu_per_hr_ft_f_to_w_per_m_k
    )
    assert steel.get_k(((1840.0 - 32.0) * 5.0 / 9.0) + 273.15) == pytest.approx(
        16.11 * btu_per_hr_ft_f_to_w_per_m_k
    )


def test_stainless_steel_347_matches_nerva_density_table() -> None:
    """NERVA DRM M-2 page 2 density values are converted to kg/m^3."""
    steel = solids.StainlessSteel347
    lb_per_in3_to_kg_per_m3 = 27679.904710191

    assert steel.get_rho(((80.6 - 32.0) * 5.0 / 9.0) + 273.15) == pytest.approx(
        0.285 * lb_per_in3_to_kg_per_m3
    )
    assert steel.get_rho(((801.0 - 32.0) * 5.0 / 9.0) + 273.15) == pytest.approx(
        0.280 * lb_per_in3_to_kg_per_m3
    )


def test_rl10_2_uses_stainless_steel_347() -> None:
    regen_sim = Path(__file__).parents[1] / "validation" / "RL10_2" / "regen_sim.py"

    assert "material = psf.common.solids.StainlessSteel347" in regen_sim.read_text()
