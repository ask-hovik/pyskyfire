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


def test_rl10_uses_stainless_steel_347() -> None:
    regen_sim = Path(__file__).parents[2] / "validation" / "RL10" / "regen_sim.py"

    assert "material = psf.common.solids.StainlessSteel347" in regen_sim.read_text()


def test_grcop42_am_follows_chen_equation_over_measured_range() -> None:
    """Chen et al. (2023) Eq. 3, T in K/1000, over 25-700 °C."""
    def chen(T):
        t = T / 1000.0
        return -39.86 * t**3 + 4.17 * t**2 + 0.97 * t + 329.14

    grcop = solids.GRCop42_AM
    for T in (298.15, 473.15, 773.15, 973.15):
        assert grcop.get_k(T) == pytest.approx(chen(T))
    assert solids.GRCop42 is grcop
    # About 85 % of pure copper at room temperature, never above it.
    assert grcop.get_k(298.15) == pytest.approx(328.7, abs=0.1)


def test_grcop84_am_fits_both_digitised_builds() -> None:
    grcop = solids.GRCop84_AM
    builds = solids.k_GRCop84_builds
    k = grcop.get_k(solids.T_GRCop84)
    assert k == pytest.approx(builds.mean(axis=0), abs=1.0)
    assert (builds[1] < k).all() and (k < builds[0]).all()
    # GRCop-42 conducts better throughout the measured range.
    assert (solids.GRCop42_AM.get_k(solids.T_GRCop84) > k).all()


@pytest.mark.parametrize("grcop", [solids.GRCop42_AM, solids.GRCop84_AM])
def test_grcop_am_is_continuous_and_silent_from_cryogenic_to_solidus(grcop) -> None:
    import warnings

    import numpy as np

    T = np.linspace(20.0, grcop.k._Tmax, 400)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        k = grcop.get_k(T)
    assert np.all(np.abs(np.diff(k)) < 0.5)            # no jump at 298 or 973 K
    assert np.all(np.diff(k[T > 298.15]) <= 1e-9)       # falls above room temperature
    assert np.all(k[T < 298.15] == grcop.get_k(298.15))  # held below the data
    assert grcop.get_alpha(293.15) == pytest.approx(
        1.51e-5 if grcop is solids.GRCop42_AM else 1.59e-5, rel=0.03
    )


def test_silica_deposit_matches_carwile_hoge_table() -> None:
    """Carwile and Hoge (1966) Table 1, in millical/(cm s K), converted to SI."""
    silica = solids.SilicaDeposit
    millical_per_cm_s_k = 0.4184

    assert silica.get_k(80.0) == pytest.approx(1.31 * millical_per_cm_s_k)
    assert silica.get_k(300.0) == pytest.approx(3.28 * millical_per_cm_s_k)
    assert silica.get_k(660.0) == pytest.approx(4.28 * millical_per_cm_s_k)


def test_silica_deposit_rises_gently_to_softening() -> None:
    """The non-radiative extension keeps rising, but stays near 2 W/m·K."""
    import warnings

    import numpy as np

    T = np.linspace(300.0, 1960.0, 200)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        k = solids.SilicaDeposit.get_k(T)
    assert np.all(np.diff(k) >= 0.0)
    assert k[-1] == pytest.approx(2.44, abs=0.01)
