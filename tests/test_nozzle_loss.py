"""Nozzle and combustion efficiency corrections, and the RL10A-3-3A case.

The validation block at the bottom reproduces the delivered performance quoted
in Binder, M., Tomsik, T., and Veres, J. P., "RL10A-3-3A Rocket Engine Modeling
Project", NASA TM-107318, 1997 (NTRS 19970010379), from the ideal CEA solution
plus the report's own tabulated efficiencies.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.common.fluids import Fluid
from pyskyfire.skycea import aerothermodynamics
from pyskyfire.skycea.nozzle_loss import (
    G0,
    PSI_TO_PA,
    RL10A_3_3A_ETA_CF,
    RL10A_3_3A_LOSS,
    EfficiencyTable,
    NozzleLoss,
)


# Table 2.5.1, "Summary of Combustion-Chamber / Nozzle Characteristics", p. 8.
# The "Throat Diameter" and "c-star" rows of that table are mislabelled (they
# are a radius and ft/sec respectively), so only the unambiguous rows are used.
LBM = 0.45359237
LBF = 4.4482216
FT = 0.3048
PAPER_P_C = 482.0 * PSI_TO_PA          # injector-face static
PAPER_MDOT = 37.36 * LBM
PAPER_MDOT_FU = 5.973 * LBM            # table 2.4.1, cooling-jacket flow
PAPER_MDOT_OX = PAPER_MDOT - PAPER_MDOT_FU
PAPER_EPS = 61.0                       # with the silver throat insert
PAPER_ISP = 440.3                      # s, delivered
PAPER_THRUST = 16452.0 * LBF           # N, gross
PAPER_ETA_CSTAR = 0.9892
PAPER_CD = 0.975                       # sec. 4.2.3.3, value adopted for the model


def _linear_table(**kwargs):
    """Table whose values are 0.9 + 0.01*MR + 0.001*(p_c/1e5), linear in both."""
    mixture_ratio = np.array([1.0, 3.0, 7.0])
    chamber_pressure = np.array([1.0e5, 5.0e5, 2.0e6])
    values = (
        0.9
        + 0.01 * mixture_ratio[:, None]
        + 0.001 * (chamber_pressure[None, :] / 1.0e5)
    )
    return EfficiencyTable(mixture_ratio, chamber_pressure, values, **kwargs)


def _fake_aero(**overrides):
    """Minimal stand-in for an Aerothermodynamics design point."""
    design = dict(
        p_c=3.0e6,
        MR=5.0,
        eps=60.0,
        mdot=17.0,
        p_amb=0.0,
        c_star=2400.0,
        CF_vac=1.9,
    )
    design.update(overrides)
    return SimpleNamespace(**design)


# ======================================================================
# EfficiencyTable
# ======================================================================

def test_table_reproduces_every_grid_node_exactly() -> None:
    table = _linear_table()
    for i, mixture_ratio in enumerate(table.mixture_ratio):
        for j, chamber_pressure in enumerate(table.chamber_pressure):
            assert table(chamber_pressure, mixture_ratio) == pytest.approx(
                table.values[i, j], rel=1e-12
            )


def test_table_interpolates_linearly_in_mixture_ratio() -> None:
    table = _linear_table(log_pressure=False)
    # Midway between the MR = 1 and MR = 3 rows at a grid pressure.
    expected = 0.5 * (table.values[0, 1] + table.values[1, 1])

    assert table(5.0e5, 2.0) == pytest.approx(expected, rel=1e-12)


def test_table_interpolates_linearly_in_pressure_when_log_disabled() -> None:
    table = _linear_table(log_pressure=False)
    # The table is exactly linear in p_c, so linear mode is exact off-node.
    assert table(1.25e6, 3.0) == pytest.approx(0.9 + 0.03 + 0.0125, rel=1e-12)


def test_table_log_pressure_mode_is_exact_for_a_log_linear_table() -> None:
    chamber_pressure = np.array([1.0e4, 1.0e5, 1.0e6])
    mixture_ratio = np.array([2.0, 4.0])
    values = 0.5 + 0.05 * np.log(chamber_pressure)[None, :] * np.ones((2, 1))
    table = EfficiencyTable(mixture_ratio, chamber_pressure, values)

    midpoint = np.sqrt(1.0e4 * 1.0e5)  # geometric mean -> midpoint in log
    assert table(midpoint, 3.0) == pytest.approx(
        0.5 + 0.05 * np.log(midpoint), rel=1e-12
    )


def test_table_log_and_linear_modes_disagree_off_node() -> None:
    """Guard against the two modes being silently identical."""
    log_table = _linear_table(log_pressure=True)
    linear_table = _linear_table(log_pressure=False)

    assert log_table(1.25e6, 3.0) != pytest.approx(
        linear_table(1.25e6, 3.0), rel=1e-6
    )


@pytest.mark.parametrize(
    "chamber_pressure, mixture_ratio, node",
    [
        (1.0e3, 1.0, (0, 0)),    # below both axes
        (1.0e8, 1.0, (0, -1)),   # above the pressure axis
        (1.0e5, 0.1, (0, 0)),    # below the mixture-ratio axis
        (1.0e5, 99.0, (-1, 0)),  # above the mixture-ratio axis
        (1.0e8, 99.0, (-1, -1)), # beyond both
    ],
)
def test_table_clamps_outside_the_grid(chamber_pressure, mixture_ratio, node) -> None:
    table = _linear_table()

    assert table(chamber_pressure, mixture_ratio) == pytest.approx(
        table.values[node], rel=1e-12
    )


def test_table_broadcasts_array_arguments() -> None:
    table = _linear_table()
    pressures = np.array([1.0e5, 5.0e5, 2.0e6])
    result = table(pressures, 3.0)

    assert result.shape == pressures.shape
    np.testing.assert_allclose(result, table.values[1], rtol=1e-12)

    grid = table(pressures[None, :], table.mixture_ratio[:, None])
    assert grid.shape == (3, 3)
    np.testing.assert_allclose(grid, table.values, rtol=1e-12)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        (dict(mixture_ratio=[1.0], chamber_pressure=[1.0e5, 2.0e5],
              values=[[0.9, 0.9]]), "at least two"),
        (dict(mixture_ratio=[3.0, 1.0], chamber_pressure=[1.0e5, 2.0e5],
              values=[[0.9, 0.9], [0.9, 0.9]]), "strictly increasing"),
        (dict(mixture_ratio=[1.0, 3.0], chamber_pressure=[2.0e5, 1.0e5],
              values=[[0.9, 0.9], [0.9, 0.9]]), "strictly increasing"),
        (dict(mixture_ratio=[1.0, 3.0], chamber_pressure=[1.0e5, 2.0e5],
              values=[[0.9, 0.9, 0.9], [0.9, 0.9, 0.9]]), "shape"),
        (dict(mixture_ratio=[1.0, 3.0], chamber_pressure=[1.0e5, 2.0e5],
              values=[[0.9, np.nan], [0.9, 0.9]]), "finite"),
        (dict(mixture_ratio=[1.0, 3.0], chamber_pressure=[1.0e5, 2.0e5],
              values=[[0.9, -0.1], [0.9, 0.9]]), "positive"),
        (dict(mixture_ratio=[1.0, 3.0], chamber_pressure=[0.0, 2.0e5],
              values=[[0.9, 0.9], [0.9, 0.9]]), "positive chamber_pressure"),
    ],
)
def test_table_rejects_malformed_input(kwargs, message) -> None:
    with pytest.raises(ValueError, match=message):
        EfficiencyTable(**kwargs)


def test_table_rejects_non_positive_pressure_query_in_log_mode() -> None:
    table = _linear_table()

    with pytest.raises(ValueError, match="must be positive"):
        table(0.0, 3.0)


# ======================================================================
# NozzleLoss
# ======================================================================

def test_identity_loss_leaves_performance_untouched() -> None:
    loss = NozzleLoss(eta_cf=1.0, eta_cstar=1.0)

    assert loss.specific_impulse(450.0, 3.0e6, 5.0) == pytest.approx(450.0)
    assert loss.overall_efficiency(3.0e6, 5.0) == pytest.approx(1.0)


def test_constant_efficiencies_multiply() -> None:
    loss = NozzleLoss(eta_cf=0.95, eta_cstar=0.98)

    assert loss.overall_efficiency(3.0e6, 5.0) == pytest.approx(0.95 * 0.98)
    assert loss.specific_impulse(450.0, 3.0e6, 5.0) == pytest.approx(450.0 * 0.95 * 0.98)
    assert loss.characteristic_velocity(2400.0, 3.0e6, 5.0) == pytest.approx(2400.0 * 0.98)
    assert loss.thrust_coefficient(1.9, 3.0e6, 5.0) == pytest.approx(1.9 * 0.95)


def test_callable_efficiency_receives_pressure_and_mixture_ratio() -> None:
    seen = []

    def eta(p_c, MR):
        seen.append((p_c, MR))
        return 0.9

    loss = NozzleLoss(eta_cf=eta, eta_cstar=1.0)

    assert loss.thrust_coefficient_efficiency(2.5e6, 4.2) == pytest.approx(0.9)
    assert seen == [(2.5e6, 4.2)]


def test_table_backed_efficiency_matches_a_direct_table_lookup() -> None:
    table = _linear_table()
    loss = NozzleLoss(eta_cf=table, eta_cstar=1.0)

    assert loss.thrust_coefficient_efficiency(7.5e5, 2.5) == pytest.approx(
        table(7.5e5, 2.5)
    )


def test_mass_flow_applies_both_discharge_coefficient_and_cstar_efficiency() -> None:
    loss = NozzleLoss(eta_cf=1.0, eta_cstar=0.99, Cd=0.97)
    p_c, A_t, c_star_ideal = 3.0e6, 0.012, 2400.0

    expected = 0.97 * p_c * A_t / (0.99 * c_star_ideal)

    assert loss.mass_flow(p_c, A_t, c_star_ideal, 5.0) == pytest.approx(expected)


def test_discharge_coefficient_does_not_touch_specific_impulse() -> None:
    """Cd rescales flow through a throat, not the impulse per unit of flow."""
    lossy = NozzleLoss(eta_cf=0.95, eta_cstar=0.98, Cd=0.90)
    clean = NozzleLoss(eta_cf=0.95, eta_cstar=0.98, Cd=1.00)

    assert lossy.specific_impulse(450.0, 3.0e6, 5.0) == pytest.approx(
        clean.specific_impulse(450.0, 3.0e6, 5.0)
    )


@pytest.mark.parametrize("bad", [0.0, -0.5, float("nan"), float("inf")])
def test_rejects_non_physical_constant_efficiency(bad) -> None:
    with pytest.raises(ValueError, match="finite positive"):
        NozzleLoss(eta_cf=bad)


@pytest.mark.parametrize("bad", [0.0, -0.5, float("nan")])
def test_rejects_non_physical_discharge_coefficient(bad) -> None:
    with pytest.raises(ValueError, match="Cd must be finite and positive"):
        NozzleLoss(eta_cf=0.95, Cd=bad)


def test_rejects_callable_returning_a_non_physical_value() -> None:
    loss = NozzleLoss(eta_cf=lambda p_c, MR: -1.0)

    with pytest.raises(ValueError, match="finite positive"):
        loss.thrust_coefficient_efficiency(3.0e6, 5.0)


# ----------------------------------------------------------------------
# apply()
# ----------------------------------------------------------------------

def test_apply_scales_cstar_and_vacuum_thrust_coefficient() -> None:
    aero = _fake_aero()
    delivered = NozzleLoss(eta_cf=0.95, eta_cstar=0.98).apply(aero)

    assert delivered.c_star == pytest.approx(aero.c_star * 0.98)
    assert delivered.CF_vac == pytest.approx(aero.CF_vac * 0.95)
    assert delivered.Isp_vac == pytest.approx(
        aero.CF_vac * 0.95 * aero.c_star * 0.98 / G0
    )


def test_apply_leaves_the_ambient_pressure_term_unscaled() -> None:
    """The p_amb/p_c*eps term is geometry and pressure -- it carries no loss."""
    aero = _fake_aero(p_amb=1.0e5)
    eta_cf, eta_cstar = 0.95, 0.98
    delivered = NozzleLoss(eta_cf=eta_cf, eta_cstar=eta_cstar).apply(aero)

    expected = aero.CF_vac * eta_cf - aero.p_amb / aero.p_c * aero.eps
    assert delivered.CF_amb == pytest.approx(expected)

    # Scaling the whole ambient coefficient would be wrong; confirm the two
    # treatments actually differ, so this test cannot pass by coincidence.
    naive = (aero.CF_vac - aero.p_amb / aero.p_c * aero.eps) * eta_cf
    assert delivered.CF_amb != pytest.approx(naive, rel=1e-6)


def test_apply_passes_mass_flow_through_unchanged() -> None:
    aero = _fake_aero()
    delivered = NozzleLoss(eta_cf=0.95, eta_cstar=0.98).apply(aero)

    assert delivered.mdot == pytest.approx(aero.mdot)
    assert delivered.F_vac == pytest.approx(aero.mdot * delivered.Isp_vac * G0)


def test_apply_ambient_override_beats_the_stored_value() -> None:
    aero = _fake_aero(p_amb=0.0)
    loss = NozzleLoss(eta_cf=1.0, eta_cstar=1.0)

    delivered = loss.apply(aero, p_amb=1.0e5)
    expected = aero.CF_vac - 1.0e5 / aero.p_c * aero.eps

    assert delivered.CF_amb == pytest.approx(expected)
    assert loss.apply(aero).CF_amb == pytest.approx(aero.CF_vac)


def test_module_gravity_constant_tracks_aerothermodynamics() -> None:
    """G0 is duplicated to keep this module CEA-free; catch any drift."""
    assert G0 == aerothermodynamics.G0


# ======================================================================
# RL10A-3-3A validation against the report
# ======================================================================

@pytest.mark.parametrize(
    "mixture_ratio, chamber_pressure_psia, published",
    [
        (1.0, 5.0, 0.9381),
        (1.0, 500.0, 0.9593),
        (4.0, 150.0, 0.9419),
        (6.0, 500.0, 0.9469),
        (8.0, 50.0, 0.8950),
        (10.0, 5.0, 0.8507),
        (20.0, 150.0, 0.9374),
        (50.0, 500.0, 0.9596),
    ],
)
def test_rl10_table_matches_published_table_e3(
    mixture_ratio, chamber_pressure_psia, published
) -> None:
    assert RL10A_3_3A_ETA_CF(
        chamber_pressure_psia * PSI_TO_PA, mixture_ratio
    ) == pytest.approx(published, rel=1e-12)


def test_rl10_loss_carries_the_reported_constants() -> None:
    assert RL10A_3_3A_LOSS.c_star_efficiency(PAPER_P_C, 5.26) == pytest.approx(
        PAPER_ETA_CSTAR
    )
    assert RL10A_3_3A_LOSS.Cd == pytest.approx(PAPER_CD)


def test_rl10_design_point_efficiency_sits_between_the_bracketing_rows() -> None:
    """O/F = 5.255 lies between the O/F = 4 and O/F = 6 rows at Pc = 500."""
    eta = RL10A_3_3A_LOSS.thrust_coefficient_efficiency(PAPER_P_C, 5.2548)

    assert 0.9469 < eta < 0.9517
    assert eta == pytest.approx(0.9482, abs=5e-4)


@pytest.fixture(scope="module")
def rl10_ideal():
    """Ideal ODE design point for the RL10A-3-3A at the reported conditions."""
    return aerothermodynamics.Aerothermodynamics.from_mass_flows_eps_Lstar(
        fu=Fluid(type="fuel", propellants=["H2"], fractions=[1.0]),
        ox=Fluid(type="oxidizer", propellants=["O2(L)"], fractions=[1.0]),
        T_fu_in=200.0,
        T_ox_in=100.0,
        mdot_fu=PAPER_MDOT_FU,
        mdot_ox=PAPER_MDOT_OX,
        p_c=PAPER_P_C,
        eps=PAPER_EPS,
        L_star=0.95,
        p_amb=0.0,
    )


def test_rl10_ideal_solution_matches_the_reported_operating_point(rl10_ideal) -> None:
    """Sanity-check the inputs before validating what is built on them."""
    assert rl10_ideal.mdot == pytest.approx(PAPER_MDOT, rel=1e-6)
    assert rl10_ideal.MR == pytest.approx(5.26, rel=2e-3)
    # The ideal solution must sit above the delivered value it will be reduced to.
    assert rl10_ideal.Isp_vac > PAPER_ISP


def test_rl10_delivered_isp_reproduces_the_reported_value(rl10_ideal) -> None:
    delivered = RL10A_3_3A_LOSS.apply(rl10_ideal)

    assert delivered.Isp_vac == pytest.approx(PAPER_ISP, rel=5e-3)


def test_rl10_delivered_thrust_reproduces_the_reported_value(rl10_ideal) -> None:
    delivered = RL10A_3_3A_LOSS.apply(rl10_ideal)

    assert delivered.F_vac == pytest.approx(PAPER_THRUST, rel=5e-3)


def test_rl10_loss_budget_is_dominated_by_the_nozzle(rl10_ideal) -> None:
    delivered = RL10A_3_3A_LOSS.apply(rl10_ideal)

    nozzle_loss = 1.0 - delivered.eta_cf
    combustion_loss = 1.0 - delivered.eta_cstar

    assert nozzle_loss == pytest.approx(0.052, abs=2e-3)
    assert combustion_loss == pytest.approx(0.0108, abs=2e-4)
    assert nozzle_loss > 4.0 * combustion_loss

    # The two together account for the whole ideal-to-delivered gap.
    assert delivered.eta_Isp == pytest.approx(
        PAPER_ISP / rl10_ideal.Isp_vac, rel=5e-3
    )
