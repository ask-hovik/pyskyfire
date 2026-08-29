"""Cross-constructor and nozzle-station contracts for the CEA wrapper."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.common import Fluid
from pyskyfire.skycea import Aerothermodynamics


@pytest.fixture(scope="module")
def equivalent_designs():
    fuel = Fluid("fuel", ["CH4"], [1.0])
    oxidizer = Fluid("oxidizer", ["O2"], [1.0])
    common = dict(
        fu=fuel,
        ox=oxidizer,
        MR=2.8,
        p_c=3.0e6,
        eps=8.0,
        L_star=1.0,
        T_fu_in=300.0,
        T_ox_in=300.0,
    )
    thrust_sized = Aerothermodynamics.from_F_eps_Lstar(F=50_000.0, **common)
    flow_sized = Aerothermodynamics.from_mdot_eps_Lstar(
        mdot=thrust_sized.mdot, **common
    )
    area_sized = Aerothermodynamics.from_At_eps_Lstar(
        A_t=thrust_sized.A_t, **common
    )
    return thrust_sized, flow_sized, area_sized


def test_equivalent_sizing_constructors_produce_the_same_design(equivalent_designs) -> None:
    reference, *others = equivalent_designs

    for design in others:
        assert design.F == pytest.approx(reference.F, rel=1e-10)
        assert design.mdot == pytest.approx(reference.mdot, rel=1e-10)
        assert design.A_t == pytest.approx(reference.A_t, rel=1e-10)
        assert design.c_star == pytest.approx(reference.c_star, rel=1e-10)
        assert design.Isp_vac == pytest.approx(reference.Isp_vac, rel=1e-10)


def test_mass_flow_constructor_closes_propellant_balance(equivalent_designs) -> None:
    reference = equivalent_designs[0]
    rebuilt = Aerothermodynamics.from_mass_flows_eps_Lstar(
        fu=reference.fu,
        ox=reference.ox,
        p_c=reference.p_c,
        mdot_fu=reference.mdot_fu,
        mdot_ox=reference.mdot_ox,
        eps=reference.eps,
        L_star=reference.L_star,
        T_fu_in=reference.T_fu_in,
        T_ox_in=reference.T_ox_in,
    )

    assert rebuilt.mdot == pytest.approx(reference.mdot, rel=1e-10)
    assert rebuilt.MR == pytest.approx(reference.MR, rel=1e-10)


def test_attached_contour_selects_subsonic_throat_and_supersonic_branches(
    equivalent_designs,
) -> None:
    aero = equivalent_designs[0]
    throat_area = aero.A_t
    contour = SimpleNamespace(
        A_t=throat_area,
        A=lambda x: throat_area * ({-1.0: 4.0, 0.0: 1.0, 1.0: 8.0}[float(x)]),
    )
    aero.attach_contour(contour)

    chamber = aero.get_state(-1.0)
    throat = aero.get_state(0.0)
    exit_state = aero.get_state(1.0)

    assert chamber.M < 1.0
    assert throat.M == pytest.approx(1.0, rel=2e-3)
    assert exit_state.M > 1.0
    assert chamber.p > throat.p > exit_state.p > 0.0
    assert all(
        value > 0.0
        for state in (chamber, throat, exit_state)
        for value in (state.T, state.rho, state.cp, state.mu, state.k, state.Pr)
    )


def test_total_enthalpy_temperature_is_not_below_static_temperature(
    equivalent_designs,
) -> None:
    aero = equivalent_designs[0]
    assert aero.get_T0(1.0) >= aero.get_T(1.0)
