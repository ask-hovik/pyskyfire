"""Pure-fluid phase and enthalpy-state contracts."""

import numpy as np
import pytest

from pyskyfire.common import Fluid
from pyskyfire.regen.coolant_state import (
    PHASE_LIQUID,
    PHASE_SUPERCRITICAL,
    PHASE_TWO_PHASE,
    PHASE_VAPOUR,
    CoolantState,
    probe_coolant,
)
from pyskyfire.skycea import CoolantTransport


def _water_transport():
    return CoolantTransport(Fluid("coolant", ["Water"], [1.0]))


def test_pure_water_supports_enthalpy_march_and_saturation() -> None:
    state = CoolantState(_water_transport(), 1.0e5)
    saturation = state.saturation(1.0e5)

    assert state.capability.enthalpy_march
    assert state.capability.saturation
    assert saturation is not None
    assert saturation.h_fg > 0.0
    assert saturation.rho_f > saturation.rho_g > 0.0


def test_enthalpy_flash_classifies_liquid_two_phase_and_vapour() -> None:
    state = CoolantState(_water_transport(), 1.0e5)
    saturation = state.saturation(1.0e5)
    assert saturation is not None

    liquid = state.from_hp(saturation.h_f - 10_000.0, 1.0e5)
    two_phase = state.from_hp(0.75 * saturation.h_f + 0.25 * saturation.h_g, 1.0e5)
    vapour = state.from_hp(saturation.h_g + 10_000.0, 1.0e5)

    assert liquid.phase == PHASE_LIQUID
    assert two_phase.phase == PHASE_TWO_PHASE
    assert two_phase.quality == pytest.approx(0.25)
    assert saturation.rho_g < two_phase.rho < saturation.rho_f
    assert vapour.phase == PHASE_VAPOUR


def test_supercritical_state_has_no_quality_or_saturation_dome() -> None:
    pressure = 25.0e6
    state = CoolantState(_water_transport(), pressure)
    enthalpy = state.enthalpy(700.0, pressure)
    bulk = state.from_hp(enthalpy, pressure)

    assert state.saturation(pressure) is None
    assert bulk.phase == PHASE_SUPERCRITICAL
    assert np.isnan(bulk.quality)
    assert bulk.T == pytest.approx(700.0, rel=1e-8)


def test_mixtures_explicitly_fall_back_from_enthalpy_march() -> None:
    capability = probe_coolant(
        "HEOS::Water[0.500]&Ethanol[0.500]",
        1.0e5,
    )

    assert capability.is_mixture
    assert not capability.enthalpy_march
    assert not capability.saturation
    assert "mixture" in capability.reason


def test_hot_film_properties_clamp_to_saturated_liquid() -> None:
    state = CoolantState(_water_transport(), 1.0e5)
    saturation = state.saturation(1.0e5)
    assert saturation is not None

    properties = state.film_properties(
        saturation.T_sat + 50.0,
        1.0e5,
        getters=(
            lambda T, p: pytest.fail("conductivity getter crossed saturation"),
            lambda T, p: pytest.fail("cp getter crossed saturation"),
            lambda T, p: pytest.fail("viscosity getter crossed saturation"),
        ),
    )

    assert properties == pytest.approx(
        (saturation.k_f, saturation.cp_f, saturation.mu_f)
    )
