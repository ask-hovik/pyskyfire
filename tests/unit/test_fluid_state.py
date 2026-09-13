"""Pure-fluid phase and enthalpy-state contracts."""

import numpy as np
import pytest

from pyskyfire.regen.coolant_properties import film_properties
from pyskyfire.thermo import (
    PHASE_LIQUID,
    PHASE_SUPERCRITICAL,
    PHASE_TWO_PHASE,
    PHASE_VAPOUR,
    Fluid,
    FluidStateResolver,
    FluidTransport,
    probe_fluid,
)


def _water_transport():
    return FluidTransport(Fluid("coolant", ["Water"], [1.0]))


def test_pure_water_supports_enthalpy_march_and_saturation() -> None:
    resolver = FluidStateResolver(_water_transport(), 1.0e5)
    saturation = resolver.saturation(1.0e5)

    assert resolver.capability.enthalpy_march
    assert resolver.capability.saturation
    assert saturation is not None
    assert saturation.h_fg > 0.0
    assert saturation.rho_f > saturation.rho_g > 0.0


def test_enthalpy_flash_classifies_liquid_two_phase_and_vapour() -> None:
    resolver = FluidStateResolver(_water_transport(), 1.0e5)
    saturation = resolver.saturation(1.0e5)
    assert saturation is not None

    liquid = resolver.from_hp(saturation.h_f - 10_000.0, 1.0e5)
    two_phase = resolver.from_hp(
        0.75 * saturation.h_f + 0.25 * saturation.h_g, 1.0e5
    )
    vapour = resolver.from_hp(saturation.h_g + 10_000.0, 1.0e5)

    assert liquid.phase == PHASE_LIQUID
    assert two_phase.phase == PHASE_TWO_PHASE
    assert two_phase.quality == pytest.approx(0.25)
    assert saturation.rho_g < two_phase.rho < saturation.rho_f
    assert vapour.phase == PHASE_VAPOUR


def test_supercritical_state_has_no_quality_or_saturation_dome() -> None:
    pressure = 25.0e6
    resolver = FluidStateResolver(_water_transport(), pressure)
    enthalpy = resolver.enthalpy(700.0, pressure)
    fluid_state = resolver.from_hp(enthalpy, pressure)

    assert resolver.saturation(pressure) is None
    assert fluid_state.phase == PHASE_SUPERCRITICAL
    assert np.isnan(fluid_state.quality)
    assert fluid_state.T == pytest.approx(700.0, rel=1e-8)


def test_mixtures_explicitly_fall_back_from_enthalpy_march() -> None:
    capability = probe_fluid(
        "HEOS::Water[0.500]&Ethanol[0.500]",
        1.0e5,
    )

    assert capability.is_mixture
    assert not capability.enthalpy_march
    assert not capability.saturation
    assert "mixture" in capability.reason


def test_hot_film_properties_clamp_to_saturated_liquid() -> None:
    resolver = FluidStateResolver(_water_transport(), 1.0e5)
    saturation = resolver.saturation(1.0e5)
    assert saturation is not None

    properties = film_properties(
        resolver,
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
