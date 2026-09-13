"""Manufactured full-solver cases with closed-form heat and pressure marches."""

from types import SimpleNamespace

import fluids
import numpy as np
import pytest

from pyskyfire.regen import coupled_solver
from pyskyfire.thermo.state import FluidState, PHASE_UNKNOWN


T_GAS = 1_000.0
H_HOT = 100.0
H_COLD = 300.0
CP = 1_000.0
MDOT = 1.0
FRICTION_GRADIENT = -1_000.0


def _transition_weight(reynolds, critical):
    end = critical * 3_500.0 / 2_300.0
    fraction = np.clip((reynolds - critical) / (end - critical), 0.0, 1.0)
    return fraction**2 * (3.0 - 2.0 * fraction)


def _straight_factor(reynolds, diameter, roughness):
    laminar = fluids.friction.friction_laminar(reynolds)
    weight = _transition_weight(reynolds, 2_040.0)
    if weight == 0.0:
        return laminar
    turbulent = fluids.friction.Clamond(reynolds, roughness / diameter)
    return (1.0 - weight) * laminar + weight * turbulent


def _curved_factor(reynolds, diameter, radius, roughness):
    coil_diameter = 2.0 * abs(radius)
    critical = fluids.friction.helical_transition_Re_Schmidt(
        diameter, coil_diameter,
    )
    weight = _transition_weight(reynolds, critical)
    laminar = fluids.friction.helical_laminar_fd_Schmidt(
        reynolds, diameter, coil_diameter,
    )
    if weight == 0.0:
        return laminar
    turbulent = fluids.friction.helical_turbulent_fd_Schmidt(
        reynolds, diameter, coil_diameter, roughness,
    )
    return (1.0 - weight) * laminar + weight * turbulent


class _FluidStateResolver:
    """Constant-property fluid-state resolver for an independent oracle."""

    def __init__(self, transport, p_ref):
        self.capability = SimpleNamespace(
            enthalpy_march=False,
            reason="manufactured constant-property coolant.",
        )

    def from_tp(self, temperature, pressure):
        return FluidState(
            T=float(temperature),
            p=float(pressure),
            h=np.nan,
            quality=np.nan,
            rho=1_000.0,
            phase=PHASE_UNKNOWN,
            T_sat=np.nan,
        )

    def saturation(self, pressure):
        return None


class _Circuit:
    name = "manufactured"
    n_channels = 1
    walls = [object()]
    is_helical = False

    def __init__(self, direction):
        self.direction = direction
        self.x_domain = np.array([0.0, 1.0])
        self.coolant_transport = SimpleNamespace(
            get_cp=lambda T, p: CP,
            get_rho=lambda T, p: 1_000.0,
            get_mu=lambda T, p: 1.0e-3,
            get_k=lambda T, p: 0.5,
        )

    def A_coolant(self, x):
        return 1.0e-3

    def dA_dx_thermal_exhaust(self, x):
        return 1.0


def _manufactured_chamber(direction):
    return SimpleNamespace(
        cooling_circuits=[_Circuit(direction)],
        combustion_transport=SimpleNamespace(
            get_T=lambda x: T_GAS,
            clear_cache=lambda: None,
        ),
    )


def _install_manufactured_physics(monkeypatch) -> None:
    helper = coupled_solver.CoupledHeatExchangerPhysics
    monkeypatch.setattr(
        coupled_solver, "FluidStateResolver", _FluidStateResolver
    )
    monkeypatch.setattr(helper, "bulk_density", lambda self, *args: 1_000.0)
    monkeypatch.setattr(
        helper,
        "hot_side_coefficients",
        lambda self, x, T_hw, **kwargs: {
            "h_hot": H_HOT,
            "h_g": H_HOT,
            "T_aw": T_GAS,
            "qpp_hot": H_HOT * (T_GAS - T_hw),
            "regime": "gas",
            "q_rad": 0.0,
        },
    )
    monkeypatch.setattr(
        helper,
        "cold_side_coefficients",
        lambda self, *args, **kwargs: {"h_cold": H_COLD},
    )
    monkeypatch.setattr(
        helper,
        "dQ_hot_dx",
        lambda self, x, T_hw: H_HOT * (T_GAS - T_hw),
    )
    monkeypatch.setattr(
        helper,
        "dQ_cond_dx",
        # Manufactured zero-resistance wall: the extra term makes a nonzero
        # hot-to-cold wall jump produce a residual while retaining the exact
        # hot-side heat rate at the physical dT=0 solution.
        lambda self, x, T_hw, T_cw: (
            H_HOT * (T_GAS - T_hw) + 200.0 * (T_hw - T_cw)
        ),
    )
    monkeypatch.setattr(
        helper,
        "dQ_cold_dx",
        lambda self, x, T_cw, T_cool, **kwargs: H_COLD * (T_cw - T_cool),
    )
    monkeypatch.setattr(
        helper,
        "coolant_friction_rate",
        lambda self, *args, **kwargs: FRICTION_GRADIENT,
    )
    monkeypatch.setattr(
        helper,
        "interface_temperatures",
        lambda self, x, T_hw, T_cw: [T_hw, T_cw],
    )


@pytest.mark.parametrize("direction", [1, -1])
def test_full_regen_solver_matches_constant_property_heat_exchanger(
    monkeypatch,
    direction,
) -> None:
    _install_manufactured_physics(monkeypatch)
    nodes = 11
    inlet_temperature = 300.0
    inlet_pressure = 2.0e6
    boundary = coupled_solver.BoundaryConditions(
        T_coolant_in=inlet_temperature,
        p_coolant_in=inlet_pressure,
        mdot_coolant=MDOT,
    )

    with pytest.warns(RuntimeWarning, match="manufactured constant-property"):
        result = coupled_solver.coupled_steady_heating_analysis(
            _manufactured_chamber(direction),
            boundary_conditions=boundary,
            nodes=nodes,
            output=False,
            heat_curvature_correction=False,
            pressure_curvature_correction=False,
        )

    overall_conductance = 1.0 / (1.0 / H_HOT + 1.0 / H_COLD)
    step = 1.0 / (nodes - 1)
    expected_outlet_temperature = T_GAS - (T_GAS - inlet_temperature) * (
        1.0 - overall_conductance * step / (MDOT * CP)
    ) ** (nodes - 1)

    assert result.T_static[-1] == pytest.approx(
        expected_outlet_temperature, rel=2e-8
    )
    assert result.p_stagnation[-1] == pytest.approx(
        inlet_pressure + FRICTION_GRADIENT
    )
    assert np.all(np.diff(result.T_static) > 0.0)
    assert np.all(np.diff(result.p_stagnation) < 0.0)
    np.testing.assert_allclose(result.T[:, -1], result.T[:, 1], atol=1e-7)
    assert np.all(result.T[:, 1] > result.T[:, 0])
    assert np.all(result.wall_converged)


@pytest.mark.parametrize("direction", [1, -1])
@pytest.mark.parametrize("curvature_enabled", [False, True])
@pytest.mark.parametrize("roughness", [0.0, 2e-5], ids=["smooth", "rough"])
@pytest.mark.parametrize("mdot", [0.4, 10.0], ids=["Re4000", "Re100000"])
def test_full_regen_solver_matches_fluids_pressure_march(
    monkeypatch,
    direction,
    curvature_enabled,
    roughness,
    mdot,
) -> None:
    helper = coupled_solver.CoupledHeatExchangerPhysics
    real_friction_rate = helper.coolant_friction_rate
    _install_manufactured_physics(monkeypatch)
    monkeypatch.setattr(helper, "coolant_friction_rate", real_friction_rate)

    chamber = _manufactured_chamber(direction)
    circuit = chamber.cooling_circuits[0]
    diameter, radius, arc_length_scale = 0.01, 0.1, 1.75
    circuit.Dh_coolant = lambda x: diameter
    circuit.radius_of_curvature = lambda x: radius
    circuit.roughness = roughness
    circuit.ds_dx = lambda x: arc_length_scale
    nodes = 11
    inlet_pressure = 2.0e6
    boundary = coupled_solver.BoundaryConditions(
        T_coolant_in=300.0,
        p_coolant_in=inlet_pressure,
        mdot_coolant=mdot,
    )

    with pytest.warns(RuntimeWarning, match="manufactured constant-property"):
        result = coupled_solver.coupled_steady_heating_analysis(
            chamber,
            boundary_conditions=boundary,
            nodes=nodes,
            output=False,
            heat_curvature_correction=False,
            pressure_curvature_correction=curvature_enabled,
            friction_model="fluids",
        )

    density, viscosity, area = 1_000.0, 1e-3, 1e-3
    velocity = mdot / (density * area)
    reynolds = density * velocity * diameter / viscosity
    # The oracle evaluates the branch correlations and smoothstep independently
    # rather than inheriting the generic fluids dispatcher's abrupt transition.
    if curvature_enabled:
        factor = _curved_factor(reynolds, diameter, radius, roughness)
    else:
        factor = _straight_factor(reynolds, diameter, roughness)
    dynamic_head = 0.5 * density * velocity**2
    gradient = -factor / diameter * dynamic_head * arc_length_scale
    distance = np.linspace(0.0, 1.0, nodes)
    expected_x = distance if direction == 1 else 1.0 - distance
    expected_pressure = inlet_pressure + gradient * distance

    np.testing.assert_allclose(result.x_coolant, expected_x, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(
        result.p_stagnation, expected_pressure, rtol=0.0, atol=1e-7,
    )
    # Constant density and area eliminate acceleration; static pressure differs
    # from stagnation pressure by the same dynamic head at every station.
    np.testing.assert_allclose(
        result.p_static, expected_pressure - dynamic_head, rtol=0.0, atol=1e-7,
    )
    assert result.p_stagnation[-1] == pytest.approx(
        inlet_pressure + gradient, rel=0.0, abs=1e-7,
    )
    np.testing.assert_allclose(result.velocity, velocity, rtol=1e-12)
    assert np.all(np.diff(result.p_stagnation) < 0.0)
    assert np.all(np.diff(result.T_static) > 0.0)
    assert np.all(result.wall_converged)
