"""Cooling-channel centerline and RL10/RTE curvature corrections."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen import physics
from pyskyfire.regen import coupled_solver
from pyskyfire.regen.channel_placement import SurfacePlacement
from pyskyfire.regen.cross_section import ChannelSection
from pyskyfire.regen.coupled_solver import CoupledHeatExchangerPhysics
from pyskyfire.regen.thrust_chamber import (
    CoolingCircuit,
    ThrustChamber,
    Wall,
    cylindrical_curve_arc_length_scale,
    radius_of_curvature,
)


def test_straight_centerline_has_infinite_radius() -> None:
    x = np.linspace(-1.0, 1.0, 151)
    points = np.column_stack([x, np.full_like(x, 2.0), np.zeros_like(x)])

    assert np.all(np.isinf(radius_of_curvature(points)))


def test_irregularly_sampled_circle_has_constant_radius() -> None:
    expected_radius = 0.4
    angle = np.array([-0.7, -0.41, -0.18, 0.0, 0.09, 0.31, 0.62])
    points = np.column_stack(
        [
            expected_radius * np.sin(angle),
            1.0 + expected_radius * (1.0 - np.cos(angle)),
            np.zeros_like(angle),
        ]
    )

    np.testing.assert_allclose(
        radius_of_curvature(points),
        expected_radius,
        rtol=1e-12,
        atol=1e-12,
    )


def test_constant_radius_helix_uses_full_3d_curvature() -> None:
    x = np.linspace(0.0, 2.0, 1001)
    radius = 2.0
    dtheta_dx = 0.7
    points = np.column_stack(
        [x, np.full_like(x, radius), dtheta_dx * x]
    )

    computed_radius = radius_of_curvature(points)
    expected_radius = -(
        1.0 + (radius * dtheta_dx) ** 2
    ) / (radius * dtheta_dx**2)
    expected_ds_dx = np.sqrt(1.0 + (radius * dtheta_dx) ** 2)

    np.testing.assert_allclose(
        computed_radius[10:-10],
        expected_radius,
        rtol=3e-6,
    )
    np.testing.assert_allclose(
        cylindrical_curve_arc_length_scale(points)[10:-10],
        expected_ds_dx,
        rtol=1e-10,
    )


@pytest.mark.parametrize("signed_radius, exponent", [(0.05, 0.05), (-0.05, -0.05)])
def test_rl10_heat_factor_preserves_bend_direction(
    signed_radius: float,
    exponent: float,
) -> None:
    reynolds = 100_000.0
    diameter = 0.01
    group = reynolds * (0.5 * diameter / abs(signed_radius)) ** 2

    assert physics.phi_curv(reynolds, diameter, signed_radius) == pytest.approx(
        group**exponent
    )


def test_curvature_factors_have_straight_and_applicability_limits() -> None:
    assert physics.phi_curv(100_000.0, 0.01, np.inf) == 1.0
    assert physics.phi_curv(100_000.0, 0.01, 100.0) == 1.0
    assert physics.phi_curv_friction(100_000.0, 0.01, np.inf) == 1.0

    # RTE Eq. 22 is inactive at group <= 6 and active above it.
    diameter = 0.01
    reynolds = 100_000.0
    radius_at_six = 0.25 * diameter * np.sqrt(reynolds / 6.0)
    assert physics.phi_curv_friction(reynolds, diameter, radius_at_six) == 1.0

    active_radius = 0.01
    group = reynolds * (0.25 * diameter / active_radius) ** 2
    assert physics.phi_curv_friction(
        reynolds, diameter, active_radius
    ) == pytest.approx(group**0.05)


class _Contour:
    xs = np.array([-1.0, 0.0, 1.0])
    x_t = 0.0

    @staticmethod
    def r(x):
        return 2.0 + 0.1 * x

    @staticmethod
    def dr_dx(x):
        return 0.1


class _CrossSection(ChannelSection):
    @staticmethod
    def P_thermal(profiles):
        return np.ones_like(profiles.h)

    @staticmethod
    def P_coolant(profiles):
        return np.ones_like(profiles.h)

    @staticmethod
    def A_coolant(profiles):
        return np.ones_like(profiles.h)

    @staticmethod
    def Dh_coolant(profiles):
        return np.ones_like(profiles.h)


class _Material:
    @staticmethod
    def get_k(temperature):
        return 1.0


def test_surface_helix_follows_coolant_centroid_and_requested_angle() -> None:
    gamma = np.deg2rad(30.0)
    wall_thickness = 0.02
    channel_height = 0.04
    circuit = CoolingCircuit(
        name="test",
        contour=_Contour(),
        cross_section=_CrossSection(),
        span=[-1.0, 1.0],
        placement=SurfacePlacement(8, helix_angle=gamma),
        channel_height=lambda x: channel_height,
        walls=[Wall(_Material(), wall_thickness)],
        coolant_transport=SimpleNamespace(),
        roughness=0.0,
    )

    ThrustChamber(
        contour=_Contour(),
        cooling_circuits=[circuit],
        combustion_transport=SimpleNamespace(),
        n_nodes=101,
        compute_gas=False,
    )

    flow_centerline = circuit.coolant_centerlines[0]
    expected_radius = _Contour.r(flow_centerline[:, 0]) + wall_thickness + 0.5 * channel_height
    expected_meridional_scale = np.sqrt(1.0 + 0.1**2)

    np.testing.assert_allclose(flow_centerline[:, 1], expected_radius)
    np.testing.assert_allclose(
        circuit.ds_dx_vals,
        expected_meridional_scale / np.cos(gamma),
        rtol=2e-7,
    )
    assert circuit.is_helical


def test_cold_side_colburn_receives_curvature_factor(monkeypatch) -> None:
    placement = SimpleNamespace(n_channel_positions=2, n_channels_per_leaf=1)
    circuit = SimpleNamespace(
        placement=placement,
        n_channels=2,
        coolant_transport=SimpleNamespace(
            get_k=lambda temperature, pressure: 1.0,
            get_cp=lambda temperature, pressure: 2.0,
            get_mu=lambda temperature, pressure: 2.0,
        ),
        Dh_coolant=lambda x: 0.01,
        A_coolant=lambda x: 0.001,
        radius_of_curvature=lambda x: 0.05,
    )
    helper = object.__new__(CoupledHeatExchangerPhysics)
    helper.thrust_chamber = SimpleNamespace(
        cooling_circuits=[circuit],
        h_cold_corr=1.0,
    )
    helper.circuit_index = 0
    helper.boundary_conditions = SimpleNamespace(mdot_coolant=1.0)
    helper.coolant_state = SimpleNamespace(
        film_properties=lambda *args: (1.0, 2.0, 3.0)
    )
    helper.bulk_density = lambda *args: 4.0
    helper.heat_curvature_correction = True

    monkeypatch.setattr(physics, "phi_curv", lambda *args: 1.23)
    monkeypatch.setattr(
        physics,
        "h_coolant_colburn",
        lambda *args, phi_curv: phi_curv,
    )

    result = helper.cold_side_coefficients(
        0.0,
        T_cw=300.0,
        T_cool=100.0,
        p_cool=1.0e6,
    )

    assert result["phi_curv"] == 1.23
    assert result["h_cold"] == 1.23


def test_heat_curvature_correction_can_be_disabled(monkeypatch) -> None:
    placement = SimpleNamespace(n_channel_positions=2, n_channels_per_leaf=1)
    circuit = SimpleNamespace(
        placement=placement,
        n_channels=2,
        coolant_transport=SimpleNamespace(
            get_k=lambda temperature, pressure: 1.0,
            get_cp=lambda temperature, pressure: 2.0,
            get_mu=lambda temperature, pressure: 2.0,
        ),
        Dh_coolant=lambda x: 0.01,
        A_coolant=lambda x: 0.001,
    )
    helper = object.__new__(CoupledHeatExchangerPhysics)
    helper.thrust_chamber = SimpleNamespace(
        cooling_circuits=[circuit],
        h_cold_corr=1.0,
    )
    helper.circuit_index = 0
    helper.boundary_conditions = SimpleNamespace(mdot_coolant=1.0)
    helper.coolant_state = SimpleNamespace(
        film_properties=lambda *args: (1.0, 2.0, 3.0)
    )
    helper.bulk_density = lambda *args: 4.0
    helper.heat_curvature_correction = False

    monkeypatch.setattr(
        physics,
        "phi_curv",
        lambda *args: pytest.fail("disabled heat correction was evaluated"),
    )
    monkeypatch.setattr(
        physics,
        "h_coolant_colburn",
        lambda *args, phi_curv: phi_curv,
    )

    result = helper.cold_side_coefficients(
        0.0,
        T_cw=300.0,
        T_cool=100.0,
        p_cool=1.0e6,
    )

    assert result["phi_curv"] == 1.0
    assert result["h_cold"] == 1.0


def test_pressure_loss_multiplies_darcy_factor_by_curvature(monkeypatch) -> None:
    placement = SimpleNamespace(n_channel_positions=2, n_channels_per_leaf=1)
    circuit = SimpleNamespace(
        placement=placement,
        n_channels=2,
        coolant_transport=SimpleNamespace(get_mu=lambda temperature, pressure: 1.0),
        A_coolant=lambda x: 4.0,
        Dh_coolant=lambda x: 0.5,
        radius_of_curvature=lambda x: 0.25,
        roughness=0.0,
        ds_dx=lambda x: 2.0,
        dA_dx_coolant=lambda x: 0.0,
    )
    helper = object.__new__(CoupledHeatExchangerPhysics)
    helper.thrust_chamber = SimpleNamespace(cooling_circuits=[circuit])
    helper.circuit_index = 0
    helper.boundary_conditions = SimpleNamespace(mdot_coolant=4.0)
    helper.coolant_state = SimpleNamespace(saturation=lambda pressure: None)
    helper.bulk_density = lambda *args: 10.0
    helper.pressure_curvature_correction = True

    monkeypatch.setattr(physics, "f_darcy", lambda *args: 2.0)
    monkeypatch.setattr(physics, "phi_curv_friction", lambda *args: 3.0)

    dp_friction_dx = helper.coolant_friction_rate(
        0.0,
        T_cool=100.0,
        p_cool=1.0e6,
    )

    velocity = (4.0 / 2.0) / (10.0 * 4.0)
    expected = -2.0 * 3.0 / 0.5 * 10.0 * velocity**2 / 2.0 * 2.0
    assert dp_friction_dx == pytest.approx(expected)


def test_pressure_curvature_correction_can_be_disabled(monkeypatch) -> None:
    placement = SimpleNamespace(n_channel_positions=2, n_channels_per_leaf=1)
    circuit = SimpleNamespace(
        placement=placement,
        n_channels=2,
        coolant_transport=SimpleNamespace(get_mu=lambda temperature, pressure: 1.0),
        A_coolant=lambda x: 4.0,
        Dh_coolant=lambda x: 0.5,
        roughness=0.0,
        ds_dx=lambda x: 2.0,
        dA_dx_coolant=lambda x: 0.0,
    )
    helper = object.__new__(CoupledHeatExchangerPhysics)
    helper.thrust_chamber = SimpleNamespace(cooling_circuits=[circuit])
    helper.circuit_index = 0
    helper.boundary_conditions = SimpleNamespace(mdot_coolant=4.0)
    helper.coolant_state = SimpleNamespace(saturation=lambda pressure: None)
    helper.bulk_density = lambda *args: 10.0
    helper.pressure_curvature_correction = False

    monkeypatch.setattr(physics, "f_darcy", lambda *args: 2.0)
    monkeypatch.setattr(
        physics,
        "phi_curv_friction",
        lambda *args: pytest.fail("disabled pressure correction was evaluated"),
    )

    dp_friction_dx = helper.coolant_friction_rate(
        0.0,
        T_cool=100.0,
        p_cool=1.0e6,
    )

    velocity = (4.0 / 2.0) / (10.0 * 4.0)
    expected = -2.0 / 0.5 * 10.0 * velocity**2 / 2.0 * 2.0
    assert dp_friction_dx == pytest.approx(expected)


def test_public_analysis_forwards_curvature_switches(monkeypatch) -> None:
    captured = {}
    sentinel = object()

    def fake_solve(*args, **kwargs):
        captured.update(kwargs)
        return sentinel

    monkeypatch.setattr(coupled_solver, "solve_coupled_heat_exchanger", fake_solve)
    chamber = SimpleNamespace(film_cooling=None)

    result = coupled_solver.coupled_steady_heating_analysis(
        chamber,
        SimpleNamespace(),
        output=False,
        heat_curvature_correction=False,
        pressure_curvature_correction=True,
    )

    assert result is sentinel
    assert captured["heat_curvature_correction"] is False
    assert captured["pressure_curvature_correction"] is True
