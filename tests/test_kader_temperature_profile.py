"""Kader thermal-boundary-layer reconstruction and its visualization."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen import physics
from pyskyfire.viz import PlotTemperatureProfile


@pytest.mark.parametrize("prandtl", [0.01, 0.7, 1.0, 10.0, 100.0])
def test_kader_profile_meets_wall_and_edge_conditions(prandtl) -> None:
    distance, temperature, delta = physics.kader_temperature_profile(
        T_wall=400.0,
        T_edge=1200.0,
        h=5_000.0,
        rho=2.0,
        cp=2_000.0,
        mu=4.0e-5,
        prandtl=prandtl,
        u_tau=10.0,
        n_points=101,
    )

    assert distance[0] == 0.0
    assert distance[-1] == pytest.approx(delta)
    assert delta > 0.0
    assert temperature[0] == 400.0
    assert temperature[-1] == 1200.0
    assert np.all(np.diff(distance) > 0.0)
    assert np.all(np.diff(temperature) >= 0.0)


def test_kader_conductive_sublayer_has_prandtl_slope() -> None:
    y_plus = 1.0e-4
    prandtl = 7.0

    theta_plus = physics.kader_temperature_plus(y_plus, prandtl, 0.0)

    assert theta_plus == pytest.approx(prandtl * y_plus, rel=1.0e-10)


def test_friction_velocity_uses_darcy_convention() -> None:
    assert physics.friction_velocity(20.0, 0.02) == pytest.approx(
        20.0 * np.sqrt(0.02 / 8.0)
    )


class _Coolant:
    def get_rho(self, temperature, pressure):
        return 70.0

    def get_mu(self, temperature, pressure):
        return 1.0e-5

    def get_cp(self, temperature, pressure):
        return 10_000.0

    def get_Pr(self, temperature, pressure):
        return 1.2


class _Circuit:
    def __init__(self):
        self.coolant_transport = _Coolant()
        self.walls = [
            SimpleNamespace(
                material=SimpleNamespace(name="Silver"),
                thickness=lambda x: 2.0e-4,
            ),
            SimpleNamespace(
                material=SimpleNamespace(name="Steel"),
                thickness=lambda x: 3.0e-4,
            ),
        ]

    def Dh_coolant(self, x):
        return 0.01

    def roughness(self, x):
        return 1.0e-6

    def radius_of_curvature(self, x):
        return np.inf


def test_temperature_plot_uses_recovery_temperature_and_material_colors() -> None:
    gas = SimpleNamespace(
        T=1500.0,
        M=0.3,
        a=1_000.0,
        rho=2.0,
        cp=2_000.0,
        mu=4.0e-5,
        Pr=0.7,
    )
    circuit = _Circuit()
    chamber = SimpleNamespace(
        combustion_transport=SimpleNamespace(get_state=lambda x: gas),
        contour=SimpleNamespace(r=lambda x: 0.1),
        cooling_circuits=[circuit],
    )
    result = SimpleNamespace(
        x=np.array([0.0]),
        x_heat_flux=np.array([0.0]),
        x_coolant=np.array([0.0]),
        T=np.array([[100.0, 200.0, 300.0, 800.0]]),
        T_static=np.array([100.0]),
        p_static=np.array([1.0e6]),
        velocity=np.array([20.0]),
        h_hot=np.array([10_000.0]),
        h_cold=np.array([5_000.0]),
        T_aw_hot=np.array([2_000.0]),
    )

    plot = PlotTemperatureProfile(result, chamber, 0, 0.0, n_bl=31)

    profile = np.asarray(plot.fig.data[2].y, dtype=float)
    assert profile[0] == pytest.approx(2_000.0)
    assert np.max(profile) == pytest.approx(2_000.0)
    wall_shapes = plot.fig.layout.shapes[1:3]
    assert wall_shapes[0].fillcolor != wall_shapes[1].fillcolor
    assert [shape.name for shape in wall_shapes] == ["Silver", "Steel"]
    assert all(shape.showlegend for shape in wall_shapes)
