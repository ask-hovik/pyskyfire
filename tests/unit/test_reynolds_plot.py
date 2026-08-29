from types import SimpleNamespace

import numpy as np

from pyskyfire.regen import physics
from pyskyfire.viz.plot_regen import PlotReynoldsNumber


class _CombustionTransport:
    x_nodes = np.array([0.0, 1.0])

    def get_state(self, x):
        return SimpleNamespace(rho=1.0, M=1.0, a=100.0, mu=0.01)


class _CoolantTransport:
    def get_rho(self, temperature, pressure):
        return 1.0

    def get_mu(self, temperature, pressure):
        return 1.0


def _chamber(diameter=1.0):
    circuit = SimpleNamespace(
        coolant_transport=_CoolantTransport(),
        Dh_coolant=lambda x: diameter,
    )
    return SimpleNamespace(
        combustion_transport=_CombustionTransport(),
        contour=SimpleNamespace(r=lambda x: 0.5),
        cooling_circuits=[circuit],
    )


def _result(velocity):
    velocity = np.asarray(velocity, dtype=float)
    return SimpleNamespace(
        circuit_name="Test circuit",
        circuit_index=0,
        x=np.array([0.0, 1.0]),
        x_coolant=np.array([0.0, 1.0]),
        T_static=np.ones(2),
        p_stagnation=np.ones(2),
        velocity=velocity,
        coolant_quality=None,
    )


def test_reynolds_plot_contains_combustion_and_coolant_profiles():
    plot = PlotReynoldsNumber(_chamber(), _result([10.0, 20.0]))

    assert [trace.name for trace in plot.fig.data] == [
        "Combustion gas",
        "Test circuit coolant",
    ]
    np.testing.assert_allclose(plot.fig.data[0].y, [10_000.0, 10_000.0])
    np.testing.assert_allclose(plot.fig.data[1].y, [10.0, 20.0])


def test_reynolds_plot_labels_only_thresholds_crossed_by_a_profile():
    plot = PlotReynoldsNumber(
        _chamber(),
        _result([physics.ReDh_laminar - 1, physics.ReDh_turbulent + 1]),
        combustion=False,
    )

    assert [shape.y0 for shape in plot.fig.layout.shapes] == [
        physics.ReDh_laminar,
        physics.ReDh_turbulent,
    ]
    assert [annotation.text for annotation in plot.fig.layout.annotations] == [
        f"Laminar limit (Re = {physics.ReDh_laminar:g})",
        f"Turbulent limit (Re = {physics.ReDh_turbulent:g})",
    ]


def test_reynolds_plot_omits_uncrossed_thresholds():
    plot = PlotReynoldsNumber(
        _chamber(),
        _result([10_000.0, 20_000.0]),
        combustion=False,
    )

    assert not plot.fig.layout.shapes
    assert not plot.fig.layout.annotations
