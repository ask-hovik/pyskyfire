"""Plot annotations for regenerative wall-balance failures."""

from types import SimpleNamespace

import numpy as np

from pyskyfire.viz.plot_regen import PlotWallTemperature


def result_with_failed_wall_nodes():
    return SimpleNamespace(
        circuit_name="Test circuit",
        x=np.array([0.0, 0.5, 1.0]),
        x_wall=np.array([0.0, 0.5, 1.0]),
        T=np.array(
            [
                [100.0, 200.0, 300.0],
                [110.0, 210.0, 310.0],
                [120.0, 220.0, 320.0],
            ]
        ),
        wall_residual_scaled=np.array([1.0e-8, 2.0e-4, 3.0e-3]),
        wall_converged=np.array([True, False, False]),
    )


def test_wall_plot_marks_only_failed_nodes_with_x_symbols() -> None:
    plot = PlotWallTemperature(result_with_failed_wall_nodes())

    assert len(plot.fig.data) == 2
    marker_trace = plot.fig.data[1]
    assert marker_trace.name == "Non-converged wall nodes"
    assert marker_trace.mode == "markers"
    assert marker_trace.marker.symbol == "x"
    np.testing.assert_allclose(marker_trace.x, [0.5, 1.0])
    np.testing.assert_allclose(marker_trace.y, [310.0, 320.0])
    np.testing.assert_allclose(marker_trace.customdata, [2.0e-4, 3.0e-3])


def test_wall_plot_markers_can_be_disabled() -> None:
    plot = PlotWallTemperature(
        result_with_failed_wall_nodes(),
        mark_nonconverged=False,
    )

    assert len(plot.fig.data) == 1
