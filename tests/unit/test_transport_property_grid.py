"""Tests for choosing the visualization grid of transport-property plots."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.viz.plot_skycea import (
    PlotTransportProperty,
    _visualization_x,
    result_gas_x,
)


class _Transport:
    """Live-solve transport stub that counts the stations it is asked for."""

    def __init__(self, contour=None):
        self.contour = contour
        self.queried = []

    def get_T(self, x):
        self.queried.append(float(x))
        return 3000.0 - 100.0 * float(x)


def _contour(x_min=-0.2, x_max=1.0):
    return SimpleNamespace(xs=np.array([x_min, x_max]))


def test_results_grid_uses_solved_wall_and_heat_flux_stations() -> None:
    result = SimpleNamespace(
        x_wall=np.array([0.0, 0.5, 1.0]),
        x_heat_flux=np.array([0.25, 0.5]),
        x=np.array([9.0]),
    )

    assert result_gas_x(result) == pytest.approx([0.0, 0.25, 0.5, 1.0])


def test_results_grid_falls_back_to_x_and_merges_several_circuits() -> None:
    circuit_a = {"x": np.array([0.0, 0.5])}
    circuit_b = {"x": np.array([0.5, 1.0])}

    assert result_gas_x((circuit_a, circuit_b)) == pytest.approx([0.0, 0.5, 1.0])


def test_results_dict_of_circuits_is_expanded() -> None:
    results = {
        "regen_half_pass": {"x": np.array([0.0, 0.5])},
        "regen_full_pass": {"x": np.array([1.0])},
    }

    assert result_gas_x(results) == pytest.approx([0.0, 0.5, 1.0])


def test_nodes_builds_a_uniform_grid_over_the_contour() -> None:
    transport = _Transport(_contour(0.0, 1.0))

    grid = _visualization_x(transport, x=None, results=None, nodes=5)

    assert grid == pytest.approx([0.0, 0.25, 0.5, 0.75, 1.0])


def test_explicit_x_is_used_verbatim() -> None:
    transport = _Transport(_contour())

    grid = _visualization_x(transport, x=[0.3, 0.1], results=None, nodes=None)

    assert grid == pytest.approx([0.3, 0.1])


def test_only_one_grid_source_is_accepted() -> None:
    transport = _Transport(_contour())

    with pytest.raises(ValueError, match="only one of x, results, or nodes"):
        _visualization_x(transport, x=[0.0], results=None, nodes=4)


def test_missing_grid_reports_the_available_choices() -> None:
    transport = _Transport(_contour())

    with pytest.raises(ValueError, match="results=..., nodes=..., or x="):
        _visualization_x(transport, x=None, results=None, nodes=None)


def test_legacy_precomputed_objects_keep_working() -> None:
    legacy = SimpleNamespace(x_nodes=np.array([0.0, 1.0]))

    grid = _visualization_x(legacy, x=None, results=None, nodes=None)

    assert grid == pytest.approx([0.0, 1.0])


def test_plot_samples_the_result_grid() -> None:
    transport = _Transport(_contour())
    result = SimpleNamespace(x_wall=np.array([0.0, 0.5, 1.0]))

    plot = PlotTransportProperty(transport, prop="T", results=result)

    assert transport.queried == pytest.approx([0.0, 0.5, 1.0])
    trace = plot.fig.data[0]
    assert np.asarray(trace.x) == pytest.approx([0.0, 0.5, 1.0])
    assert np.asarray(trace.y) == pytest.approx([3000.0, 2950.0, 2900.0])
