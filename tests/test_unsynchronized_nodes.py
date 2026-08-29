"""Node-grid normalization for the coupled regenerative solver."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen import coupled_solver


def circuit(direction=1):
    return SimpleNamespace(
        direction=direction,
        x_domain=np.array([-1.0, 1.0]),
    )


def test_integer_nodes_expand_to_three_internal_grids() -> None:
    grids = coupled_solver._prepare_node_grids(4, circuit(direction=-1))

    assert isinstance(grids, list)
    assert len(grids) == 3
    for grid in grids:
        np.testing.assert_allclose(grid, [1.0, 1.0 / 3.0, -1.0 / 3.0, -1.0])
    assert grids[0] is not grids[1]


def test_ragged_mapping_is_sorted_in_coolant_march_order() -> None:
    grids = coupled_solver._prepare_node_grids(
        {
            "wall": [0.9, -0.9, 0.4, -0.4],
            "heat_flux": [-0.8, 0.8],
            "coolant": [-1.0, 0.0, 1.0],
        },
        circuit(direction=1),
    )

    wall, heat_flux, coolant = grids
    np.testing.assert_allclose(wall, [-0.9, -0.4, 0.4, 0.9])
    np.testing.assert_allclose(heat_flux, [-0.8, 0.8])
    np.testing.assert_allclose(coolant, [-1.0, 0.0, 1.0])
    assert [len(grid) for grid in grids] == [4, 2, 3]


def test_each_coolant_segment_requires_a_wall_node() -> None:
    with pytest.raises(ValueError, match="segments without wall nodes"):
        coupled_solver._prepare_node_grids(
            [[-0.9, -0.8], [-0.9, 0.9], [-1.0, 0.0, 1.0]],
            circuit(direction=1),
        )


def test_public_analysis_forwards_nodes(monkeypatch) -> None:
    captured = {}
    requested_nodes = [[-0.5, 0.5], [-0.4, 0.4], [-0.6, 0.6]]

    def fake_solve(*args, **kwargs):
        captured["nodes"] = args[2]
        return object()

    monkeypatch.setattr(coupled_solver, "solve_coupled_heat_exchanger", fake_solve)
    chamber = SimpleNamespace(film_cooling=None)
    coupled_solver.coupled_steady_heating_analysis(
        chamber,
        SimpleNamespace(),
        nodes=requested_nodes,
        output=False,
    )

    assert captured["nodes"] is requested_nodes
