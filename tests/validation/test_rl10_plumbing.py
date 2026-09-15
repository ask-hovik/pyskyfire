"""Check the RL10 plumbing geometry against its own drawings.

A bending table over-specifies the tube: the tangent-distance column is
derivable from the radius and angle columns, and the drawing's overall length
is derivable from all of them. Those redundancies are the only independent
check available on numbers read off a figure, so they are asserted here.
"""

import importlib.util
import math
from pathlib import Path
import sys

import numpy as np
import pytest


PLUMBING_DIR = (
    Path(__file__).resolve().parents[2] / "validation" / "RL10" / "plumbing"
)


def _load_plumbing():
    """Import the validation plumbing package without installing it."""
    module_name = "pyskyfire_test_rl10_plumbing"
    if module_name in sys.modules:
        return sys.modules[module_name]

    spec = importlib.util.spec_from_file_location(
        module_name,
        PLUMBING_DIR / "__init__.py",
        submodule_search_locations=[str(PLUMBING_DIR)],
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load RL10 plumbing from {PLUMBING_DIR}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def duct():
    return _load_plumbing().turbine_discharge_duct


def test_bore_follows_from_the_drawing_callouts(duct) -> None:
    """OD less twice the wall must equal Section A-A's average diameter."""
    assert duct.OUTER_DIAMETER - 2.0 * duct.WALL_THICKNESS == pytest.approx(
        duct.BORE_DIAMETER, abs=5.0e-4
    )


def test_resolved_tangent_distances_match_figure_g1(duct) -> None:
    """``R tan(angle / 2)`` must reproduce Figure G1's own column."""
    from_table = _load_plumbing().tangent_distances(duct.BENDS)
    assert from_table == pytest.approx(duct.TABULATED_TANGENT_DISTANCES, abs=1.0e-3)


def test_centerline_length_matches_table_g1(duct) -> None:
    """Straights plus bend arcs must come out at the stated 51.7 in duct."""
    arcs = sum(
        radius * math.radians(angle) for radius, angle in duct.BENDS
    )
    assert sum(duct.STRAIGHTS) + arcs == pytest.approx(
        duct.REFERENCE_DUCT_LENGTH, abs=0.05
    )


def test_resolved_straights_come_back_out_at_their_tabulated_lengths(duct) -> None:
    """Each corner must absorb its tangent distance from both neighbours.

    The route vertices sit at leg intersections, so the straight between two
    consecutive vertices is the tabulated straight plus the tangent distances
    of the bends at either end.
    """
    inch = _load_plumbing().INCH
    vertices = np.asarray([node.position for node in duct.build_route()]) / inch
    legs = np.linalg.norm(np.diff(vertices, axis=0), axis=1)

    tangents = _load_plumbing().tangent_distances(duct.BENDS)
    padded = (0.0, *tangents, 0.0)
    expected = [
        straight + padded[index] + padded[index + 1]
        for index, straight in enumerate(duct.STRAIGHTS)
    ]
    assert legs == pytest.approx(expected, abs=1.0e-6)


def test_pipe_resolves_one_lumped_bend_per_corner(duct) -> None:
    """Each bend must become one connection between its two tangent points."""
    pipe = duct.build_pipe()
    kinds = [type(connection).__name__ for connection in pipe.connections]

    assert kinds == ["Straight", "RoundedBend"] * len(duct.BENDS) + ["Straight"]
    assert len(pipe.stations) == 2 * len(duct.BENDS) + 2
    assert pipe.roughness == pytest.approx(duct.ROUGHNESS)
    assert all(
        diameter == pytest.approx(duct.BORE_DIAMETER * _load_plumbing().INCH)
        for diameter in pipe.diameter
    )


def test_drawn_arc_length_matches_table_g1(duct) -> None:
    """The reopened bend arcs must recover the stated length, not the chords.

    ``Pipe.centerline`` cuts each corner as a chord, because a lumped bend is
    one connection between tangent points. ``display_path`` reopens those into
    arcs, so its polyline converges on the drawing's 51.7 in from below.
    """
    inch = _load_plumbing().INCH
    points = duct.build_pipe().display_path(segments=256).points
    length = np.linalg.norm(np.diff(points, axis=0), axis=1).sum() / inch

    assert length == pytest.approx(duct.REFERENCE_DUCT_LENGTH, abs=0.05)
