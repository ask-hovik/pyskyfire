"""RL10 silver-insert geometry and simulation wiring tests."""

import ast
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import CubicSpline, PchipInterpolator

import pyskyfire as psf


ROOT = Path(__file__).resolve().parents[1]
RL10_DIR = ROOT / "validation" / "RL10_2"
REFERENCE_DIR = RL10_DIR / "reference_data"
REGEN_SIM_PATH = RL10_DIR / "regen_sim.py"


def _load_geometry_functions():
    """Load pure geometry helpers without running the RL10 simulation setup."""
    tree = ast.parse(REGEN_SIM_PATH.read_text())
    names = {
        "silver_insert_thickness_from_json",
        "contour_from_area_ratios",
        "align_positive_channel_height_x",
        "channel_height_from_json",
    }
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    namespace = {
        "json": json,
        "np": np,
        "CubicSpline": CubicSpline,
        "PchipInterpolator": PchipInterpolator,
        "psf": psf,
    }
    exec(
        compile(
            ast.Module(body=functions, type_ignores=[]),
            str(REGEN_SIM_PATH),
            "exec",
        ),
        namespace,
    )
    return namespace


def test_silver_insert_thickness_is_positive_only_inside_insert() -> None:
    functions = _load_geometry_functions()
    thickness = functions["silver_insert_thickness_from_json"](
        REFERENCE_DIR / "reference_silver_insert.json"
    )
    x_start, x_end = thickness.x_span

    assert thickness.throat_radius == pytest.approx(0.062776)
    assert thickness(x_start - 1.0e-6) == 0.0
    assert thickness(x_end + 1.0e-6) == 0.0
    assert thickness(x_start) == pytest.approx(0.0, abs=1.0e-15)
    assert thickness(x_end) == pytest.approx(0.0, abs=1.0e-15)
    assert thickness(0.0) > 0.0
    assert np.all(thickness(np.linspace(x_start, x_end, 501)) >= 0.0)


def test_only_insert_throat_radius_is_added_to_area_ratio_contour() -> None:
    functions = _load_geometry_functions()
    thickness = functions["silver_insert_thickness_from_json"](
        REFERENCE_DIR / "reference_silver_insert.json"
    )
    contour = functions["contour_from_area_ratios"](
        REFERENCE_DIR / "reference_area_ratios.json",
        thickness.throat_radius,
        46.0 * 0.0254,
        61.0,
    )

    assert contour.r_t == pytest.approx(thickness.throat_radius)
    assert contour.chamber_angle is None
    assert contour.nozzle_angle is None
    assert np.count_nonzero(contour.input_xs == 0.0) == 1
    assert contour.input_xs.size == 18
    # The coarse Table E1 near-throat stations remain in place.
    assert -0.0127 in contour.input_xs
    assert 0.012319 in contour.input_xs
    # Non-throat detail from the insert does not constrain this contour.
    assert thickness.inner_x[1] not in contour.input_xs
    assert thickness.inner_x[-2] not in contour.input_xs


def test_only_long_tube_circuit_receives_silver_layer() -> None:
    tree = ast.parse(REGEN_SIM_PATH.read_text())
    build = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "build_thrust_chamber"
    )
    assignments = {
        target.id: node.value
        for node in build.body
        if isinstance(node, ast.Assign)
        for target in node.targets
        if isinstance(target, ast.Name)
    }

    def wall_names(call):
        walls = next(keyword.value for keyword in call.keywords if keyword.arg == "walls")
        return [element.id for element in walls.elts]

    assert wall_names(assignments["half_pass"]) == ["wall"]
    assert wall_names(assignments["full_pass"]) == ["silver_wall", "wall"]


def test_hot_side_liner_is_separate_from_channel_wall_geometry() -> None:
    circuit = psf.regen.CoolingCircuit.__new__(psf.regen.CoolingCircuit)
    circuit.walls = [
        psf.regen.Wall(object(), 0.002, name="liner"),
        psf.regen.Wall(object(), 0.0003, name="tube wall"),
    ]
    circuit.channel_heights = np.array([0.003, 0.003])
    circuit.channel_local_sector = np.array([0.1, 0.1])
    centerline = np.array([[0.0, 0.06, 0.0], [0.01, 0.061, 0.0]])

    profile = circuit._prof(centerline)

    np.testing.assert_allclose(profile.t_wall, 0.0003)
    np.testing.assert_allclose(
        circuit.hot_side_layer_thickness(centerline[:, 0]),
        0.002,
    )
    np.testing.assert_allclose(
        circuit.total_thickness(centerline[:, 0]),
        0.0023,
    )


def test_channel_height_x_alignment_preserves_nozzle_length() -> None:
    functions = _load_geometry_functions()
    align = functions["align_positive_channel_height_x"]
    source_interleaving = 0.3336
    binder_interleaving = 0.24422599467379513
    source_x = np.array([-0.2, 0.0, 0.1, source_interleaving, 1.2995633187772928])

    aligned_x = align(source_x, source_interleaving, binder_interleaving)

    np.testing.assert_allclose(aligned_x[:2], source_x[:2])
    assert aligned_x[3] == pytest.approx(binder_interleaving)
    assert aligned_x[2] == pytest.approx(
        source_x[2] * binder_interleaving / source_interleaving
    )
    assert aligned_x[-1] - binder_interleaving == pytest.approx(
        source_x[-1] - source_interleaving
    )
    assert aligned_x[-1] > 46.0 * 0.0254
    assert np.all(np.diff(aligned_x) > 0.0)


def test_channel_height_x_alignment_can_hold_nozzle_exit_fixed() -> None:
    functions = _load_geometry_functions()
    align = functions["align_positive_channel_height_x"]
    source_interleaving = 0.3336
    target_interleaving = 0.24422599467379513
    source_exit = 1.2062077509969191
    target_exit = 46.0 * 0.0254
    source_x = np.array(
        [-0.2, 0.0, 0.1, source_interleaving, 0.8, source_exit]
    )

    aligned_x = align(
        source_x,
        source_interleaving,
        target_interleaving,
        source_end_x=source_exit,
        target_end_x=target_exit,
    )

    assert aligned_x[3] == pytest.approx(target_interleaving)
    assert aligned_x[-1] == pytest.approx(target_exit)
    assert np.all(np.diff(aligned_x) > 0.0)


def test_recop_od_height_is_converted_to_internal_height() -> None:
    functions = _load_geometry_functions()
    functions["contour"] = type(
        "ContourBounds",
        (),
        {"xs": np.array([-0.3048, 46.0 * 0.0254])},
    )()
    wall_thickness = 0.3302e-3
    channel_height = functions["channel_height_from_json"](
        REFERENCE_DIR / "reference_channel_height.json",
        target_interleaving_x=0.24422599467379513,
        source_interleaving_x=0.3336,
        wall_thickness=wall_thickness,
    )

    profile = channel_height.profile
    np.testing.assert_allclose(
        channel_height(channel_height.x_profile),
        np.asarray(profile["tube_height_od_m"]) - 2.0 * wall_thickness,
    )
    assert channel_height.x_profile[-1] == pytest.approx(46.0 * 0.0254)
    assert np.min(channel_height(channel_height.x_profile)) > 0.0
