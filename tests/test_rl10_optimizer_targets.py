import importlib.util
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "validation"
    / "RL10_2"
    / "old"
    / "optimize_channel_heights.py"
)
SPEC = importlib.util.spec_from_file_location("rl10_optimizer_targets", MODULE_PATH)
optimizer = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = optimizer
SPEC.loader.exec_module(optimizer)


def test_combined_target_gives_each_observable_equal_aggregate_weight():
    target = optimizer.load_fit_target("combined")
    predicted_parts = [
        component.values + component.residual_scale
        for component in target.components
    ]
    residual = optimizer.normalized_data_residual(
        target,
        np.concatenate(predicted_parts),
    )

    start = 0
    for component in target.components:
        stop = start + component.values.size
        assert np.isclose(np.linalg.norm(residual[start:stop]), 1.0)
        start = stop

    assert np.isclose(np.linalg.norm(residual), 2.0)


def test_selected_targets_are_combined_with_equal_aggregate_weight():
    target = optimizer.load_fit_targets(
        ["wall_temperature", "static_pressure"],
        "RTE_model",
    )
    predicted_parts = [
        component.values + component.residual_scale
        for component in target.components
    ]
    residual = optimizer.normalized_data_residual(
        target,
        np.concatenate(predicted_parts),
    )

    assert target.name == "wall_temperature+static_pressure"
    assert [component.name for component in target.components] == [
        "wall_temperature",
        "static_pressure",
    ]
    assert np.isclose(np.linalg.norm(residual), np.sqrt(2.0))


def test_combined_cannot_be_mixed_with_individual_targets():
    with np.testing.assert_raises_regex(ValueError, "cannot be mixed"):
        optimizer.load_fit_targets(["combined", "wall_temperature"])


def test_rte_wall_temperature_target_loads_requested_reference_model():
    target = optimizer.load_fit_target("wall_temperature", "RTE_model")

    assert target.reference_model == "RTE_model"
    assert target.values.size == 20
    assert target.values[0] == np.float64(255.55104905009307)


def test_missing_reference_model_reports_available_series():
    with np.testing.assert_raises_regex(ValueError, "choose one of:.*RTE_model"):
        optimizer.load_fit_target("wall_temperature", "missing_model")


def test_reference_stations_are_restricted_to_modeled_span():
    target = optimizer.load_fit_target("wall_temperature", "RTE_model")
    case = SimpleNamespace(
        half_pass_nodes={
            "wall": np.array([0.25, 1.05]),
            "heat_flux": np.array([0.25, 1.05]),
            "coolant": np.array([0.25, 1.05]),
        },
        full_pass_nodes={
            "wall": np.array([1.05, -0.31]),
            "heat_flux": np.array([1.05, -0.31]),
            "coolant": np.array([1.05, -0.31]),
        },
    )

    restricted = optimizer.restrict_target_to_model_span(case, target)

    assert restricted.values.size == target.values.size - 1
    assert np.max(restricted.x) <= 1.05


def test_coolant_temperature_target_follows_half_then_full_pass():
    target = optimizer.load_fit_target("coolant_temperature")
    half_pass = SimpleNamespace(
        x_coolant=target.x[:2],
        T_static=np.array([10.0, 20.0]),
    )
    full_pass = SimpleNamespace(
        x_coolant=np.concatenate(([target.x[1]], target.x[2:])),
        T_static=np.array([20.0, 30.0, 40.0, 50.0, 60.0]),
    )

    predicted = optimizer.ChannelHeightObjective._predict_component(
        target,
        [half_pass, full_pass],
    )

    np.testing.assert_allclose(predicted, [10.0, 20.0, 30.0, 40.0, 50.0, 60.0])


def test_rte_coolant_temperature_uses_full_pass_only():
    target = optimizer.load_fit_target("coolant_temperature", "RTE_model")
    half_pass = SimpleNamespace(
        x_coolant=target.x,
        T_static=np.full(target.x.size, -1.0),
    )
    full_pass = SimpleNamespace(
        x_coolant=target.x,
        T_static=np.arange(target.x.size, dtype=float),
    )

    predicted = optimizer.ChannelHeightObjective._predict_component(
        target,
        [half_pass, full_pass],
    )

    np.testing.assert_allclose(predicted, np.arange(target.x.size, dtype=float))
