"""Numerical regression guard for the RL10 regenerative-cooling case."""

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


RL10_SIMULATION = (
    Path(__file__).resolve().parents[2] / "validation" / "RL10" / "regen_sim.py"
)

# These are intentionally tight golden values, not experimental validation
# tolerances. A change here should prompt a review of why the numerical RL10
# solution moved before the new baseline is accepted.
RELATIVE_TOLERANCE = 1.0e-6
EXPECTED = {
    "peak_heat_flux": 30_902_258.786_860_1,
    "max_hot_wall_temperature": 1_028.656_784_141_646_8,
    "coolant_outlet_stagnation_temperature": 254.796_867_530_316_23,
    "coolant_outlet_stagnation_pressure": 5_903_235.081_605_2,
}


def _load_rl10_simulation():
    """Import the validation setup without requiring it to be a package."""
    module_name = "pyskyfire_test_rl10_regen_sim"
    spec = importlib.util.spec_from_file_location(module_name, RL10_SIMULATION)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load RL10 simulation from {RL10_SIMULATION}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def test_rl10_regenerative_solution_remains_numerically_stable() -> None:
    """Flag any one-part-per-million movement in the main RL10 outputs."""
    simulation = _load_rl10_simulation()
    _, cooling_passes = simulation.run_simulation(output=False)
    outlet = cooling_passes[-1]

    actual = {
        "peak_heat_flux": max(
            float(np.max(result.dQ_dA)) for result in cooling_passes
        ),
        "max_hot_wall_temperature": max(
            float(np.max(result.T[:, -1])) for result in cooling_passes
        ),
        # The validation report and the pass-to-pass handoff both define the
        # jacket outlet in terms of the stagnation state.
        "coolant_outlet_stagnation_temperature": float(
            outlet.T_stagnation[-1]
        ),
        "coolant_outlet_stagnation_pressure": float(outlet.p_stagnation[-1]),
    }

    assert actual == pytest.approx(EXPECTED, rel=RELATIVE_TOLERANCE, abs=0.0)
