#!/usr/bin/env python3
"""Diagnose CEA property kinks in the RL10 regenerative-cooling solve.

The diagnostic is deliberately separate from the production solver.  It:

1. samples the temperature-conditioned CEA state and the Bartz reference
   state at the long-tube inlet;
2. writes the sampled properties to CSV and plots the 200 K clamp and the
   273.15 K ice/liquid transition;
3. repeats the two-pass solve with the baseline CEA settings, CEA's native
   smooth-species-truncation option, and a conductivity-only smoothing
   intervention around Pyskyfire's 200 K CEA clamp.

The last comparison is a causal test: all CEA quantities except the
conductivity supplied to the Bartz correlation remain unchanged.

Run from the repository root, for example::

    .venv/bin/python tools/diagnose_rl10_cea_kink.py \
        --output-dir /tmp/rl10-cea-diagnostic
"""

from __future__ import annotations

import argparse
import contextlib
import csv
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import sys
import tempfile
import warnings

import cea

_matplotlib_config = Path(tempfile.gettempdir()) / "pyskyfire-matplotlib"
_matplotlib_config.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_matplotlib_config))
import matplotlib.pyplot as plt
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

import pyskyfire as psf  # noqa: E402
from pyskyfire.regen.coupled_solver import (  # noqa: E402
    CoupledHeatExchangerPhysics,
)
from pyskyfire.skycea.aerothermodynamics import (  # noqa: E402
    MIN_CEA_TEMPERATURE,
)


DEFAULT_RESULTS = (
    REPOSITORY_ROOT / "validation" / "RL10_2" / "regen_results.pkl"
)
ICE_LIQUID_TRANSITION_K = 273.15


@dataclass(frozen=True)
class SolverCase:
    name: str
    inlet_x_m: float
    inlet_hot_wall_temperature_K: float
    inlet_scaled_residual: float
    inlet_converged: bool
    failed_wall_nodes: int
    maximum_scaled_residual: float
    failed_x_m: list[float]


def _reference_state(transport, x: float, hot_wall_temperature: float):
    """Return the enthalpy-driven Bartz reference state and wall state."""
    gas = transport.get_state(x)
    wall = transport.get_state(x, T=hot_wall_temperature)
    dynamic_enthalpy = 0.5 * gas.M**2 * gas.a**2
    reference_enthalpy = (
        0.5 * (wall.h + gas.h) + 0.18 * dynamic_enthalpy
    )
    return wall, transport.get_state(x, h=reference_enthalpy), reference_enthalpy


def sample_properties(chamber, x: float, temperatures: np.ndarray):
    """Sample the CEA and Bartz-reference properties at one axial station."""
    transport = chamber.combustion_transport
    transport.reset_solver_state()
    rows = []
    for hot_wall_temperature in temperatures:
        wall, reference, reference_enthalpy = _reference_state(
            transport, x, float(hot_wall_temperature)
        )
        rows.append(
            {
                "requested_hot_wall_temperature_K": float(hot_wall_temperature),
                "cea_wall_temperature_K": float(wall.T),
                "wall_enthalpy_J_per_kg": float(wall.h),
                "wall_conductivity_W_per_m_K": float(wall.k),
                "wall_cp_J_per_kg_K": float(wall.cp),
                "reference_enthalpy_J_per_kg": float(reference_enthalpy),
                "reference_temperature_K": float(reference.T),
                "reference_conductivity_W_per_m_K": float(reference.k),
                "reference_cp_J_per_kg_K": float(reference.cp),
                "reference_viscosity_Pa_s": float(reference.mu),
                "reference_Prandtl": float(reference.Pr),
                "solid_water_mole_fraction": float(
                    wall.X.get("H2O(cr)", 0.0)
                ),
                "liquid_water_mole_fraction": float(
                    wall.X.get("H2O(L)", 0.0)
                ),
            }
        )
    return rows


def _write_property_csv(rows, path: Path):
    with path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def _plot_properties(rows, path: Path):
    values = {key: np.asarray([row[key] for row in rows], dtype=float)
              for key in rows[0]}
    temperature = values["requested_hot_wall_temperature_K"]
    figure, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

    axes[0].plot(
        temperature,
        values["wall_enthalpy_J_per_kg"] / 1.0e6,
        color="tab:blue",
    )
    axes[0].set_ylabel("CEA wall enthalpy (MJ/kg)")

    axes[1].plot(
        temperature,
        values["wall_conductivity_W_per_m_K"],
        label="Temperature-conditioned CEA state",
    )
    axes[1].plot(
        temperature,
        values["reference_conductivity_W_per_m_K"],
        label="Bartz reference state",
    )
    axes[1].set_ylabel("Conductivity (W/m/K)")
    axes[1].legend()

    axes[2].plot(
        temperature,
        values["solid_water_mole_fraction"],
        label="H2O(cr)",
    )
    axes[2].plot(
        temperature,
        values["liquid_water_mole_fraction"],
        label="H2O(L)",
    )
    axes[2].set_ylabel("Mole fraction")
    axes[2].set_xlabel("Trial hot-wall temperature (K)")
    axes[2].legend()

    for axis in axes:
        axis.axvline(
            MIN_CEA_TEMPERATURE,
            color="tab:orange",
            linestyle="--",
            label="Pyskyfire CEA clamp",
        )
        axis.axvline(
            ICE_LIQUID_TRANSITION_K,
            color="tab:red",
            linestyle=":",
            label="H2O(cr)/H2O(L)",
        )
        axis.grid(alpha=0.25)

    figure.suptitle("RL10 long-tube inlet: CEA property transitions")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def _configure_native_smooth_truncation(transport, width: float):
    """Rebuild both native CEA solvers with smooth species truncation."""
    reactants, products = transport._make_mixtures()
    options = {
        "reactants": reactants,
        "transport": True,
        "smooth_truncation": True,
        "truncation_width": float(width),
    }
    transport._eq_solver = cea.EqSolver(products, **options)
    transport._eq_solution = cea.EqSolution(transport._eq_solver)
    transport._rocket_solver = cea.RocketSolver(products, **options)
    transport._rocket_solution = cea.RocketSolution(
        transport._rocket_solver
    )
    transport.clear_cache()


def _smoothstep(value: float) -> float:
    return value * value * (3.0 - 2.0 * value)


@contextlib.contextmanager
def _smooth_bartz_conductivity_at_clamp(
    transport,
    *,
    half_width_K: float,
):
    """Smooth only Bartz's conductivity across the 200 K clamp.

    The original gas-side calculation is evaluated first.  Inside the blend
    interval its heat-transfer quantities are scaled by ``k_smooth/k_raw``;
    the Bartz correlation is linear in conductivity.  Enthalpy, cp, viscosity,
    Prandtl number, wall conduction, and coolant properties are untouched.
    """
    original = CoupledHeatExchangerPhysics._gas_side_coefficients
    endpoint_cache: dict[float, tuple[float, float]] = {}
    lower = MIN_CEA_TEMPERATURE - half_width_K
    upper = MIN_CEA_TEMPERATURE + half_width_K

    def reference_conductivity(x: float, hot_wall_temperature: float) -> float:
        _, reference, _ = _reference_state(
            transport, x, hot_wall_temperature
        )
        return float(reference.k)

    def smoothed(self, x, T_hw, T_gr_mode="reference"):
        coefficients = original(self, x, T_hw, T_gr_mode=T_gr_mode)
        if not lower < T_hw < upper:
            return coefficients

        cache_key = round(float(x), 12)
        if cache_key not in endpoint_cache:
            endpoint_cache[cache_key] = (
                reference_conductivity(x, lower),
                reference_conductivity(x, upper),
            )
        k_lower, k_upper = endpoint_cache[cache_key]
        blend = _smoothstep((float(T_hw) - lower) / (upper - lower))
        k_smooth = k_lower + blend * (k_upper - k_lower)
        k_raw = reference_conductivity(x, float(T_hw))
        scale = k_smooth / k_raw
        for key in ("h_hot", "h_g", "h_gr", "qpp_hot"):
            coefficients[key] *= scale
        return coefficients

    CoupledHeatExchangerPhysics._gas_side_coefficients = smoothed
    try:
        yield
    finally:
        CoupledHeatExchangerPhysics._gas_side_coefficients = original


def _run_two_pass_simulation(
    results_path: Path,
    case_name: str,
    *,
    native_smooth_width: float | None = None,
    conductivity_blend_width: float | None = None,
) -> SolverCase:
    """Repeat both passes so the production CEA call history is preserved."""
    saved = psf.common.Results.load(results_path)
    chamber = saved.thrust_chamber
    transport = chamber.combustion_transport
    transport.reset_solver_state()
    if native_smooth_width is not None:
        _configure_native_smooth_truncation(
            transport, native_smooth_width
        )

    inlet_boundary = psf.regen.BoundaryConditions(
        T_coolant_in=float(saved.params["T_coolant_in"]),
        p_coolant_in=float(saved.params["p_coolant_in"]),
        mdot_coolant=float(transport.mdot_fu),
    )
    half_node_count = len(saved.cooling_data[0].x_wall)
    full_node_count = len(saved.cooling_data[1].x_wall)
    smoothing = (
        _smooth_bartz_conductivity_at_clamp(
            transport, half_width_K=conductivity_blend_width
        )
        if conductivity_blend_width is not None
        else contextlib.nullcontext()
    )
    with smoothing, warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        half_result = psf.regen.coupled_steady_heating_analysis(
            chamber,
            nodes=half_node_count,
            circuit_index=0,
            boundary_conditions=inlet_boundary,
            solver="newton",
            output=False,
            heat_curvature_correction=True,
        )
        full_boundary = psf.regen.BoundaryConditions(
            T_coolant_in=float(half_result.T_stagnation[-1]),
            p_coolant_in=float(half_result.p_stagnation[-1]),
            mdot_coolant=float(transport.mdot_fu),
        )
        result = psf.regen.coupled_steady_heating_analysis(
            chamber,
            nodes=full_node_count,
            circuit_index=1,
            boundary_conditions=full_boundary,
            solver="newton",
            output=False,
            heat_curvature_correction=True,
        )

    failed = np.flatnonzero(~np.asarray(result.wall_converged, dtype=bool))
    return SolverCase(
        name=case_name,
        inlet_x_m=float(result.x_wall[0]),
        inlet_hot_wall_temperature_K=float(result.T[0, -1]),
        inlet_scaled_residual=float(result.wall_residual_scaled[0]),
        inlet_converged=bool(result.wall_converged[0]),
        failed_wall_nodes=int(failed.size),
        maximum_scaled_residual=float(
            np.nanmax(result.wall_residual_scaled)
        ),
        failed_x_m=[float(result.x_wall[index]) for index in failed],
    )


def run_diagnostic(
    results_path: Path,
    output_dir: Path,
    *,
    sample_min_K: float = 180.0,
    sample_max_K: float = 320.0,
    sample_step_K: float = 0.05,
    conductivity_blend_half_width_K: float = 10.0,
    native_truncation_width: float = 0.25,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    saved = psf.common.Results.load(results_path)
    inlet_x = float(saved.cooling_data[1].x_wall[0])
    temperatures = np.arange(
        sample_min_K,
        sample_max_K + 0.5 * sample_step_K,
        sample_step_K,
    )
    rows = sample_properties(saved.thrust_chamber, inlet_x, temperatures)
    _write_property_csv(rows, output_dir / "cea_property_sweep.csv")
    _plot_properties(rows, output_dir / "cea_property_kinks.png")

    cases = [
        _run_two_pass_simulation(results_path, "baseline"),
        _run_two_pass_simulation(
            results_path,
            "cea_smooth_species_truncation",
            native_smooth_width=native_truncation_width,
        ),
        _run_two_pass_simulation(
            results_path,
            "conductivity_only_clamp_smoothing",
            conductivity_blend_width=conductivity_blend_half_width_K,
        ),
    ]
    summary = {
        "results_path": str(results_path.resolve()),
        "inlet_x_m": inlet_x,
        "pyskyfire_minimum_cea_temperature_K": MIN_CEA_TEMPERATURE,
        "ice_liquid_transition_K": ICE_LIQUID_TRANSITION_K,
        "cases": [asdict(case) for case in cases],
        "causal_checks": {
            "conductivity_only_smoothing_resolves_inlet": (
                not cases[0].inlet_converged and cases[2].inlet_converged
            ),
            "cea_smooth_species_truncation_resolves_inlet": (
                not cases[0].inlet_converged and cases[1].inlet_converged
            ),
        },
    }
    with (output_dir / "solver_comparison.json").open(
        "w", encoding="utf-8"
    ) as output_file:
        json.dump(summary, output_file, indent=2)
        output_file.write("\n")
    return summary


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "tools" / "rl10_cea_diagnostic",
    )
    parser.add_argument("--sample-step-K", type=float, default=0.05)
    parser.add_argument(
        "--conductivity-blend-half-width-K", type=float, default=10.0
    )
    parser.add_argument("--native-truncation-width", type=float, default=0.25)
    return parser.parse_args()


def main():
    args = _parse_args()
    summary = run_diagnostic(
        args.results,
        args.output_dir,
        sample_step_K=args.sample_step_K,
        conductivity_blend_half_width_K=(
            args.conductivity_blend_half_width_K
        ),
        native_truncation_width=args.native_truncation_width,
    )
    print(json.dumps(summary, indent=2))
    print(f"\nWrote diagnostics to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
