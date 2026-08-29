"""Compare independently selectable cooling-channel curvature corrections.

The example runs the same helical cooling circuit three times:

1. no curvature corrections;
2. RTE/Ito pressure-loss correction only;
3. RL10 coolant heat-transfer correction only.

All cases start from identical inlet conditions. The resulting coolant
temperature, pressure, heat-transfer coefficient, and hot-wall temperature
are overlaid in one HTML figure.
"""

import argparse
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import pyskyfire as psf


CASES = (
    ("No curvature correction", False, False),
    ("Pressure correction only", False, True),
    ("Heat-transfer correction only", True, False),
)


def build_thrust_chamber():
    """Build one curvature-rich helical cooling-channel example."""
    fuel = psf.common.Fluid(
        type="fuel",
        propellants=["C2H5OH"],
        fractions=[1.0],
    )
    oxidizer = psf.common.Fluid(
        type="oxidizer",
        propellants=["N2O"],
        fractions=[1.0],
    )
    coolant = psf.common.Fluid(
        type="fuel",
        propellants=["ethanol"],
        fractions=[1.0],
    )

    aerothermodynamics = psf.skycea.Aerothermodynamics.from_F_eps_Lstar(
        fu=fuel,
        ox=oxidizer,
        MR=2.8,
        p_c=20.0e5,
        F=10.0e3,
        eps=6.0,
        L_star=1.2,
    )

    xs, rs = psf.regen.contour.get_contour(
        V_c=aerothermodynamics.V_c,
        AR_c=1.5,
        r_t=aerothermodynamics.r_t,
        area_ratio=6.0,
        nozzle="rao",
        R_1f=1.0,
        R_2f=2.0,
        R_3f=0.3,
    )
    contour = psf.regen.Contour(xs, rs, name="Curvature comparison contour")

    placement = psf.regen.SurfacePlacement(
        # A small number of relatively deep channels keeps the RTE curvature
        # group above its applicability threshold in part of the passage, so
        # the pressure-only trace is visibly distinct.
        n_channel_positions=10,
        helix_angle=np.deg2rad(45.0),
    )
    circuit = psf.regen.CoolingCircuit(
        name="Helical cooling pass",
        contour=contour,
        coolant_transport=psf.skycea.CoolantTransport(coolant),
        cross_section=psf.regen.CrossSectionSquared(blockage_ratio=0.1),
        span=[1.0, -1.0],
        placement=placement,
        walls=[
            psf.regen.Wall(
                material=psf.common.solids.StainlessSteel304,
                thickness=0.5e-3,
            )
        ],
        roughness=10.0e-6,
        channel_height=lambda x: 4.0e-3,
    )

    chamber = psf.regen.ThrustChamber(
        contour=contour,
        combustion_transport=aerothermodynamics,
        cooling_circuits=[circuit],
        n_nodes=150,
    )
    boundary_conditions = psf.regen.BoundaryConditions(
        T_coolant_in=273.15,
        p_coolant_in=23.0e5,
        mdot_coolant=aerothermodynamics.mdot_fu,
    )
    return chamber, boundary_conditions


def run_cases(thrust_chamber, boundary_conditions):
    """Run all cases with independent public curvature switches."""
    results = {}
    for label, heat_correction, pressure_correction in CASES:
        print(f"\nRunning: {label}")
        results[label] = psf.regen.coupled_steady_heating_analysis(
            thrust_chamber,
            boundary_conditions=boundary_conditions,
            nodes=80,
            circuit_index=0,
            output=True,
            heat_curvature_correction=heat_correction,
            pressure_curvature_correction=pressure_correction,
        )
    return results


def make_comparison_figure(results):
    """Overlay all three cases in a single four-panel Plotly figure."""
    figure = make_subplots(
        rows=2,
        cols=2,
        shared_xaxes=True,
        subplot_titles=(
            "Bulk coolant temperature",
            "Coolant stagnation pressure",
            "Coolant-side heat-transfer coefficient",
            "Hot-wall temperature",
        ),
    )
    colors = ("#444444", "#0072B2", "#D55E00")

    for color, (label, result) in zip(colors, results.items()):
        common = dict(
            x=result.x,
            mode="lines",
            name=label,
            legendgroup=label,
            line=dict(color=color, width=2.5),
        )
        figure.add_trace(
            go.Scatter(y=result.T_static, showlegend=True, **common),
            row=1,
            col=1,
        )
        figure.add_trace(
            go.Scatter(
                y=np.asarray(result.p_stagnation) / 1.0e5,
                showlegend=False,
                **common,
            ),
            row=1,
            col=2,
        )
        figure.add_trace(
            go.Scatter(
                y=np.asarray(result.h_cold) / 1.0e3,
                showlegend=False,
                **common,
            ),
            row=2,
            col=1,
        )
        figure.add_trace(
            go.Scatter(
                y=np.asarray(result.T)[:, -1],
                showlegend=False,
                **common,
            ),
            row=2,
            col=2,
        )

    figure.update_xaxes(title_text="Axial position (m)", row=2, col=1)
    figure.update_xaxes(title_text="Axial position (m)", row=2, col=2)
    figure.update_yaxes(title_text="Temperature (K)", row=1, col=1)
    figure.update_yaxes(title_text="Pressure (bar)", row=1, col=2)
    figure.update_yaxes(title_text="h_cold (kW m^-2 K^-1)", row=2, col=1)
    figure.update_yaxes(title_text="Temperature (K)", row=2, col=2)
    figure.update_layout(
        title="Cooling-channel curvature-correction comparison",
        template="plotly_white",
        height=800,
        legend=dict(orientation="h", yanchor="bottom", y=1.08),
        margin=dict(l=70, r=30, t=130, b=60),
    )
    return figure


def print_summary(results) -> None:
    print("\nCase summary")
    for label, result in results.items():
        pressure_drop = result.p_stagnation[0] - result.p_stagnation[-1]
        print(
            f"{label:31s} "
            f"T_out={result.T_static[-1]:8.2f} K, "
            f"pressure drop={pressure_drop / 1.0e3:8.2f} kPa, "
            f"max T_hot_wall={np.max(result.T[:, -1]):8.2f} K"
        )


def main(output_dir: Path | None = None) -> None:
    if output_dir is None:
        output_dir = Path(__file__).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    thrust_chamber, boundary_conditions = build_thrust_chamber()
    results = run_cases(thrust_chamber, boundary_conditions)
    print_summary(results)

    output_path = output_dir / "curvature-correction-comparison.html"
    make_comparison_figure(results).write_html(output_path)
    print(f"Comparison plot saved to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for the generated comparison plot.",
    )
    arguments = parser.parse_args()
    main(output_dir=arguments.output_dir)
