"""Mixture-ratio optimisation with Pyskyfire."""

import argparse
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from scipy.optimize import minimize_scalar
import pyskyfire as psf


def make_aerothermodynamics(MR: float) -> psf.skycea.Aerothermodynamics:
    fu = psf.common.Fluid(
        type="fuel",
        propellants=["CH4"],  # Ethanol for CEA
        fractions=[1.0],
    )

    ox = psf.common.Fluid(
        type="oxidizer",
        propellants=["O2"],  # Nitrous oxide for CEA
        fractions=[1.0],
    )

    return psf.skycea.Aerothermodynamics.from_F_eps_Lstar(
        fu=fu,
        ox=ox,
        MR=MR,
        p_c=100e5,      # Chamber pressure [Pa]
        F=100e3,         # Thrust [N]
        eps=60.0,      # Nozzle area ratio [-]
        L_star=1.2,    # Characteristic length [m]
    )


def vacuum_isp(MR: float) -> float:
    aero = make_aerothermodynamics(MR)
    return aero.Isp_vac


def objective(MR: float) -> float:
    return -vacuum_isp(MR)


def main(output_dir: Path | None = None) -> None:
    if output_dir is None:
        output_dir = Path(__file__).parent

    output_dir.mkdir(parents=True, exist_ok=True)

    # tutorial:start:optimisation
    result = minimize_scalar(
        objective,
        bounds=(1.5, 6.0),
        method="bounded",
        options={"xatol": 1e-3},
    )

    MR_opt = result.x
    aero_opt = make_aerothermodynamics(MR_opt)
    Isp_opt = aero_opt.Isp_vac
    # tutorial:end:optimisation

    print(f"Optimal mixture ratio: {MR_opt:.4f}")
    print(f"Vacuum specific impulse: {Isp_opt:.2f} s")
    print(f"Fuel mass flow: {aero_opt.mdot_fu:.5f} kg/s")
    print(f"Oxidizer mass flow: {aero_opt.mdot_ox:.5f} kg/s")
    print(f"Total mass flow: {aero_opt.mdot:.5f} kg/s")

    # tutorial:start:sweep
    MR_values = np.linspace(1.5, 6.0, 60)
    Isp_values = [vacuum_isp(MR) for MR in MR_values]
    # tutorial:end:sweep

    # tutorial:start:plot
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=MR_values,
            y=Isp_values,
            mode="lines",
            name="Vacuum specific impulse",
        )
    )

    fig.add_trace(
        go.Scatter(
            x=[MR_opt],
            y=[Isp_opt],
            mode="markers",
            name="Optimum",
            marker=dict(size=10),
        )
    )

    fig.update_layout(
        title="Mixture-ratio optimisation",
        xaxis_title="Mixture ratio, O/F [-]",
        yaxis_title="Vacuum specific impulse [s]",
    )

    fig.write_html(output_dir / "mixture-ratio-optimisation.html")
    # tutorial:end:plot

    # tutorial:start:contour
    xs, rs = psf.regen.contour.get_contour(
        V_c=aero_opt.V_c,
        AR_c=1.8,
        r_t=aero_opt.r_t,
        area_ratio=aero_opt.eps,
        nozzle="rao",
        R_1f=1,
        R_2f=2,
        R_3f=0.3,
    )

    contour = psf.regen.Contour(xs, rs, name="Optimised mixture-ratio contour")
    plot = psf.viz.PlotContour(contour)

    contour_fig = plot.fig
    contour_fig.write_html(output_dir / "optimized-contour.html")

    # tutorial:end:contour


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for generated HTML outputs.",
    )

    args = parser.parse_args()
    main(output_dir=args.output_dir)