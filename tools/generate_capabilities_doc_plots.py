"""Generate the interactive plots used by the capabilities documentation page."""

from pathlib import Path

import plotly.io as pio
import pyskyfire as psf


ROOT = Path(__file__).resolve().parents[1]
RESULTS_PATH = ROOT / "validation" / "RL10" / "regen_results.pkl"
DEMOS_DIR = ROOT / "docs" / "_static" / "demos"

ENGINE_3D_PATH = DEMOS_DIR / "engine_3d.html"
VISCOSITY_FIELD_PATH = DEMOS_DIR / "rl10_viscosity_field.html"


def write_plotly_demo(fig_or_plot, out_path: Path) -> None:
    """Write a Plotly figure as a standalone, responsive HTML demo."""
    fig = getattr(fig_or_plot, "fig", fig_or_plot)
    fig.update_layout(
        autosize=True,
        margin=dict(l=60, r=30, t=40, b=50),
    )

    pio.write_html(
        fig,
        file=str(out_path),
        full_html=True,
        include_plotlyjs="directory",
        default_width="100%",
        default_height="100%",
        auto_open=False,
        config={
            "responsive": True,
            "displaylogo": False,
        },
    )


def generate_engine_3d(thrust_chamber) -> None:
    """Generate the interactive 3D cooling-channel model for the RL10 case."""
    viewer = psf.viz.make_engine_3d(thrust_chamber, show=False)
    try:
        viewer.save_html(ENGINE_3D_PATH)
    finally:
        viewer.close()

    print(f"Wrote {ENGINE_3D_PATH.relative_to(ROOT)}")


def generate_viscosity_field(thrust_chamber) -> None:
    """Generate the hot-gas viscosity field for the RL10 case."""
    plot = psf.viz.PlotTransportPropertyField(
        thrust_chamber.combustion_transport,
        prop="mu",
        mode="surface",
    )
    write_plotly_demo(plot, VISCOSITY_FIELD_PATH)
    print(f"Wrote {VISCOSITY_FIELD_PATH.relative_to(ROOT)}")


def main() -> None:
    DEMOS_DIR.mkdir(parents=True, exist_ok=True)

    results = psf.common.Results.load(RESULTS_PATH)
    thrust_chamber = results.thrust_chamber

    generate_engine_3d(thrust_chamber)
    generate_viscosity_field(thrust_chamber)


if __name__ == "__main__":
    main()
