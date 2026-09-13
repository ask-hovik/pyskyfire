"""Generate the interactive charts used by the engineering explanations."""

import argparse
from pathlib import Path

import pyskyfire as psf


OUTPUT_DIRECTORY = (
    Path(__file__).resolve().parents[1]
    / "docs"
    / "_static"
    / "explanation-artifacts"
)


def main() -> None:
    """Write all engineering charts, or only the selected chart."""
    charts = {
        "rao": ("rao-nozzle-angles.html", psf.viz.PlotThetaVsEpsilon),
        "moody": ("moody-diagram.html", psf.viz.PlotMoodyDiagram),
    }
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--chart", choices=charts,
        help="Generate only this chart (default: generate all charts).",
    )
    args = parser.parse_args()
    selected = [args.chart] if args.chart else charts
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

    for name in selected:
        filename, chart_factory = charts[name]
        output_path = OUTPUT_DIRECTORY / filename
        chart_factory().save_html(output_path)
        print(f"Chart saved to {output_path}")


if __name__ == "__main__":
    main()
