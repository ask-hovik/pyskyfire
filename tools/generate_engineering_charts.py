"""Generate the interactive charts used by the engineering explanations."""

from pathlib import Path

import pyskyfire as psf


OUTPUT_DIRECTORY = (
    Path(__file__).resolve().parents[1]
    / "docs"
    / "_static"
    / "explanation-artifacts"
)


def main() -> None:
    """Write each engineering chart to its own HTML file."""
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

    charts = {
        "rao-nozzle-angles.html": psf.viz.PlotThetaVsEpsilon(),
        "moody-diagram.html": psf.viz.PlotMoodyDiagram(),
    }

    for filename, chart in charts.items():
        output_path = OUTPUT_DIRECTORY / filename
        chart.save_html(output_path)
        print(f"Chart saved to {output_path}")


if __name__ == "__main__":
    main()
