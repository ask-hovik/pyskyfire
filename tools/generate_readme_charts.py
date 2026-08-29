"""Export the README charts from the RL10 validation case.

GitHub strips ``<script>`` and ``<iframe>`` from rendered markdown, so a Plotly
figure cannot stay interactive in the README. This writes two artifacts from the
same figure object:

* a PNG for the README itself, and
* a standalone interactive HTML published with the docs, which the README image
  links to.

The figure is built by ``validation/RL10/post_process.py`` so the README chart
and the validation report cannot diverge.
"""

import argparse
import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import pyskyfire as psf  # noqa: E402

RL10_DIR = REPO_ROOT / "validation" / "RL10"
REGEN_RESULTS = RL10_DIR / "regen_results.pkl"

PNG_OUTPUT = REPO_ROOT / "images" / "RL10_wall_temperature.png"
HTML_OUTPUT = (
    REPO_ROOT / "docs" / "_static" / "readme-artifacts" / "rl10-wall-temperature.html"
)


def load_post_process():
    """Import the RL10 post-processor without requiring it to be a package."""

    path = RL10_DIR / "post_process.py"
    spec = importlib.util.spec_from_file_location("rl10_post_process", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--height", type=int, default=620)
    parser.add_argument("--scale", type=float, default=2.0,
                        help="Pixel density multiplier for the PNG.")
    args = parser.parse_args()

    if not REGEN_RESULTS.is_file():
        raise FileNotFoundError(
            f"Missing {REGEN_RESULTS}. Run validation/RL10/regen_sim.py first."
        )

    post_process = load_post_process()
    results = psf.common.Results.load(str(REGEN_RESULTS))
    cooling_data_a, cooling_data_b = results.cooling_data[0], results.cooling_data[1]

    figure = post_process.build_wall_temperature_comparison(
        cooling_data_a, cooling_data_b
    )

    HTML_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.save_html(HTML_OUTPUT)
    print(f"Wrote {HTML_OUTPUT}")

    PNG_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.write_image(
        str(PNG_OUTPUT),
        width=args.width,
        height=args.height,
        scale=args.scale,
    )
    size_kb = PNG_OUTPUT.stat().st_size / 1024
    print(f"Wrote {PNG_OUTPUT} ({size_kb:.0f} KB, "
          f"{int(args.width * args.scale)}x{int(args.height * args.scale)})")


if __name__ == "__main__":
    main()
