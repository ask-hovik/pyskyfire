"""Export cooling-channel CAD input from a saved PySkyFire simulation.

This is intentionally an external example rather than a PySkyFire export API:
the shared library exposes centerlines, layouts, and lazy section planes, while
this script owns one particular Onshape JSON contract.

Squared sections use transverse construction planes and projected rib faces.
Rounded sections are transformed into transverse elliptical arcs. Optional
hot-wall profiles use the same planes. Simulation geometry remains unchanged.

Hot-wall profiles also carry ``runs``: spans of consecutive sections whose
curves are the same channels in the same order. The FeatureScript lofts those
curve by curve, so the CAD surface never has to guess which curve of one
profile belongs with which curve of the next.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pyskyfire as psf
from pyskyfire.regen.cad_geometry import build_hot_wall, transverse_construction


SUPPORTED_SECTIONS = (
    psf.regen.CrossSectionSquared,
    psf.regen.CrossSectionRounded,
    psf.regen.CrossSectionRoundedInternal,
)

RESULTS_FILENAME = "results.pkl"
OUTPUT_FILENAME = "onshape-cooling-channels.json"
METERS_TO_MILLIMETERS = 1000.0
RADIANS_TO_DEGREES = 180.0 / np.pi


def _onshape_vector(vector) -> list[float]:
    """Map PySkyFire [axial, y, z] into Onshape [x, y, axial]."""

    axial, radial_y, radial_z = np.asarray(vector, dtype=float)
    return [float(radial_y), float(radial_z), float(axial)]


def select_circuit_indices(thrust_chamber, selectors=None) -> list[int]:
    """Resolve zero-based indices or exact circuit names, preserving order."""

    circuits = thrust_chamber.cooling_circuits
    if not selectors:
        return list(range(len(circuits)))

    names = {circuit.name: index for index, circuit in enumerate(circuits)}
    selected = []
    for selector in selectors:
        try:
            index = int(selector)
        except ValueError:
            if selector not in names:
                available = ", ".join(
                    f"{index}: {circuit.name}" for index, circuit in enumerate(circuits)
                )
                raise ValueError(
                    f"Unknown cooling circuit {selector!r}. Available: {available}"
                ) from None
            index = names[selector]
        if not 0 <= index < len(circuits):
            raise ValueError(
                f"Cooling-circuit index {index} is outside 0..{len(circuits) - 1}"
            )
        if index not in selected:
            selected.append(index)
    return selected


def _section_data(construction) -> list[dict]:
    """Serialize the same transverse arcs sampled by the native viewer."""
    transverse = transverse_construction(construction)
    planes = transverse.planes
    if len(planes) < 2:
        raise ValueError("An Onshape channel loft requires at least two sections")
    sections = []
    for index in range(len(planes)):
        section = {
            "name": f"section_{index}",
            "innerArc": transverse.inner_arcs[index].to_data(),
            "outerArc": transverse.outer_arcs[index].to_data(),
            "circular": bool(construction.circular[index]),
            "plane": {
                "origin_mm": _onshape_vector(
                    planes.origins[index] * METERS_TO_MILLIMETERS
                ),
                "normal": _onshape_vector(planes.normals[index]),
                "xAxis": _onshape_vector(planes.x_axes[index]),
            },
        }
        sections.append(section)
    return sections


def build_onshape_data(
    thrust_chamber,
    *,
    layout: str = "analytical",
    circuit_selectors=None,
    include_hot_wall=False,
    hot_wall_trim_fraction=0.1,
    hot_wall_transition_fraction=0.02,
    hot_wall_transition_length_m=None,
) -> dict:
    """Return the JSON-serializable Onshape cooling-channel document."""

    selected_indices = select_circuit_indices(
        thrust_chamber,
        circuit_selectors,
    )
    if layout == "analytical":
        resolved = psf.regen.resolve_analytical_layout(thrust_chamber)
    elif layout == "numerical":
        resolved = psf.regen.resolve_numerical_layout(thrust_chamber)
    else:
        raise ValueError("layout must be 'analytical' or 'numerical'")
    resolved = resolved.for_circuits(selected_indices)

    channels = []
    for channel in resolved.channels:
        circuit = thrust_chamber.cooling_circuits[channel.circuit_index]
        # The FeatureScript only knows how to sketch the two section families
        # it has cases for, so refuse anything else here rather than exporting
        # arcs it would silently mis-sketch.
        if not isinstance(circuit.cross_section, SUPPORTED_SECTIONS):
            raise TypeError(
                "Unsupported channel cross section for the Onshape export: "
                f"{type(circuit.cross_section).__name__}"
            )
        construction = psf.regen.build_channel_construction(
            circuit, channel.reference_path
        )
        channels.append(
            {
                "name": channel.name,
                "circuitName": circuit.name,
                "circuitIndex": channel.circuit_index,
                "sectionType": construction.section_type,
                "pattern": {
                    "count": channel.pattern_count,
                    "angularStep_deg": float(
                        channel.pattern_angle_rad * RADIANS_TO_DEGREES
                    ),
                },
                "sections": _section_data(construction),
            }
        )

    data = {
        "schema": "pyskyfire.onshape.cooling-channels",
        "schemaVersion": 6,
        "units": {"length": "millimeter", "angle": "degree"},
        "coordinateSystem": {
            "description": (
                "PySkyFire axial X maps to Onshape Z; the mapping is right-handed"
            ),
            "patternAxis": {
                "origin_mm": [0.0, 0.0, 0.0],
                "direction": [0.0, 0.0, 1.0],
            },
        },
        "layout": {
            "method": resolved.method,
            "selectedCircuitIndices": selected_indices,
        },
        "channels": channels,
    }

    if include_hot_wall:
        data["hotWall"] = build_hot_wall(
            thrust_chamber,
            resolved,
            trim_fraction=hot_wall_trim_fraction,
            transition_fraction=hot_wall_transition_fraction,
            transition_length_m=hot_wall_transition_length_m,
        )
    return data


def export_onshape_json(
    input_path: Path,
    output_path: Path,
    *,
    layout: str = "analytical",
    circuit_selectors=None,
    include_hot_wall=False,
    hot_wall_trim_fraction=0.1,
    hot_wall_transition_fraction=0.02,
    hot_wall_transition_length_m=None,
) -> Path:
    """Load a saved simulation and write its Onshape JSON document."""

    results = psf.common.Results.load(input_path)
    if "thrust_chamber" not in results:
        raise ValueError(f"{input_path} does not contain a thrust_chamber result")
    data = build_onshape_data(
        results["thrust_chamber"],
        layout=layout,
        circuit_selectors=circuit_selectors,
        include_hot_wall=include_hot_wall,
        hot_wall_trim_fraction=hot_wall_trim_fraction,
        hot_wall_transition_fraction=hot_wall_transition_fraction,
        hot_wall_transition_length_m=hot_wall_transition_length_m,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(data, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return output_path


def parse_args(
    *,
    script_dir=None,
    results_filename=RESULTS_FILENAME,
    output_filename=OUTPUT_FILENAME,
):
    parser = argparse.ArgumentParser(
        description="Export PySkyFire cooling channels for the included Onshape FeatureScript."
    )
    script_dir = (
        Path(script_dir) if script_dir is not None else Path(__file__).resolve().parent
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=script_dir / results_filename,
        help=f"Saved PySkyFire results pickle (default: {results_filename}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir,
        help="Directory for the generated JSON file.",
    )
    parser.add_argument(
        "--output-name",
        default=output_filename,
        help=f"Generated JSON filename (default: {output_filename}).",
    )
    parser.add_argument(
        "--layout",
        choices=("analytical", "numerical"),
        default="analytical",
        help="Layout resolver to use (default: analytical).",
    )
    parser.add_argument(
        "--circuit",
        action="append",
        dest="circuits",
        metavar="INDEX_OR_NAME",
        help=(
            "Export one circuit by zero-based index or exact name. Repeat for a "
            "subset; omit to export all circuits."
        ),
    )
    parser.add_argument(
        "--hot-wall", action="store_true", help="Include continuous hot-wall profiles."
    )
    parser.add_argument(
        "--hot-wall-trim-fraction",
        type=float,
        default=0.1,
        help="Fraction of each arc parameter span trimmed per blended side (default: 0.1).",
    )
    transition = parser.add_mutually_exclusive_group()
    transition.add_argument(
        "--hot-wall-transition-fraction",
        type=float,
        default=0.02,
        help="Transition length / full contour axial length (default: 0.02).",
    )
    transition.add_argument(
        "--hot-wall-transition-length-m",
        type=float,
        help="Absolute transition length in metres instead of a fraction.",
    )
    return parser.parse_args()


def main(**defaults) -> None:
    args = parse_args(**defaults)
    output_path = export_onshape_json(
        args.input,
        args.output_dir / args.output_name,
        layout=args.layout,
        circuit_selectors=args.circuits,
        include_hot_wall=args.hot_wall,
        hot_wall_trim_fraction=args.hot_wall_trim_fraction,
        hot_wall_transition_fraction=args.hot_wall_transition_fraction,
        hot_wall_transition_length_m=args.hot_wall_transition_length_m,
    )
    print(f"Onshape cooling-channel data saved to {output_path}")


if __name__ == "__main__":
    main()
