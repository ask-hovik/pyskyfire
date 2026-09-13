"""Small, geometry-only engine for debugging the Onshape hot-wall export.

Run from the repository root:
    python examples/hot_wall_export/export_hot_wall.py

No combustion calculation, cooling solve, saved simulation or viewer is needed.
"""

import argparse
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path

import numpy as np
import pyskyfire as psf


EXAMPLE_DIR = Path(__file__).resolve().parent
EXPORTER_PATH = EXAMPLE_DIR.parent / "advanced" / "export_onshape.py"


# Edit these settings, then run this file. Command-line options are optional
# overrides; these values are the single source of defaults for every entry point.
params = dict(
    section="squared",  # "squared", "rounded", or "mixed"
    n_channels=50,  # Total count; mixed splits it equally
    n_stations=7,  # Axial loft profiles per circuit
    helix_angle_deg=15.0,  # Helix angle [degrees]
    channel_height=2e-3,  # Height [m]; see rounded_height_input
    rounded_height_input="base",  # "base" excludes lobes; "full" includes them
    rib_width=0.3e-3,  # Squared channel rib width [m]
    wall_thickness=0.3e-3,  # Main wall thickness [m]
    material=psf.common.solids.StainlessSteel304,
    circuit_span=[-1.0, 1.0],  # Signed contour fractions; 0 = throat
    contour_xs=np.array([-90, -60, -30, 0, 30, 60, 90]) * 1e-3,  # Axial positions [m]
    contour_rs=np.array([30, 30, 24, 18, 22, 29, 36]) * 1e-3,  # Radii [m] contour_rs=np.full(7, 30e-3),#
    chamber_angle_deg=0.0,  # Chamber-side slope angle at throat
    nozzle_angle_deg=0.0,  # Nozzle-side slope angle at throat
    trim_fraction=0.3,  # Trim per blended side of hot-wall arc
    transition_fraction=0.02,  # Fraction of full contour axial length
    transition_length_m=None,  # Optional length [m], overrides fraction
    output_dir=EXAMPLE_DIR,
    output_name=None,  # None -> hot-wall-{section}.json
)


def build_engine(**overrides):
    """Build geometry using ``params``, optionally overriding individual settings."""
    settings = params | overrides
    section = settings["section"]
    n_channels = settings["n_channels"]
    n_stations = settings["n_stations"]
    helix_angle_deg = settings["helix_angle_deg"]
    if section not in ("squared", "rounded", "mixed"):
        raise ValueError("section must be squared, rounded, or mixed")
    if n_channels < 4 or (section == "mixed" and n_channels % 2):
        raise ValueError(
            "Use at least four channels, with an even total for mixed sections"
        )
    if n_stations < 3:
        raise ValueError("Use at least three axial stations")
    if not np.isfinite(helix_angle_deg) or abs(helix_angle_deg) >= 90:
        raise ValueError(
            "Helix angle must be finite and strictly between -90 and 90 degrees"
        )

    # Hand-specified, gently varying contour: 180 mm long, throat at X=0.
    # Seven default loft stations coincide with these contour points.
    contour = psf.regen.Contour(
        xs=settings["contour_xs"],
        rs=settings["contour_rs"],
        name="Hot-wall export debug engine",
        chamber_angle=np.radians(settings["chamber_angle_deg"]),
        nozzle_angle=np.radians(settings["nozzle_angle_deg"]),
    )
    families = ("squared", "rounded") if section == "mixed" else (section,)
    circuits = []
    for family in families:
        cross_section = (
            psf.regen.CrossSectionSquared(rib_width=settings["rib_width"])
            if family == "squared"
            else psf.regen.CrossSectionRounded(
                height_input=settings["rounded_height_input"]
            )
        )
        # Base-height input keeps the rounded construction valid with very few
        # channels. Its full height includes the lobes and is consequently larger.
        circuits.append(
            psf.regen.CoolingCircuit(
                name=f"{family.capitalize()} channels",
                contour=contour,
                cross_section=cross_section,
                span=settings["circuit_span"],
                placement=psf.regen.SurfacePlacement(
                    n_channel_positions=n_channels // len(families),
                    helix_angle=np.radians(helix_angle_deg),
                ),
                channel_height=settings["channel_height"],
                walls=[
                    psf.regen.Wall(
                        settings["material"], thickness=settings["wall_thickness"]
                    )
                ],
                coolant_transport=None,
                roughness=0.0,
            )
        )
    return psf.regen.ThrustChamber(
        contour=contour,
        cooling_circuits=circuits,
        combustion_transport=None,
        compute_gas=False,
        n_nodes=n_stations,
    )


def export_example(**overrides):
    """Export ``params``; keyword arguments override just the supplied settings."""
    settings = params | overrides
    section = settings["section"]
    n_channels = settings["n_channels"]
    n_stations = settings["n_stations"]
    output_name = settings["output_name"]
    # Reuse the production example exporter rather than duplicate its geometry
    # or JSON schema. Load by absolute path so this script works from any cwd.
    spec = spec_from_file_location("hot_wall_example_exporter", EXPORTER_PATH)
    exporter = module_from_spec(spec)
    spec.loader.exec_module(exporter)
    engine = build_engine(**settings)
    data = exporter.build_onshape_data(
        engine,
        include_hot_wall=True,
        hot_wall_trim_fraction=settings["trim_fraction"],
        hot_wall_transition_fraction=settings["transition_fraction"],
        hot_wall_transition_length_m=settings["transition_length_m"],
    )
    output_dir = Path(settings["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / (output_name or f"hot-wall-{section}.json")
    path.write_text(
        json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    wall = data["hotWall"]["sections"]
    print(f"Saved {path} ({path.stat().st_size / 1024:.1f} KiB)")
    print(
        f"{n_channels} channels; {len(data['channels'])} coolant loft seed(s); "
        f"{n_stations} profiles per seed"
    )
    print(
        f"Hot wall: {len(wall)} profiles, "
        f"{sum(len(s['curves']) for s in wall)} total sketch curves"
    )
    print(
        "Import the JSON with examples/advanced/onshape_cooling_channels.fs, "
        "then enable Construct hot wall."
    )
    return path


def main():
    # Suppress absent arguments so they cannot override edits to params.
    parser = argparse.ArgumentParser(
        description=__doc__, argument_default=argparse.SUPPRESS
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--output-name", help="Override the generated JSON filename.")
    parser.add_argument("--section", choices=("squared", "rounded", "mixed"))
    parser.add_argument(
        "--channels", dest="n_channels", type=int, help="Override params['n_channels']."
    )
    parser.add_argument(
        "--stations", dest="n_stations", type=int, help="Override params['n_stations']."
    )
    parser.add_argument("--helix-angle-deg", type=float)
    parser.add_argument("--trim-fraction", type=float)
    export_example(**vars(parser.parse_args()))


if __name__ == "__main__":
    main()
