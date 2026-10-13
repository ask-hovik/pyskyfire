"""Hold the coolant area constant by root-finding the helix angle."""

import argparse
from pathlib import Path

import numpy as np
import pyskyfire as psf


def main(output_dir: Path | None = None) -> None:
    if output_dir is None:
        output_dir = Path(__file__).parent

    output_dir.mkdir(parents=True, exist_ok=True)

    # tutorial:start:geometry
    xs, rs = psf.regen.contour.generate_contour(
        V_c=1.5e-3,         # Chamber volume [m³]
        AR_c=1.5,           # Chamber aspect ratio [-]
        r_t=0.025,          # Throat radius [m]
        area_ratio=8,       # Nozzle area ratio [-]
        nozzle="rao",
        R_1f=1, R_2f=2, R_3f=0.3,
    )
    contour = psf.regen.Contour(xs, rs, name="Contour")
    wall = psf.regen.Wall(material=psf.common.solids.StainlessSteel304, thickness=0.8e-3)

    def build_chamber(helix_angle):
        circuit = psf.regen.CoolingCircuit(
            name="Cooling Pass",
            contour=contour,
            coolant_transport=None,
            cross_section=psf.regen.CrossSectionSquared(sidewall_half_thickness=0.5e-3),
            span=[1.0, -1.0],
            placement=psf.regen.SurfacePlacement(n_channel_positions=40, helix_angle=helix_angle),
            walls=[wall],
            roughness=10e-6,
            channel_height=2e-3,
        )
        return psf.regen.ThrustChamber(
            contour=contour,
            combustion_transport=None,
            cooling_circuits=[circuit],
            n_nodes=400,
        )

    x_stations = np.asarray(build_chamber(0.0).cooling_circuits[0].x_domain)

    def tabulated_helix(cos_gamma):
        """Helix law interpolating cos(gamma) tabulated on the stations [rad]."""
        def helix_angle(x):
            return float(np.arccos(np.interp(x, x_stations, cos_gamma)))
        return helix_angle

    def coolant_area(cos_gamma):
        """Coolant area at every station for a tabulated helix [m²]."""
        chamber = build_chamber(tabulated_helix(cos_gamma))
        return np.asarray(chamber.cooling_circuits[0].A_coolant_vals)
    # tutorial:end:geometry

    # tutorial:start:root-find
    cos_a = np.ones_like(x_stations)             # axial channel
    A_a = coolant_area(cos_a)
    A_target = A_a.min()                         # narrowest axial area, at the throat

    cos_b = np.full_like(x_stations, np.cos(np.radians(30.0)))
    A_b = coolant_area(cos_b)

    while np.abs(A_b - A_target).max() > 1e-6 * A_target:
        # One secant step per station, all stations from the same build.
        # A station that has stopped moving takes no step.
        step = (A_b - A_target) * (cos_b - cos_a) / np.where(A_b != A_a, A_b - A_a, np.inf)
        cos_a, A_a = cos_b, A_b
        cos_b = np.clip(cos_b - step, 0.0, 1.0)
        A_b = coolant_area(cos_b)
        print(f"max area error: {np.abs(A_b - A_target).max() / A_target:.1e}")
    # tutorial:end:root-find

    # tutorial:start:scene
    thrust_chamber = build_chamber(tabulated_helix(cos_b))
    gamma = np.degrees(thrust_chamber.cooling_circuits[0].helix_angle_vals)
    print(f"target area {A_target * 1e6:.3f} mm², helix angle {gamma.min():.1f}° to {gamma.max():.1f}°")

    scene = psf.viz.Scene()
    scene.add(psf.viz.ThrustChamberComponent(thrust_chamber, stride=1), name="chamber")
    scene.save_html(output_dir / "engine-3d.html")
    # tutorial:end:scene


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for generated HTML outputs.",
    )

    args = parser.parse_args()
    main(output_dir=args.output_dir)
