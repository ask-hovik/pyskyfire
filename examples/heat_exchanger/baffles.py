"""Regeneratively cooled paired ring and radial-spoke baffles.

The two rings lean toward one another downstream, making a revolved V.
Each spoke has a straight downflow bank and a back-to-back return bank.
Series headers are ideal connections; their end-turn geometry/loss is omitted.
The dimensions and constant metal conductivity target an easy numerical
example, not a qualified engine design.
"""

import argparse
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import pyskyfire as psf


def build_baffles(nodes=24):
    hx = psf.heat_exchanger
    fuel = psf.thermo.Fluid(type="fuel", propellants=["H2"], fractions=[1.0])
    oxidizer = psf.thermo.Fluid(type="oxidizer", propellants=["O2"], fractions=[1.0])
    gas = psf.thermo.Aerothermodynamics.from_F_eps_Lstar(
        fu=fuel, ox=oxidizer, MR=5.5, p_c=3e6, F=100e3, eps=8, L_star=1.2
    )
    xs, rs = psf.regen.contour.generate_contour(
        V_c=gas.V_c,
        AR_c=1.0,
        r_t=gas.r_t,
        area_ratio=8,
        nozzle="rao",
        R_1f=1,
        R_2f=2,
        R_3f=0.3,
    )
    contour = psf.regen.Contour(xs, rs, name="Prescribed hydrogen / oxygen chamber")
    coolant = psf.thermo.FluidTransport(
        psf.thermo.Fluid(type="fuel", propellants=["Hydrogen"], fractions=[1.0])
    )
    metal = psf.common.Material(
        "Demonstration copper", k=psf.common.ConstantModel(250.0)
    )
    chamber_cooling = psf.regen.CoolingCircuit(
        name="Chamber cooling",
        contour=contour,
        cross_section=psf.regen.CrossSectionSquared(sidewall_half_thickness=0.0005),
        span=[-1.0, 1.0],
        placement=psf.regen.SurfacePlacement(n_channel_positions=96, helix_angle=0.0),
        channel_height=0.003,
        walls=[psf.regen.Wall(metal, thickness=0.001)],
        coolant_transport=coolant,
        roughness=2e-6,
    )
    chamber = psf.regen.ThrustChamber(
        contour=contour,
        cooling_circuits=[chamber_cooling],
        combustion_transport=gas,
        n_nodes=100,
    )
    thickness = 0.0004
    wall = psf.regen.Wall(metal, thickness)
    shared_spoke_wall = psf.regen.Wall(metal, 2 * thickness)
    height = 0.003
    ring_count = 48
    length = min(0.12, 0.30 * abs(float(xs[0])))
    x0 = 0.002  # injector datum; +x points downstream
    inner_radius = 0.45 * contour.r_c
    outer_radius = inner_radius + 0.008
    lean = 0.002
    angles = np.linspace(0, 2 * np.pi, ring_count, endpoint=False)
    arc_section = psf.regen.CrossSectionSquared(
        thickness,
        floor="arc",
        roof="arc",
        exposure={"floor": 1.0, "roof": 1.0},
        enable_fin=True,
    )
    # Identical flow areas in both banks keep the ideal header endpoints equal.
    target_area = (
        height * (2 * np.pi * inner_radius / ring_count - 2 * thickness) * 0.94
    )

    def ring(name, start_x, start_r, dx, dr):
        path = psf.ducts.Path.from_nodes(
            [psf.ducts.PathNode((0, 0, 0)), psf.ducts.PathNode((dx, dr, 0))]
        ).refined(np.hypot(dx, dr) / nodes)
        total = path.s[-1]

        def radius(s):
            return start_r + dr * s / total - height / 2

        def angle(s):
            def residual(theta):
                prof = psf.regen.SectionProfiles.from_dimensions(
                    height=height,
                    wall_thickness=thickness,
                    radius=radius(s),
                    angle=theta,
                )
                return float(arc_section.A_coolant(prof)[0]) - target_area

            return brentq(
                residual, 2.01 * thickness / radius(s), 2 * np.pi / ring_count
            )

        placements = [
            psf.ducts.PathPlacement(
                (start_x, start_r * np.cos(phi), start_r * np.sin(phi)),
                rotation=psf.ducts.rotation_about_axis((1, 0, 0), phi),
                initial_floor_direction=(0, -np.cos(phi), -np.sin(phi)),
            )
            for phi in angles
        ]
        return hx.HeatExchangerCircuit(
            name=name,
            path=path,
            cross_section=arc_section,
            wall=wall,
            transport=coolant,
            dimensions=dict(height=height, radius=radius, angle=angle),
            placements=placements,
            roughness=2e-6,
        )

    inner = ring("inner ring down", x0, inner_radius, length, lean)
    outer = ring("outer ring return", x0 + length, outer_radius - lean, -length, lean)
    # Six radial blades, five channels per blade, with a return layer behind it.
    spoke_count = 6
    radial_count = 5
    rmin = outer_radius + 0.006
    rmax = 0.9 * contour.r_c
    pitch = (rmax - rmin) / radial_count
    width = pitch - 2 * thickness
    rect_section = psf.regen.CrossSectionSquared(
        thickness,
        floor="line",
        roof="line",
        exposure={"floor": 1.0, "roof": 1.0},
        enable_fin=True,
    )

    def spokes(name, returning):
        sign = -1 if returning else 1
        start_x = x0 + length if returning else x0
        path = psf.ducts.Path.from_nodes(
            [psf.ducts.PathNode((0, 0, 0)), psf.ducts.PathNode((sign * length, 0, 0))]
        ).refined(length / nodes)
        placements = []
        for phi in np.linspace(0, 2 * np.pi, spoke_count, endpoint=False):
            radial = np.array([0, np.cos(phi), np.sin(phi)])
            tangent = np.array([0, -np.sin(phi), np.cos(phi)])
            for r in rmin + (np.arange(radial_count) + 0.5) * pitch:
                start = (
                    np.array([start_x, 0, 0])
                    + r * radial
                    + sign * (height / 2 + thickness) * tangent
                )
                placements.append(
                    psf.ducts.PathPlacement(
                        tuple(start), initial_floor_direction=tuple(sign * tangent)
                    )
                )
        return hx.HeatExchangerCircuit(
            name=name,
            path=path,
            cross_section=rect_section,
            wall=wall,
            transport=coolant,
            dimensions=dict(height=height, floor_width=width, roof_width=width),
            placements=placements,
            roughness=2e-6,
        )

    down = spokes("spokes down", False)
    up = spokes("spokes return", True)
    n = len(down.s) - 1
    exchanger = hx.HeatExchanger(
        [inner, outer, down, up],
        flow_connections=[
            hx.FlowConnection(inner.name, outer.name),
            hx.FlowConnection(down.name, up.name),
        ],
        thermal_links=[
            hx.ThermalLink(
                down.name,
                up.name,
                tuple((i, n - 1 - i) for i in range(n)),
                shared_spoke_wall,
                surface_a="roof",
                surface_b="roof",
            )
        ],
    )
    inlets = {
        inner.name: hx.Inlet(150.0, 6e6, 0.25),
        down.name: hx.Inlet(150.0, 6e6, 0.12),
    }
    # Prescribed exterior coefficient deliberately makes the demonstration
    # reproducible. Gas driving temperatures still come from the CEA solution.
    installation = chamber.install_heat_exchanger(
        exchanger, position=(float(xs[0]), 0.0, 0.0), coefficient=1000.0
    )
    return chamber, installation, inlets


def build_scene(chamber, installation, result=None):
    """Use the simulation placement for the baffles in the chamber scene."""
    scene = psf.viz.Scene(
        title="Cooled chamber and injector baffles",
        color_by="T_coolant_wall" if result is not None else None,
        show_axes=True,
    )
    chamber_node = scene.add(
        psf.viz.ThrustChamberComponent(chamber, layout="analytical", stride=2),
        name="chamber",
    )
    for group in chamber_node.children:
        for channel in group.children:
            channel.opacity = 0.16
    scene.add(
        psf.viz.HeatExchangerComponent(installation, result, samples=4), name="baffles"
    )
    return scene


def main(output_dir=None, nodes=24):
    chamber, installation, inlets = build_baffles(nodes)
    result = installation.solve(inlets)
    output_dir = Path(output_dir or Path(__file__).parent / "output")
    output_dir.mkdir(parents=True, exist_ok=True)
    scene = build_scene(chamber, installation, result)
    scene.save_html(output_dir / "baffle-scene.html")
    report = psf.viz.heat_exchanger_report(
        installation.exchanger,
        result,
        title="Cooled ring and spoke baffles",
        installation=installation,
        scene=scene,
    )
    report.add_tab("Assumptions").add_markdown(
        "The ring pair forms a shallow revolved V. Each radial spoke has two straight banks joined by an ideal series header. "
        "All sections are squared: arced ring floors/roofs and rectangular spokes. "
        "The scene includes 96 straight chamber cooling channels from injector to nozzle exit, "
        "shown transparently. Their geometry is constructed without a separate chamber-cooling solve. "
        "Baffles start 2 mm downstream of the injector in both simulation and scene. "
        "A fixed exterior coefficient of 1000 W/m²/K is used with CEA recovery temperature. "
        "Symmetric rib enhancement is enabled. Metal conductivity is fixed at 250 W/m/K. "
        "Header geometry, header losses, radiation, and gas-flow feedback are omitted. "
        "These choices demonstrate convergence, not a qualified wall-temperature design."
    )
    report.save_html(output_dir / "baffle-report.html")
    for name, r in result.circuits.items():
        walls = [w.temperature_a for w in result.walls if w.circuit_a == name] + [
            w.temperature_b
            for w in result.walls
            if w.circuit_b == name or (w.circuit_a == name and w.circuit_b is None)
        ]
        print(
            f"{name}: {r.T[0]:.2f} -> {r.T[-1]:.2f} K; max wall {max(walls):.2f} K; Q={r.total_heat:.1f} W"
        )
    print(
        f"Converged in {result.iterations} iterations; energy residual {result.energy_balance_error:.6g} W"
    )
    print(
        f"Chamber volume displaced: {installation.volume_change.chamber_displaced:.6g} m³"
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--nodes", type=int, default=24)
    args = parser.parse_args()
    main(args.output_dir, args.nodes)
