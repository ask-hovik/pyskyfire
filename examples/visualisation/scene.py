"""Assemble saved chamber results, an impeller, and a solved feed pipe.

Run after examples/advanced/sizer_sim.py, or pass --input to saved results.
The arrangement is illustrative; it does not create fluid connections.
"""
import argparse
from pathlib import Path

import numpy as np
import pyskyfire as psf
from pyskyfire.pump.impeller_new import Impeller


def build_feed_pipe():
    """A short bent feed line, resolved and pressure-loss marched.

    The route is stated as control vertices; `generate_pipe_stations` turns
    each corner into its two tangent stations and one lumped bend loss.
    """
    nodes = [
        psf.ducts.PathNode((0.0, 0.0, 0.0), diameter=0.012),
        psf.ducts.PathNode((0.09, 0.0, 0.0), diameter=0.012,
                           bend=psf.ducts.RoundedBend(radius=0.030)),
        psf.ducts.PathNode((0.16, 0.06, 0.0), diameter=0.012),
    ]
    stations, connections = psf.ducts.generate_pipe_stations(nodes)
    return psf.ducts.Pipe(stations, connections, roughness=15e-6)


def solve_feed_pipe(pipe):
    """Colour the duct with its own solution rather than invented numbers."""
    fluid = psf.thermo.Fluid(type="fuel", propellants=["Methanol"], fractions=[1.0])
    transport = psf.thermo.FluidTransport(fluid)
    resolver = psf.thermo.FluidStateResolver(transport, p_ref=5.0e6)
    return pipe.solve(
        inlet_state=resolver.from_tp(300.0, 5.0e6),
        mass_flow=1.5,
        transport=transport,
    )


def build_scene(thrust_chamber, cooling_data):
    scene = psf.viz.Scene(color_by="p_static")
    chamber_view = scene.add(
        psf.viz.ThrustChamberComponent(thrust_chamber, cooling_data=cooling_data,
                                      layout="analytical", stride=3),
        name="chamber",
    )
    impeller = Impeller(Q=0.01, H=100, n=10000)
    scene.add(
        psf.viz.ImpellerComponent(impeller), name="fuel_pump",
        transform=psf.viz.Transform(translation=(-0.3, 0.15, 0), rotation=(0, np.pi/2, 0)),
    )
    # The duct carries its own geometry and solution; the scene only places it.
    pipe = build_feed_pipe()
    scene.add(psf.viz.PipeComponent(pipe, solve_feed_pipe(pipe)), name="fuel_feed",
              anchor="outlet", attach_to=chamber_view.anchor("circuit_0.inlet"))
    return scene


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(__file__).parents[1]/"advanced"/"results.pkl")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    results = psf.common.Results.load(args.input)
    scene = build_scene(results["thrust_chamber"], results["cooling_data"])
    scene.save_html(args.output_dir/"engine.html")
    if args.show:
        scene.show()


if __name__ == "__main__":
    main()
