"""Assemble saved chamber results, an impeller, and an explicitly shaped pipe.

Run after examples/advanced/sizer_sim.py, or pass --input to saved results.
The arrangement is illustrative; it does not create fluid connections.
"""
import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyskyfire as psf
from pyskyfire.pump.impeller_new import Impeller


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
    # Route and radii are supplied explicitly; no routing or pressure loss is solved.
    pipe = SimpleNamespace(
        centerline=np.array([[0., 0., 0.], [0.04, 0., 0.], [0.08, 0.01, 0.], [0.12, 0.03, 0.]]),
        radius=0.003,
    )
    scene.add(psf.viz.PipeComponent(pipe), name="fuel_feed", anchor="outlet",
              attach_to=chamber_view.anchor("circuit_0.inlet"))
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
