"""Counterflow hydrogen/oxygen exchanger; both fluids use CoolProp.

Run with: python examples/heat_exchanger/two_fluid.py --output-dir /tmp/hx
"""

import argparse
from pathlib import Path
import numpy as np
import pyskyfire as psf


def build_exchanger(nodes=24):
    hx = psf.heat_exchanger
    length = 0.3
    wall = psf.regen.Wall(
        psf.common.Material("Example metal", k=psf.common.ConstantModel(80.0)), 0.0005
    )
    section = psf.regen.CrossSectionSquared(
        0.0005, floor="line", roof="line", exposure={"floor": 1.0}
    )
    down = psf.ducts.Path.from_nodes(
        [psf.ducts.PathNode((0, 0, 0)), psf.ducts.PathNode((length, 0, 0))]
    ).refined(length / nodes)
    up = psf.ducts.Path.from_nodes(
        [psf.ducts.PathNode((0, 0, 0)), psf.ducts.PathNode((-length, 0, 0))]
    ).refined(length / nodes)

    def bank(name, path, fluid, start, normal):
        return hx.HeatExchangerCircuit(
            name=name,
            path=path,
            cross_section=section,
            wall=wall,
            transport=psf.thermo.FluidTransport(
                psf.thermo.Fluid(type="fuel", propellants=[fluid], fractions=[1.0])
            ),
            dimensions=dict(height=0.003, floor_width=0.005, roof_width=0.005),
            placements=[
                psf.ducts.PathPlacement(
                    (start[0], start[1], i * 0.006), initial_floor_direction=normal
                )
                for i in range(12)
            ],
            roughness=2e-6,
        )

    cold = bank("oxygen", down, "Oxygen", (0, 0.002), (0, -1, 0))
    hot = bank("hydrogen", up, "Hydrogen", (length, -0.002), (0, 1, 0))
    n = len(down.s) - 1
    exchanger = hx.HeatExchanger(
        [cold, hot],
        thermal_links=[
            hx.ThermalLink(
                "oxygen",
                "hydrogen",
                tuple((i, n - 1 - i) for i in range(n)),
                psf.regen.Wall(wall.material, 0.001),
            )
        ],
    )
    inlets = {
        "oxygen": hx.Inlet(130.0, 8e6, 0.25),
        "hydrogen": hx.Inlet(300.0, 8e6, 0.05),
    }
    return exchanger, inlets


def main(output_dir=None, nodes=24):
    exchanger, inlets = build_exchanger(nodes)
    result = exchanger.solve(inlets)
    output_dir = Path(output_dir or Path(__file__).parent / "output")
    output_dir.mkdir(parents=True, exist_ok=True)
    psf.viz.heat_exchanger_report(
        exchanger, result, title="Hydrogen / oxygen counterflow exchanger"
    ).save_html(output_dir / "two-fluid-report.html")
    for name, r in result.circuits.items():
        print(f"{name}: {r.T[0]:.2f} -> {r.T[-1]:.2f} K; Q={r.total_heat:.1f} W")
    print(
        f"Converged in {result.iterations} iterations; energy residual {result.energy_balance_error:.6g} W"
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--nodes", type=int, default=24)
    args = parser.parse_args()
    main(args.output_dir, args.nodes)
