"""Command-line preliminary impeller sizing, with explicit design inputs."""
import argparse

from pyskyfire.pump import Impeller


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Q", type=float, required=True, help="Delivered flow [m³/s]")
    parser.add_argument("--H", type=float, required=True, help="Head [m]")
    parser.add_argument("--n", type=float, required=True, help="Shaft speed [rpm]")
    parser.add_argument("--shaft-diameter", type=float, required=True, help="Shaft diameter [m]")
    parser.add_argument("--inlet-hub-diameter", type=float, required=True, help="Wetted inlet hub diameter [m]")
    parser.add_argument("--blade-count", type=int, required=True)
    parser.add_argument("--beta2-deg", type=float, required=True, help="Outlet metal angle [deg]")
    parser.add_argument("--blade-thickness", type=float, required=True, help="Total body thickness [m]")
    parser.add_argument("--phi1", type=float, default=0.20)
    parser.add_argument("--volumetric-efficiency", type=float, default=1.0)
    parser.add_argument("--hydraulic-efficiency", type=float)
    parser.add_argument("--mechanical-efficiency", type=float)
    parser.add_argument("--rho", type=float, help="Liquid density [kg/m³], for shaft power")
    impeller = Impeller(**vars(parser.parse_args()))
    print(impeller)
    if impeller.shaft_power is not None:
        print(f"Shaft power: {impeller.shaft_power:g} W")


if __name__ == "__main__":
    main()
