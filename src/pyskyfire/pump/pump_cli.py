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
    parser.add_argument("--blade-thickness", type=float, required=True, help="Total body thickness [m]")
    parser.add_argument("--hydraulic-efficiency", type=float, required=True,
                        help="The impeller's own hydraulic efficiency")
    size = parser.add_mutually_exclusive_group()
    size.add_argument("--beta2-deg", type=float, help="Outlet metal angle [deg]; d2 is then solved")
    size.add_argument("--psi", type=float, help="Head coefficient; beta2 is then solved")
    size.add_argument("--d2", type=float, help="Outlet diameter [m]; beta2 is then solved")
    size.add_argument("--u2", type=float, help="Tip speed [m/s]; beta2 is then solved")
    parser.add_argument("--fT", type=float, default=1.0, help="Eq. (3.26) band factor, when no size is given")
    parser.add_argument("--phi1", type=float, default=0.20)
    parser.add_argument("--volumetric-efficiency", type=float, default=1.0)
    parser.add_argument("--rho", type=float, help="Liquid density [kg/m³], for blade power")
    impeller = Impeller(**vars(parser.parse_args()))
    print(impeller)
    checks = impeller.diagnostics
    print(f"Solved {checks.solved}: beta2 = {impeller.velocity_triangles.outlet[1].beta_metal_deg:.4g} deg, "
          f"psi = {checks.psi:.4g} (Eq. 3.26: {checks.psi_correlation:.4g})")
    print(f"Slip factor {checks.slip.slip_factor:.4g}, cu2 = {checks.slip.cu:.4g} m/s, "
          f"deviation {checks.deviation_deg:.3g} deg, w2/w1 = {checks.deceleration_ratio:.3g}")
    if impeller.blade_power is not None:
        print(f"Blade power: {impeller.blade_power:g} W")


if __name__ == "__main__":
    main()
