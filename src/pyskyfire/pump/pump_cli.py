# pump_cli.py
import argparse
import sys

# If "centrifugal_pump" is a package (with __init__.py):
#  - Option A: import the Impeller class directly if you've exposed it in __init__.py
#       from centrifugal_pump import Impeller
#
#  - Option B: if not exposed, import from the submodule:
#       from centrifugal_pump.centrifugal_pump import Impeller
#
# The exact import depends on how your package is structured.
from centrifugal_pump import Impeller
from centrifugal_pump import plot


def main():
    parser = argparse.ArgumentParser(
        description="Command-line interface to run centrifugal pump (impeller) calculations."
    )
    parser.add_argument("--Q", type=float, required=True,
                        help="Flow rate in m^3/s (e.g., 0.05)")
    parser.add_argument("--H", type=float, required=True,
                        help="Head in meters (e.g., 20)")
    parser.add_argument("--n", type=int, required=True,
                        help="Rotational speed in rpm (e.g., 3500)")
    
    # Optionally, you could add more arguments to control plotting, geometry, etc.
    # parser.add_argument("--plot", action="store_true",
    #                     help="Show 3D plot of the impeller (if applicable).")
    
    args = parser.parse_args()
    
    # Create an Impeller object with the user-specified arguments
    imp = Impeller(Q=args.Q, H=args.H, n=args.n)
    print(imp)
    plot.plot_impeller_views(imp)
    plot.plot_impeller_3D(imp)
    
    # Print some results
    # The Impeller class in your code places results in imp.results
    # which is a list. We can loop over them and print.
    #print("=== Pump Calculation Results ===")
    #for idx, result_dict in enumerate(imp.results):
    #    print(f"Result set {idx+1}:")
    #    for key, val in result_dict.items():
    #        print(f"  {key}: {val}")
    #print("=================================")


if __name__ == "__main__":
    # Entry point for the script
    main()
