"""Minimal regenerative-cooling simulation for Pyskyfire."""

import argparse
from pathlib import Path

import pyskyfire as psf


def main(output_dir: Path | None = None) -> None:
    # tutorial:start:engine-inputs
    params = dict(
        p_c=20e5,                                                                               # Chamber Pressure (Pa)
        F=5e3,                                                                                  # Thrust (N)
        eps=10,                                                                                 # Nozzle area ratio
        L_star=1.2,                                                                             # Combustion chamber characteristic length 
        MR=2.8,                                                                                 # Mixture ratio
        AR_c=1.8,                                                                               # Chamber aspect ratio
        cea_fu=psf.common.Fluid(type="fuel", propellants=["C2H5OH"], fractions=[1.0]),          # Ethanol
        cea_ox=psf.common.Fluid(type="oxidizer", propellants=["N2O"], fractions=[1.0]),         # Nitrous oxide
        coolprop_fu=psf.common.Fluid(type="fuel", propellants=["ethanol"], fractions=[1.0]),    # Ethanol
        T_coolant_in=298.15,                                                                    # Coolant inlet temperature
        p_coolant_in=23e5,                                                                      # Coolant inlet pressure
        material=psf.common.solids.StainlessSteel304,                                           # Wall material
        wall_thickness=0.5e-3,                                                                  # Wall thickness
        n_channels=60,                                                                          # Number of cooling channels
        blockage_ratio=0.1,                                                                     # Fraction of cooling channel cross section filled with ribs
        roughness_height=10e-6,                                                                 # Cooling channel roughness parameter
        helix_angle=45,                                                                         # Cooling channel helix angle
        channel_height=2e-3                                                                     # Cooling channel height
    )
    # tutorial:end:engine-inputs

    # tutorial:start:aerothermodynamics
    # Create hot gas property object:
    aerothermodynamics = psf.skycea.Aerothermodynamics.from_F_eps_Lstar(
        fu=params["cea_fu"],        # Fuel input
        ox=params["cea_ox"],        # Ox input
        MR=params["MR"],            # Mixture ratio
        p_c=params["p_c"],          # Chamber pressure
        F=params["F"],              # Thrust
        eps=params["eps"],          # Area ratio
        L_star=params["L_star"],    # Characteristic length 
    )

    params["V_c"] = aerothermodynamics.V_c # Retrieve calculated chamber volume
    params["r_t"] = aerothermodynamics.r_t # Retrieve calculated throat radius
    params["mdot_fu"] = aerothermodynamics.mdot_fu # retrieve fuel mass flow rate

    # Create coolant property object:
    coolant_transport = psf.skycea.CoolantTransport(params["coolprop_fu"]) # 
    # tutorial:end:aerothermodynamics

    # tutorial:start:contour
    # Calculate contour points
    xs, rs = psf.regen.contour.get_contour(
        V_c=params["V_c"],          # Chamber volume
        AR_c=params["AR_c"],        # Chamber aspect ratio
        r_t=params["r_t"],          # Throat radius
        area_ratio=params["eps"],   # Area ratio
        nozzle="rao",               # Type of nozzle
        R_1f=1,                     # Throat radius parameter
        R_2f=2,                     # Chamber-to-contraction radius parameter
        R_3f=0.3,                   # Throat-to-nozzle radius parameter
    )

    # Create contour object
    contour = psf.regen.Contour(xs, rs, name="Minimal Contour")
    # tutorial:end:contour

    # tutorial:start:walls
    wall = psf.regen.Wall(
        material=params["material"],
        thickness=params["wall_thickness"],
    )
    # tutorial:end:walls

    # tutorial:start:cooling-circuit
    cross_section = psf.regen.CrossSectionSquared(blockage_ratio=params["blockage_ratio"])

    def channel_height_function(x):
        return params["channel_height"]

    def helix_angle_function(x):
        return params["helix_angle"] * 3.14 / 180

    placement = psf.regen.SurfacePlacement(
        n_channel_positions=params["n_channels"],
        helix_angle=helix_angle_function,
    )

    cooling_circuit = psf.regen.CoolingCircuit(
        name="Cooling Pass",
        contour=contour,
        coolant_transport=coolant_transport,
        cross_section=cross_section,
        span=[1.0, -1.0],
        placement=placement,
        walls=[wall],
        roughness=params["roughness_height"],
        channel_height=channel_height_function,
    )
    # tutorial:end:cooling-circuit

    # tutorial:start:thrust-chamber
    thrust_chamber = psf.regen.ThrustChamber(
        contour=contour,
        combustion_transport=aerothermodynamics,
        cooling_circuits=[cooling_circuit],
        n_nodes=150,
    )
    # tutorial:end:thrust-chamber

    # tutorial:start:simulation
    boundary_conditions = psf.regen.BoundaryConditions(
        T_coolant_in=params["T_coolant_in"],
        p_coolant_in=params["p_coolant_in"],
        mdot_coolant=params["mdot_fu"],
    )

    cooling_data = psf.regen.steady_heating_analysis(
        thrust_chamber,
        n_nodes=100,
        circuit_index=0,
        boundary_conditions=boundary_conditions,
        output=True,
    )
    # tutorial:end:simulation

    # tutorial:start:report
    if output_dir is None:
        output_dir = Path(__file__).parent

    output_dir.mkdir(parents=True, exist_ok=True)

    report = psf.viz.Report("Minimal Engine")

    tab_params = report.add_tab("Parameters")
    tab_params.add_table(
        params,
        caption="Input Parameters",
        key_title="Parameter",
        value_title="Value",
        precision=3,
    )
    tab_params.add_table(
        thrust_chamber.combustion_transport.optimum,
        caption="Optimal Values",
        key_title="Parameter",
        value_title="Value",
        precision=3,
    )

    tab_overview = report.add_tab("Engine Overview")

    # Engine contour
    contour_plot = psf.viz.PlotContour(thrust_chamber.contour)
    contour_plot.save_html(output_dir / "contour.html")
    tab_overview.add_figure(contour_plot)

    # 3D engine
    engine_viewer = psf.viz.make_engine_3d(thrust_chamber, show=False)
    engine_viewer.save_html(output_dir / "engine-3d.html")
    tab_overview.add_iframe(
        engine_viewer.data_url,
        caption="Engine 3D",
    )
    engine_viewer.close()

    tab_cooling_data = report.add_tab("Cooling Data")
    tab_cooling_data.add_figure(
        psf.viz.PlotWallTemperature(
            cooling_data,
            plot_hot=True,
            plot_coolant_wall=True,
        )
    )

    heat_flux = psf.viz.PlotHeatFlux(cooling_data)
    heat_flux.save_html(output_dir / "heat-flux.html")
    tab_cooling_data.add_figure(psf.viz.PlotCoolantTemperature(cooling_data))
    tab_cooling_data.add_figure(psf.viz.PlotCoolantPressure(cooling_data))
    tab_cooling_data.add_figure(heat_flux)
    tab_cooling_data.add_figure(psf.viz.PlotVelocity(cooling_data))

    tab_thrust_chamber = report.add_tab("Thrust Chamber Properties")
    tab_thrust_chamber.add_figure(psf.viz.PlotCoolantArea(thrust_chamber))
    tab_thrust_chamber.add_figure(psf.viz.PlotHydraulicDiameter(thrust_chamber))
    tab_thrust_chamber.add_figure(psf.viz.PlotRadiusOfCurvature(thrust_chamber))
    tab_thrust_chamber.add_figure(psf.viz.PlotdAdxThermalHotGas(thrust_chamber))
    tab_thrust_chamber.add_figure(psf.viz.PlotdAdxThermalCoolant(thrust_chamber))
    tab_thrust_chamber.add_figure(psf.viz.PlotdAdxCoolantArea(thrust_chamber))

    tab_combustion = report.add_tab("Combustion")
    for prop in ["M", "gamma", "T", "p", "h", "cp", "k", "mu", "Pr", "rho", "a"]:
        tab_combustion.add_figure(
            psf.viz.PlotTransportProperty(
                thrust_chamber.combustion_transport,
                prop=prop,
            )
        )

    tab_thermal_gradient = report.add_tab("Thermal Gradient")
    tab_thermal_gradient.add_figure(
        psf.viz.PlotTemperatureProfile(cooling_data, thrust_chamber, 0, -0.1)
    )
    tab_thermal_gradient.add_figure(
        psf.viz.PlotTemperatureProfile(cooling_data, thrust_chamber, 0, 0)
    )
    tab_thermal_gradient.add_figure(
        psf.viz.PlotTemperatureProfile(cooling_data, thrust_chamber, 0, 0.05)
    )

    report_path = output_dir / "minimal-report.html"
    report.save_html(report_path)
    print(f"Report saved to {report_path}")
    # tutorial:end:report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for generated HTML outputs.",
    )

    args = parser.parse_args()
    main(output_dir=args.output_dir)
