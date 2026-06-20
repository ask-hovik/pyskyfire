"""Minimal regenerative-cooling simulation for Pyskyfire."""

from pathlib import Path
import pyskyfire as psf


def main() -> None:
    # tutorial:start:engine-inputs
    params = dict(
        p_c=50e5,
        F=5e3,
        eps=10,
        L_star=1.2,
        MR=2.8,
        AR_c=1.8,
        cea_fu=psf.common.Fluid(
            type="fuel", propellants=["C2H5OH"], fractions=[1.0]
        ),
        cea_ox=psf.common.Fluid(
            type="oxidizer", propellants=["N2O"], fractions=[1.0]
        ),
        coolprop_fu=psf.common.Fluid(
            type="fuel", propellants=["ethanol"], fractions=[1.0]
        ),
        T_coolant_in=298.15,
        p_coolant_in=23e5,
        material=psf.common.solids.StainlessSteel304,
        wall_thickness=0.5e-3,
        n_channels=60,
        blockage_ratio=0.0,
        roughness_height=10e-6,
    )
    # tutorial:end:engine-inputs

    # tutorial:start:aerothermodynamics
    aerothermodynamics = psf.skycea.Aerothermodynamics.from_F_eps_Lstar(
        fu=params["cea_fu"],
        ox=params["cea_ox"],
        MR=params["MR"],
        p_c=params["p_c"],
        F=params["F"],
        eps=params["eps"],
        L_star=params["L_star"],
    )

    params["V_c"] = aerothermodynamics.V_c
    params["r_t"] = aerothermodynamics.r_t

    coolant_transport = psf.skycea.CoolantTransport(params["coolprop_fu"])
    # tutorial:end:aerothermodynamics

    # tutorial:start:contour-and-wall
    xs, rs = psf.regen.contour.get_contour(
        V_c=params["V_c"],
        AR_c=params["AR_c"],
        r_t=params["r_t"],
        area_ratio=params["eps"],
        nozzle="rao",
        R_1f=1,
        R_2f=2,
        R_3f=0.3,
    )

    contour = psf.regen.Contour(xs, rs, name="Minimal Contour")
    wall = psf.regen.Wall(
        material=params["material"],
        thickness=params["wall_thickness"],
    )
    # tutorial:end:contour-and-wall

    # tutorial:start:cooling-circuit
    cross_section = psf.regen.CrossSectionSquared(
        blockage_ratio=params["blockage_ratio"]
    )

    def channel_height_function(x):
        return 2e-3

    def helix_angle_function(x):
        return 45 * 3.14 / 180

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

    # tutorial:start:visualisation
    plot_3d, _ = psf.viz.make_engine_3d(thrust_chamber)
    plot_3d.show()
    del plot_3d
    # tutorial:end:visualisation>

    # tutorial:start:simulation
    boundary_conditions = psf.regen.BoundaryConditions(
        T_coolant_in=params["T_coolant_in"],
        p_coolant_in=params["p_coolant_in"],
        mdot_coolant=aerothermodynamics.mdot_fu,
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
    out_path = Path(__file__).with_name("minimal_report.html")
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
    tab_overview.add_figure(psf.viz.PlotContour(thrust_chamber.contour))

    tab_cooling_data = report.add_tab("Cooling Data")
    tab_cooling_data.add_figure(
        psf.viz.PlotWallTemperature(
            cooling_data,
            plot_hot=True,
            plot_coolant_wall=True,
        )
    )
    tab_cooling_data.add_figure(psf.viz.PlotCoolantTemperature(cooling_data))
    tab_cooling_data.add_figure(psf.viz.PlotCoolantPressure(cooling_data))
    tab_cooling_data.add_figure(psf.viz.PlotHeatFlux(cooling_data))
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

    report.save_html(out_path)
    print(f"Report saved to {out_path}")
    # tutorial:end:report


if __name__ == "__main__":
    main()
