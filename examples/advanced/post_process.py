"""
Generate an HTML report from the results written by sizer_sim.py.

The report always includes thrust-chamber and regenerative-cooling content.
It adds an Engine Cycle tab only when RUN_MODE = "full_cycle".
"""

from pathlib import Path

import pyskyfire as psf


RESULTS_FILENAME = "results.pkl"
REPORT_FILENAME = "methane_engine_report.html"


# tutorial:start:common-report
def add_common_report_content(report, params, thrust_chamber, cooling_data):
    """Add report tabs available for both regen-only and full-cycle results."""

    # Parameters
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

    # Engine overview
    tab_overview = report.add_tab("Engine Overview")
    engine_viewer = psf.viz.make_engine_3d(thrust_chamber, show=False)
    tab_overview.add_iframe(engine_viewer.data_url, caption="Engine 3D")
    engine_viewer.close()
    tab_overview.add_figure(psf.viz.PlotContour(thrust_chamber.contour))

    # Cooling data
    ordered_cooling_data = list(cooling_data.values())
    tab_cooling = report.add_tab("Cooling Data")
    tab_cooling.add_figure(
        psf.viz.PlotWallTemperature(
            *ordered_cooling_data,
            plot_hot=True,
            plot_coolant_wall=True,
            plot_interfaces=True,
        )
    )
    tab_cooling.add_figure(psf.viz.PlotCoolantTemperature(*ordered_cooling_data))
    tab_cooling.add_figure(psf.viz.PlotCoolantPressure(*ordered_cooling_data))
    tab_cooling.add_figure(psf.viz.PlotHeatFlux(*ordered_cooling_data))
    tab_cooling.add_figure(psf.viz.PlotVelocity(*ordered_cooling_data))

    # Thrust-chamber properties
    tab_chamber = report.add_tab("Thrust Chamber Properties")
    tab_chamber.add_figure(psf.viz.PlotCoolantArea(thrust_chamber))
    tab_chamber.add_figure(psf.viz.PlotHydraulicDiameter(thrust_chamber))
    tab_chamber.add_figure(
        psf.viz.PlotRadiusOfCurvature(thrust_chamber),
        caption="Radius-of-curvature computation is still experimental.",
    )
    tab_chamber.add_figure(psf.viz.PlotdAdxThermalHotGas(thrust_chamber))
    tab_chamber.add_figure(psf.viz.PlotdAdxThermalCoolant(thrust_chamber))
    tab_chamber.add_figure(psf.viz.PlotdAdxCoolantArea(thrust_chamber))

    # Combustion transport
    tab_combustion = report.add_tab("Combustion")
    transport = thrust_chamber.combustion_transport
    for prop in ("M", "gamma", "T", "p", "h", "cp", "k", "mu", "Pr", "rho", "a"):
        tab_combustion.add_figure(
            psf.viz.PlotTransportProperty(transport, prop=prop)
        )

    # Through-wall temperatures at three axial locations
    tab_gradient = report.add_tab("Thermal Gradient")
    fuel_throat_data = cooling_data["fuel_throat"]
    for x in (-0.1, 0.0, 0.05):
        tab_gradient.add_figure(
            psf.viz.PlotTemperatureProfile(
                fuel_throat_data,
                thrust_chamber,
                0,
                x,
            )
        )
# tutorial:end:common-report


# tutorial:start:cycle-report
def add_full_cycle_report_content(report, results):
    """Add network and station plots that require a full-cycle result."""

    stations = results["stations"]
    residuals = results["residuals"]

    tab_cycle = report.add_tab("Engine Cycle")
    tab_cycle.add_figure(
        psf.viz.PlotResidualHistory(residuals),
        caption="Maximum relative residual per fixed-point iteration.",
    )

    fuel_stations = [
        "fu_engine_in",
        "fu_pump_in",
        "fu_pump_out",
        "fu_regen_in",
        "fu_regen_out",
        "fu_turbine_in",
        "fu_turbine_out",
        "fu_injector_plenum",
        "fu_chamber_in",
    ]
    oxidizer_stations = [
        "ox_engine_in",
        "ox_pump_in",
        "ox_pump_out",
        "ox_regen_in",
        "ox_regen_1_in",
        "ox_regen_1_out",
        "ox_regen_out",
        "ox_turbine_in",
        "ox_turbine_out",
        "ox_injector_plenum",
        "ox_chamber_in",
    ]

    for property_name, caption in (
        ("p", "Fuel-side pressure"),
        ("T", "Fuel-side temperature"),
        ("mdot", "Fuel-side mass flow"),
    ):
        tab_cycle.add_figure(
            psf.viz.PlotStationProperty(
                station_dicts=stations,
                station_list=fuel_stations,
                property_name=property_name,
            ),
            caption=caption,
        )

    for property_name, caption in (
        ("p", "Oxidizer-side pressure"),
        ("T", "Oxidizer-side temperature"),
        ("mdot", "Oxidizer-side mass flow"),
    ):
        tab_cycle.add_figure(
            psf.viz.PlotStationProperty(
                station_dicts=stations,
                station_list=oxidizer_stations,
                property_name=property_name,
            ),
            caption=caption,
        )

    tab_cycle.add_figure(
        psf.viz.PlotPTDiagram(
            station_dicts=[stations],
            station_list=fuel_stations,
            fluid_name="methane",
            title="Fuel-side P-T path",
            scale="linear",
        )
    )
    tab_cycle.add_figure(
        psf.viz.PlotPTDiagram(
            station_dicts=[stations],
            station_list=oxidizer_stations,
            fluid_name="oxygen",
            title="Oxidizer-side P-T path",
            scale="linear",
        )
    )
# tutorial:end:cycle-report


# tutorial:start:generate-report
def main():
    script_dir = Path(__file__).resolve().parent
    input_path = script_dir / RESULTS_FILENAME
    output_path = script_dir / REPORT_FILENAME

    results = psf.common.Results.load(input_path)

    required_keys = {"mode", "params", "thrust_chamber", "cooling_data"}
    missing_keys = required_keys.difference(results)
    if missing_keys:
        missing = ", ".join(sorted(missing_keys))
        raise ValueError(
            f"{input_path.name} does not use the expected result format. "
            f"Missing: {missing}."
        )

    mode = results["mode"]
    if mode not in {"regen_only", "full_cycle"}:
        raise ValueError(f"Unknown result mode: {mode!r}")

    print(f"Generating {mode} report from {input_path}")

    report = psf.viz.Report("Methane Engine")
    add_common_report_content(
        report=report,
        params=results["params"],
        thrust_chamber=results["thrust_chamber"],
        cooling_data=results["cooling_data"],
    )

    if mode == "full_cycle":
        required_cycle_keys = {"stations", "residuals"}
        missing_cycle_keys = required_cycle_keys.difference(results)
        if missing_cycle_keys:
            missing = ", ".join(sorted(missing_cycle_keys))
            raise ValueError(
                f"Full-cycle result is missing required data: {missing}."
            )
        add_full_cycle_report_content(report, results)

    report.save_html(output_path)
    print(f"Report saved to {output_path}")
# tutorial:end:generate-report


if __name__ == "__main__":
    main()
