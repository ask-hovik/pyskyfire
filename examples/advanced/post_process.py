"""
Generate an HTML report from the results written by sizer_sim.py.

The report always includes thrust-chamber and regenerative-cooling content.
It adds an Engine Cycle tab only when RUN_MODE = "full_cycle".
"""

import argparse
from pathlib import Path

import pyskyfire as psf


RESULTS_FILENAME = "results.pkl"
REPORT_FILENAME = "methane_engine_report.html"


# tutorial:start:common-report
def add_common_report_content(output_dir, report, params, thrust_chamber, cooling_data):
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
    engine_viewer = psf.viz.make_engine_3d(
        thrust_chamber,
        stride=3,
        show=False,
    )
    engine_viewer.save_html(output_dir / "engine-3d.html")
    tab_overview.add_iframe(engine_viewer.data_url, caption="Engine 3D")
    engine_viewer.close()

    contour_plot = psf.viz.PlotContour(thrust_chamber.contour)
    tab_overview.add_figure(contour_plot)
    contour_plot.save_html(output_dir / "contour.html")

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
    tab_cooling.add_figure(
        psf.viz.PlotCoolantPressure(cooling_data["fuel_chamber"]),
        caption="Fuel coolant pressure.",
    )
    tab_cooling.add_figure(
        psf.viz.PlotCoolantPressure(
            cooling_data["oxidizer_nozzle_outbound"],
            cooling_data["oxidizer_nozzle_return"],
        ),
        caption="Oxidizer coolant pressure.",
    )
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
            psf.viz.PlotTransportProperty(
                transport, prop=prop, results=ordered_cooling_data
            )
        )

    # Through-wall temperatures at three axial locations
    tab_gradient = report.add_tab("Thermal Gradient")
    fuel_chamber_data = cooling_data["fuel_chamber"]
    profile_x = sorted(float(x) for x in fuel_chamber_data.x)
    for x in (profile_x[0], 0.0, profile_x[-1]):
        tab_gradient.add_figure(
            psf.viz.PlotTemperatureProfile(
                fuel_chamber_data,
                thrust_chamber,
                0,
                x,
            )
        )
# tutorial:end:common-report


# tutorial:start:cycle-report
def add_full_cycle_report_content(output_dir, report, results):
    """Add network and station plots that require a full-cycle result."""

    stations = results["stations"]
    residuals = results["residuals"]
    net = results["net"]

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
        "ox_stage1_pump_out",
        "ox_direct_branch",
        "ox_stage2_pump_in",
        "ox_regen_duct_in",
        "ox_regen_in",
        "ox_regen_out",
        "ox_turbine_in",
        "ox_turbine_out",
        "ox_main_flow_merge",
        "ox_injector_plenum",
        "ox_chamber_in",
    ]

    for property_name, caption in (
        ("p", "Fuel-side pressure"),
        ("T", "Fuel-side temperature"),
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
    ):
        tab_cycle.add_figure(
            psf.viz.PlotStationProperty(
                station_dicts=stations,
                station_list=oxidizer_stations,
                property_name=property_name,
            ),
            caption=caption,
        )

    sankey = psf.viz.PlotMassFlowSankey(engine_network=net, title="Methane Engine Mass Flow",)
    tab_cycle.add_figure(sankey)

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

    NETWORK_LAYOUT_PATH = (
    Path(__file__).resolve()
    .with_name("methane-engine-cycle.layout.json")
)
    tab_network = report.add_tab("Engine Network")
    viewer = psf.viz.make_network_viz(
        results["net"],
        title="Methane engine cycle",
        layout=NETWORK_LAYOUT_PATH,
    )
    viewer.save_html(path=output_dir / "network.html")
    tab_network.add_iframe(
        viewer.data_url,
        caption="Editable engine-cycle schematic",
        height="900px",
    )

# tutorial:end:cycle-report


# tutorial:start:generate-report
def main(
    output_dir: Path | None = None,
    input_path: Path | None = None,
):
    script_dir = Path(__file__).resolve().parent

    if output_dir is None:
        output_dir = script_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / REPORT_FILENAME

    if input_path is None:
        input_path = script_dir / RESULTS_FILENAME

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
        output_dir=output_dir,
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
        add_full_cycle_report_content(output_dir, report, results)

    report.save_html(output_path)
    print(f"Report saved to {output_path}")
# tutorial:end:generate-report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        dest="input_path",
        help="Results file generated by sizer_sim.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for generated HTML outputs.",
    )

    args = parser.parse_args()
    main(output_dir=args.output_dir, input_path=args.input_path)
