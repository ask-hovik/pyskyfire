import os
import pyskyfire as psf
import matplotlib.pyplot as plt
import webview

script_dir = os.path.dirname(os.path.abspath(__file__))
in_path = os.path.join(script_dir, "methane_engine_sizer_res.pkl")

res = psf.common.Results.load(in_path)

net = res.net
input_params = res.input_params

stations = net.stations
signals = net.signals
residuals = net.residuals
thrust_chamber = input_params["thrust_chamber"]

# post process results: 
fu_cooling_data = net.block_results["fu_regen"]
ox_cooling_data_1 = net.block_results["ox_regen_1"]
ox_cooling_data_2 = net.block_results["ox_regen_2"]

psf.regen.plot.plot_contour([thrust_chamber.contour], colors=["black"], linestyles=["-"], show_labels=False, title=None)

psf.common.plot.plot_residual_history(residuals, color="black", linestyle="-", marker="^")

print(f"Isp: {thrust_chamber.optimal_values.Isp}")

#psf.regen.plot_engine_3D(thrust_chamber, show_axis=False)

# station list now available both from reference data and model
fu_stations = [
    "fu_engine_in", "fu_pump_in", "fu_pump_out", "fu_regen_in", "fu_regen_out", 
    "fu_turbine_in", "fu_turbine_housing_out", "fu_injector_plenum", 
    "fu_chamber_in"
]

ox_stations = [
    "ox_engine_in", "ox_pump_in", "ox_pump_out", "ox_regen_in", "ox_regen_out",  
    "ox_turbine_in", "ox_turbine_housing_out", "ox_injector_plenum", 
    "ox_chamber_in"
]

# fuel plots
fig, ax = psf.common.plot.plot_station_property(
    station_dicts=[stations],  
    station_list=fu_stations,
    property_name="p",
    ylabel="p  [Pa]", 
    color_cycle=["tab:red"],
    title= None
)

fig, ax = psf.common.plot.plot_station_property(
    station_dicts=[stations],  
    station_list=fu_stations,
    property_name="T",
    ylabel="T  [K]", 
    color_cycle=["tab:red"],
    title = None
)

fig, ax = psf.common.plot.plot_station_property(
    station_dicts=[stations],  
    station_list=fu_stations,
    property_name="mdot",
    ylabel="mdot  [kg/s]", 
    color_cycle=["tab:red"],
    title = None
)

# ox plots
fig, ax = psf.common.plot.plot_station_property(
    station_dicts=[stations],  
    station_list=ox_stations,
    property_name="p",
    ylabel="p  [Pa]", 
    color_cycle=["tab:blue"], 
    title = None
)

fig, ax = psf.common.plot.plot_station_property(
    station_dicts=[stations],  
    station_list=ox_stations,
    property_name="T",
    ylabel="T  [K]",
    color_cycle=["tab:blue"], 
    title = None
)

fig, ax = psf.common.plot.plot_station_property(
    station_dicts=[stations],  
    station_list=ox_stations,
    property_name="mdot",
    ylabel="mdot  [kg/s]",
    color_cycle=["tab:blue"],
    title = None
)

fu_PT_stations = [
    "fu_engine_in", "fu_pump_out", 
    "fu_turbine_in", "fu_chamber_in"
]
#subset = ['fu_tank', 'stage1_pump_in', 'pump_interstage2']
fig, ax = psf.common.plot.plot_PT_diagram(
    station_dicts=[stations],
    station_list=fu_PT_stations,
    fluid_name='Methane',
    title=None,
    color_cycle=["tab:red"],
    labels=["Fuel Path"],
    legend_loc="lower right",
    sat_points=300, 
    scale='log',
    annotate_ha=["left", "left", "center", "left"],
    annotate_va=["top", "bottom", "bottom", "top"]
)

ox_PT_stations = [
    "ox_engine_in", "ox_pump_out", 
    "ox_turbine_in", "ox_chamber_in"
]
fig, ax = psf.common.plot.plot_PT_diagram(
    station_dicts=[stations],
    station_list=ox_PT_stations,
    fluid_name='Oxygen',
    title=None,
    color_cycle=["tab:blue"],
    labels=["Oxidizer Path"],
    legend_loc="lower right",
    sat_points=300, 
    scale='log',
    annotate_ha=["left", "left", "center", "left"],
    annotate_va=["top", "bottom", "bottom", "top"]
)

psf.regen.plot.plot_temperature_profile(fu_cooling_data, thrust_chamber, circuit_index=0, x_query=0.5)
psf.regen.plot.plot_temperature_profile(ox_cooling_data_1, thrust_chamber, circuit_index=1, x_query=-0.5)
psf.regen.plot.plot_coolant_pressure(fu_cooling_data, ox_cooling_data_1, ox_cooling_data_2, static=True, stagnation=True)
psf.regen.plot.plot_coolant_temperature(fu_cooling_data, ox_cooling_data_1, ox_cooling_data_2)
psf.regen.plot.plot_wall_temperature(fu_cooling_data, ox_cooling_data_1, ox_cooling_data_2, 
                                     circuits=[thrust_chamber.cooling_circuit_group.circuits[0], thrust_chamber.cooling_circuit_group.circuits[1], thrust_chamber.cooling_circuit_group.circuits[2]],
                                     label_interfaces=[True, False, False], 
                                     plot_hot=True, 
                                     plot_interfaces=True, 
                                     plot_coolant_wall=True,
                                     annotate_va=["bottom", "bottom", "bottom"],
                                     annotate_ha=["left", "center", "center"]
                                     )
psf.regen.plot.plot_velocity(fu_cooling_data, ox_cooling_data_1, ox_cooling_data_2)
psf.regen.plot.plot_heat_flux(fu_cooling_data, ox_cooling_data_1, ox_cooling_data_2)

plt.show()

network_name = "methane_engine_network.html"
network_path = os.path.join(script_dir, network_name)

psf.common.plot.plot_engine_network(
        net,
        output_file=network_name,
        output_path=script_dir,
        notebook=False,
        height="1080px",
        width="100%", 
        mass_flow_based_arrows=True, 
        station_mode="both")

webview.create_window('Engine Network', network_path)
webview.start()

