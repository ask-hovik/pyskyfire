import os
import json
import pyskyfire as psf
import numpy as np

# ===============
# Loading results
# ===============

script_dir = os.path.dirname(os.path.abspath(__file__))
in_path = os.path.join(script_dir, "sizer_results.pkl")

res = psf.common.Results.load(in_path)

net = res.net
stations = res.stations
signals = res.signals
residuals = res.residuals
params = res.input_params
thrust_chamber = params["thrust_chamber"]
cooling_data_b = res.block_results["regen_full_pass"]
cooling_data_a = res.block_results["regen_half_pass"]

cooling_data_a["name"] = "Pyskyfire"
cooling_data_b["name"] = "Pyskyfire" # TODO: maybe bake this into the simulation process itself?



# fetching interesting values from aerothermo TODO: It seems a bit unpractical that these values are only stored in aerothermodynamics self

optimal_values = dict(
    MR              = thrust_chamber.combustion_transport.MR,
    p_c             = thrust_chamber.combustion_transport.p_c,
    T_c             = thrust_chamber.combustion_transport.T_c,
    F               = thrust_chamber.combustion_transport.F,
    eps             = thrust_chamber.combustion_transport.eps,
    L_star          = thrust_chamber.combustion_transport.L_star,
    c_star          = thrust_chamber.combustion_transport.c_star,
    p_amb           = thrust_chamber.combustion_transport.p_amb,
    Isp_ideal_amb   = thrust_chamber.combustion_transport.Isp_ideal_amb,
    Isp_vac         = thrust_chamber.combustion_transport.Isp_vac,
    Isp_amb         = thrust_chamber.combustion_transport.Isp_amb,
    Isp_SL          = thrust_chamber.combustion_transport.Isp_SL,
    CF_vac          = thrust_chamber.combustion_transport.CF_vac,
    CF_amb          = thrust_chamber.combustion_transport.CF_amb,
    CF_SL           = thrust_chamber.combustion_transport.CF_SL,
    mdot            = thrust_chamber.combustion_transport.mdot,
    mdot_fu         = thrust_chamber.combustion_transport.mdot_fu,
    mdot_ox         = thrust_chamber.combustion_transport.mdot_ox,
    t_stay          = thrust_chamber.combustion_transport.t_stay,
    A_t             = thrust_chamber.combustion_transport.A_t,
    A_e             = thrust_chamber.combustion_transport.A_e,
    r_t             = thrust_chamber.combustion_transport.r_t,
    r_e             = thrust_chamber.combustion_transport.r_e,
    V_c             = thrust_chamber.combustion_transport.V_c,
    npts            = thrust_chamber.combustion_transport.npts,
)


# ========
# Plotting
# ========

# Generate report
report = psf.viz.Report("RL10 engine")

# -------------- Engine Sizer Results ----------------
tab_sizer = report.add_tab("Engine Sizer")
fig_residual_history = psf.viz.PlotResidualHistory(residuals)
tab_sizer.add_figure(fig_residual_history, caption="Residual development over iterations")

# load the reference json
ref_path = os.path.join(script_dir, "reference_station_data.json")
with open(ref_path, "r") as f:          # or just with open(ref_path) as f:
    fu_data = json.load(f)["engine_condition_fuel"]

with open(ref_path, "r") as g:          # or just with open(ref_path) as f:
    ox_data = json.load(g)["engine_condition_oxidizer"]

# mapping between json station names and model names
fu_paper2code = {
    "fuel_tank"            : "fu_engine_in",
    "fuel_pump_inlet"      : "stage1_pump_in",
    "fuel_pump_interstage" : "pump_interstage2",
    "fuel_pump_discharge"  : "regen_duct_in",
    "cooling_jacket_inlet" : "regen_in",
    "cooling_jacket_exit"  : "regen_out",     
    "turbine_inlet"        : "turbine_in",
    "turbine_housing_exit" : "turbine_out",
    "fuel_injector_plenum" : "fu_injector_plenum",

}

ox_paper2code = {
    "oxygen_tank"           : "ox_engine_in",
    "ox_pump_inlet" : "ox_pump_in",
    "ox_pump_discharge" : "ox_duct_in",
    "ox_injector_plenum" : "ox_injector_plenum"

}

# rename json stations and build compatible dictionary
ref_fu_stations = {}
ref_ox_stations = {}
for paper_name, code_name in fu_paper2code.items():
    d = fu_data[paper_name]
    ref_fu_stations[code_name] = {
        "p"    : d["total_pressure"],    # Pa
        "T"    : d["total_temperature"], # K
        "mdot" : d["mass_flow"],         # kg/s  (None will just propagate)
    }

for paper_name, code_name in ox_paper2code.items():
    e = ox_data[paper_name]
    ref_ox_stations[code_name] = {
        "p"    : e["total_pressure"],    # Pa
        "T"    : e["total_temperature"], # K
        "mdot" : e["mass_flow"],         # kg/s  (None will just propagate)
    }


# station list now available both from reference data and model
fu_stations = [
    "fu_engine_in", "stage1_pump_in", "pump_interstage2",
    "regen_duct_in", "regen_in", "regen_out",
    "turbine_in", "turbine_out", "fu_injector_plenum"
]

ox_stations = [
    "ox_engine_in", "ox_pump_in", "ox_duct_in", "ox_injector_plenum"
]

fig_fuel_side_pressure = psf.viz.PlotStationProperty(station_dicts=[stations, ref_fu_stations],
                            station_list=fu_stations,
                            property_name="p",
                            labels=["Pyskyfire", "NASA, Binder"]
                            )

fig_fuel_side_temperature = psf.viz.PlotStationProperty(station_dicts=[stations, ref_fu_stations],
                            station_list=fu_stations,
                            property_name="T",
                            labels=["Pyskyfire", "NASA, Binder"]
                            )

fig_fuel_side_mdot = psf.viz.PlotStationProperty(station_dicts=[stations, ref_fu_stations],
                            station_list=fu_stations,
                            property_name="mdot",
                            labels=["Pyskyfire", "NASA, Binder"]
                            )

fig_ox_side_pressure = psf.viz.PlotStationProperty(station_dicts=[stations, ref_ox_stations],
                            station_list=ox_stations,
                            property_name="p",
                            labels=["Pyskyfire", "NASA, Binder"]
                            )

fig_ox_side_temperature = psf.viz.PlotStationProperty(station_dicts=[stations, ref_ox_stations],
                            station_list=ox_stations,
                            property_name="T",
                            labels=["Pyskyfire", "NASA, Binder"]
                            )

fig_ox_side_mdot = psf.viz.PlotStationProperty(station_dicts=[stations, ref_ox_stations],
                            station_list=ox_stations,
                            property_name="mdot",
                            labels=["Pyskyfire", "NASA, Binder"]
                            )

tab_sizer.add_figure(fig_fuel_side_pressure)
tab_sizer.add_figure(fig_fuel_side_temperature)
tab_sizer.add_figure(fig_fuel_side_mdot)
tab_sizer.add_figure(fig_ox_side_pressure)
tab_sizer.add_figure(fig_ox_side_temperature)
tab_sizer.add_figure(fig_ox_side_mdot)

fig_fuel_side_PT_diagram = psf.viz.PlotPTDiagram(station_dicts=[stations],
                                                 station_list=fu_stations,
                                                 fluid_name="hydrogen",
                                                 title='Fuel-side P-T path',
                                                 scale="linear") #TODO: annotations in log mode don't work properly

fig_ox_side_PT_diagram = psf.viz.PlotPTDiagram(station_dicts=[stations],
                                                 station_list=ox_stations,
                                                 fluid_name="oxygen",
                                                 title='Ox-side P-T path',
                                                 scale="linear") 

tab_sizer.add_figure(fig_fuel_side_PT_diagram)
tab_sizer.add_figure(fig_ox_side_PT_diagram)

tab_network = report.add_tab("Engine Network")
"""fig_network = psf.viz.PlotEngineNetwork(engine_network=net, 
                                        title="Engine Network",
                                        station_mode="values", 
                                        mass_flow_based_arrows=True)"""

html_network = psf.viz.render_engine_network(
    net,   # or thrust_chamber.network
    station_mode="values",           # "name" | "values" | "both" | "hidden"
    mass_flow_based_arrows=True,
    edge_length=20,
    physics_settings="default",      # or pass a dict of barnes_hut params
    #height="900px",
)

tab_network.add_raw_html(html_network)



# --------------  Engine parameters -----------------
tab_params = report.add_tab("Parameters")
tab_params.add_table(params, caption="Input Parameters", key_title="Parameter", value_title="Value", precision = 3)
tab_params.add_table(optimal_values, caption="Optimal Values", key_title="Parameters", value_title="Value", precision = 3)

# ---------------- Engine Overview ------------------
tab1 = report.add_tab("Engine Overview")
tab1.add_text(
"This is the report outlining the Pyskyfire simulation of the regenerative cooling of the RL10 engine. " \
"Reference data is largely collected from NASA, Binder, 1997, RL10 Modelling Project.")

# --- RL10 contour from JSON ---
json_path = os.path.join(script_dir, "reference_RL10_contour.json")
with open(json_path) as f:
    d = json.load(f)
xs = np.asarray(d["xs"], dtype=float)
rs = np.asarray(d["rs"], dtype=float)
actual_RL10_contour = psf.regen.Contour(xs, rs, name="Actual RL10 contour")
fig_engine_contour = psf.viz.PlotContour(thrust_chamber.contour, actual_RL10_contour)

# Generate STL (a bit heavy so if nothing changed could be an idea to comment out TODO: lighten the precompute paths logic)
#stl_path = os.path.join(script_dir, "engine_channels") # TODO: clunky
#psf.viz.make_engine_gmsh(thrust_chamber, filename=stl_path, display_channels=15)
# TODO: Why is the engine oriented along the z-axis? Where on earth did I reorient it from the x-axis?? Fuck
stl_in = os.path.join(script_dir, "engine_channels.stl")
fig_cooling_channel_stl = psf.viz.EmbedSTL(stl_in, color="#8ab7ff", opacity=1.0, show_wireframe=False)
tab1.add_figure(fig_cooling_channel_stl)

tab1.add_figure(fig_engine_contour, caption="RL10 engine contour")


def load_reference(path: str, prop_key: str):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    out = []
    for model, ds in data.items():
        x = np.asarray(ds["x"], dtype=float)
        y = np.asarray(ds[prop_key], dtype=float)
        #order = np.argsort(x)
        #x, y = x[order], y[order]

        cd = {"x": x, "name": ds.get("name", model)}
        if prop_key == "T_hot_wall":
            cd["T"] = y.reshape(-1, 1)
        else:
            cd[prop_key] = y
        out.append(cd)

    return out

# ======== tab 2 - cooling data =========
tab2 = report.add_tab("Cooling Data")

# -------- Wall temperature comparison ----------
reference_wall_temp_path = os.path.join(script_dir, "reference_wall_temperature.json")
reference_wall_temperature_data = load_reference(path=reference_wall_temp_path, prop_key="T_hot_wall")

fig_wall_temp_comparison = psf.viz.PlotWallTemperature(*reference_wall_temperature_data,
                                            cooling_data_b, 
                                            template="plotly_white"
                                            )

# Update traces and legends 
fig_wall_temp_comparison.update_traces(
    name="New System Model", # Name on legend
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="diamond", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(legendgroup="New System Model", name="Hot Wall"), # which one are you editing
    )

fig_wall_temp_comparison.update_traces(
    name="P&W Simulation",
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="triangle-up", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(legendgroup="P&W Simulation", name="Hot Wall"),
    )

fig_wall_temp_comparison.update_traces(
    name="RTE Model",
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="square", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(legendgroup="RTE Model", name="Hot Wall"),
    )

fig_wall_temp_comparison.update_traces(
    name="Pyskyfire",
    line=dict(color="red"),
    mode="lines",
    selector=dict(legendgroup="Pyskyfire", name="Hot Wall"),
    )

# add figure into report
tab2.add_figure(fig_wall_temp_comparison, caption="Comparison between simulations presented in NASA, Binder 1997, and Pyskyfire.")

fig_wall_temp_complete = psf.viz.PlotWallTemperature(cooling_data_b, plot_coolant_wall=True)
tab2.add_figure(fig_wall_temp_complete, caption="Hot wall and cold wall temperature predicted by Pyskyfire.")

# ------------ coolant temperature comparison ------------
reference_coolant_temp_path = os.path.join(script_dir, "reference_coolant_static_temperature.json")
reference_coolant_temperature_data = load_reference(path=reference_coolant_temp_path, prop_key="T_static")

fig_coolant_temperature = psf.viz.PlotCoolantTemperature(*reference_coolant_temperature_data, cooling_data_a, cooling_data_b, )

#for i, t in enumerate(fig_coolant_temperature.fig.data):
#    print(i, "type=", t.type, "name=", t.name, "legendgroup=", getattr(t, "legendgroup", None))

# Update traces and legends 
fig_coolant_temperature.update_traces(
    name="New System Model", # Name on legend
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="diamond", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="New System Model"), # which one are you editing
    )

fig_coolant_temperature.update_traces(
    name="P&W Simulation",
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="triangle-up", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="P&W Simulation"),
    )

fig_coolant_temperature.update_traces(
    name="RTE Model",
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="square", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="RTE Model"),
    )

fig_coolant_temperature.update_traces(
    name="Pyskyfire",
    line=dict(color="red"),
    mode="lines",
    selector=dict(name="Pyskyfire"),
    )

fig_coolant_temperature.update_traces(
    name="Test Data",
    mode="markers",
    marker=dict(symbol="circle", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="Test Data"),
    )

tab2.add_figure(fig_coolant_temperature)

# ----------- coolant pressure ------------
reference_coolant_pressure_path = os.path.join(script_dir, "reference_coolant_static_pressure.json")
reference_coolant_pressure_data = load_reference(path=reference_coolant_pressure_path, prop_key="p_static")
fig_coolant_pressure = psf.viz.PlotCoolantPressure(*reference_coolant_pressure_data, cooling_data_b, static = True, stagnation=False)


# Update traces and legends 
fig_coolant_pressure.update_traces(
    name="New System Model", # Name on legend
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="diamond", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="New System Model — Static"), # which one are you editing
    )

fig_coolant_pressure.update_traces(
    name="P&W Simulation",
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="triangle-up", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="P&W Simulation — Static"),
    )

fig_coolant_pressure.update_traces(
    name="RTE Model",
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="square", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="RTE Model — Static"),
    )

fig_coolant_pressure.update_traces(
    name="Pyskyfire",
    line=dict(color="red"),
    mode="lines",
    selector=dict(name="Pyskyfire — Static"),
    )

fig_coolant_pressure.update_traces(
    name="Test Data",
    mode="markers",
    marker=dict(symbol="circle", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="Test Data — Static"),
    )
tab2.add_figure(fig_coolant_pressure)

fig_coolant_pressure_complete = psf.viz.PlotCoolantPressure(cooling_data_b, static = True, stagnation=True)
tab2.add_figure(fig_coolant_pressure_complete)

# --------- heat flux ----------
reference_heat_flux_path = os.path.join(script_dir, "reference_heat_flux.json")
reference_heat_flux_data = load_reference(path=reference_heat_flux_path, prop_key="dQ_dA")
fig_heat_flux = psf.viz.PlotHeatFlux(*reference_heat_flux_data, cooling_data_b)


# Update traces and legends 
fig_heat_flux.update_traces(
    name="New System Model", # Name on legend
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="diamond", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="New System Model"), # which one are you editing
    )

fig_heat_flux.update_traces(
    name="P&W Simulation",
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="triangle-up", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="P&W Simulation"),
    )

fig_heat_flux.update_traces(
    name="RTE Model",
    line=dict(color="black"),
    mode="lines+markers",
    marker=dict(symbol="square", size=7, color="white", line=dict(color="black", width=1)),
    selector=dict(name="RTE Model"),
    )

fig_heat_flux.update_traces(
    name="Pyskyfire",
    line=dict(color="red"),
    mode="lines",
    selector=dict(name="Pyskyfire"),
    )

tab2.add_figure(fig_heat_flux)

# ----------- Velocity ----------
fig_velocity = psf.viz.PlotVelocity(cooling_data_b)

#for i, t in enumerate(fig_velocity.fig.data):
#    print(i, "type=", t.type, "name=", t.name, "legendgroup=", getattr(t, "legendgroup", None))

fig_velocity.update_traces(
    name="Pyskyfire",
    line=dict(color="red"),
    mode="lines",
    selector=dict(name="Pyskyfire"),
    )
tab2.add_figure(fig_velocity)

# Precomputed properties
tab3 = report.add_tab("Thrust Chamber Properties")
p1 = psf.viz.PlotCoolantArea(thrust_chamber, 1)
p2 = psf.viz.PlotHydraulicDiameter(thrust_chamber, 1)
p3 = psf.viz.PlotRadiusOfCurvature(thrust_chamber, 1)
p4 = psf.viz.PlotdAdxThermalHotGas(thrust_chamber, 1)
p5 = psf.viz.PlotdAdxThermalCoolant(thrust_chamber, 1)
p6 = psf.viz.PlotdAdxCoolantArea(thrust_chamber, 1)

tab3.add_figure(p1)
tab3.add_figure(p2)
tab3.add_figure(p3)
tab3.add_figure(p4)
tab3.add_figure(p5)
tab3.add_figure(p6)

# -------------- Built in charts ----------------
tab4 = report.add_tab("Lookup Tables")
theta_v_epsilon = psf.viz.PlotThetaVsEpsilon()
tab4.add_figure(theta_v_epsilon)

moody_diagram = psf.viz.PlotMoodyDiagram()
tab4.add_figure(moody_diagram)

# -------------- Thermal Gradient Charts ---------------
tab5 = report.add_tab("Thermal Gradient")
chamber_gradient = psf.viz.PlotTemperatureProfile(cooling_data_b, thrust_chamber, 1, -0.2)
throat_gradient = psf.viz.PlotTemperatureProfile(cooling_data_b, thrust_chamber, 1, 0)
exit_gradient = psf.viz.PlotTemperatureProfile(cooling_data_b, thrust_chamber, 1, 1.05)
tab5.add_figure(chamber_gradient)
tab5.add_figure(throat_gradient)
tab5.add_figure(exit_gradient)

# -------------- Combustion Transport ----------------
tab6 = report.add_tab("Combustion")
M = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="M")
gamma = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="gamma")
T = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="T")
p = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="p")
h = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="h")
cp = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="cp")
k = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="k")
mu = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="mu")
Pr = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="Pr")
rho = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="rho")
a = psf.viz.PlotTransportProperty(thrust_chamber.combustion_transport, prop="a")

tab6.add_figure(M)
tab6.add_figure(gamma)
tab6.add_figure(T)
tab6.add_figure(p)
tab6.add_figure(h)
tab6.add_figure(cp)
tab6.add_figure(k)
tab6.add_figure(mu)
tab6.add_figure(Pr)
tab6.add_figure(rho)
tab6.add_figure(a)

# Save report to the same directory as the script
out_path = os.path.join(script_dir, "sizer_report.html")
report.save_html(out_path)
print(f"Report saved to {out_path}")

