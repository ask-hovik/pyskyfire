import os
import json
import pyskyfire as psf
import numpy as np

# ===============
# Loading results
# ===============

script_dir = os.path.dirname(os.path.abspath(__file__)) 
in_path = os.path.join(script_dir, "regen_results.pkl")

res = psf.common.Results.load(in_path)
params = res.params
thrust_chamber = res.thrust_chamber
cooling_data_a = res.cooling_data[0]
cooling_data_b = res.cooling_data[1]

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
out_path = os.path.join(script_dir, "regen_report.html")
report.save_html(out_path)
print(f"Report saved to {out_path}")

