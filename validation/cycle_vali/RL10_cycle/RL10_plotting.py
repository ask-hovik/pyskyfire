# -------------------------------------------------------------
# RL10_plotting.py
# -------------------------------------------------------------

import os
import json
import pyskyfire as psf
import matplotlib.pyplot as plt
import webview
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
in_path = os.path.join(script_dir, "engine_sizer_res.pkl")

res = psf.common.Results.load(in_path)

net = res.net
stations = res.stations
signals = res.signals
residuals = res.residuals
input_params = res.input_params
thrust_chamber = input_params["thrust_chamber"]
cooling_data = res.block_results["regen_full_pass"]
cooling_data_2 = res.block_results["regen_half_pass"]

psf.regen.plot.plot_contour([thrust_chamber.contour])
psf.common.plot.plot_residual_history(residuals)

# load the reference json
ref_path = os.path.join(script_dir, "RL10_station_data.json")
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


fig1, axes1 = plt.subplots(3, 1, sharex=True, figsize=(8, 10))

# fuel plots
psf.common.plot.plot_station_property(
    station_dicts=[stations, ref_fu_stations],  
    station_list=fu_stations,
    property_name="p",
    ax=axes1[0],
    title = "Fuel-Side Property vs. Station",
    labels=["Pyskyfire", "NASA, Binder"],
    ylabel="p  (Pa)"
)

psf.common.plot.plot_station_property(
    station_dicts=[stations, ref_fu_stations],  
    station_list=fu_stations,
    property_name="T",
    title=False,
    ax=axes1[1],
    labels=False,
    ylabel="T  (K)"
)

psf.common.plot.plot_station_property(
    station_dicts=[stations, ref_fu_stations],  
    station_list=fu_stations,
    property_name="mdot",
    ax=axes1[2],
    title=False,
    labels=False,
    ylabel=r"$\mathrm{\dot{m}}$  (kg/s)"
)

fig2, axes2 = plt.subplots(3, 1, sharex=True, figsize=(8, 10))

# ox plots
fig, ax = psf.common.plot.plot_station_property(
    station_dicts=[stations, ref_ox_stations],  
    station_list=ox_stations,
    property_name="p",
    title="Ox-Side Property vs. Station",
    ax=axes2[0],
    labels=["Pyskyfire", "NASA, Binder"],
    ylabel="p (Pa)"
)

fig, ax = psf.common.plot.plot_station_property(
    station_dicts=[stations, ref_ox_stations],  
    station_list=ox_stations,
    property_name="T",
    ax=axes2[1],
    title=False,
    labels=False,
    ylabel="T (K)"
)

fig, ax = psf.common.plot.plot_station_property(
    station_dicts=[stations, ref_ox_stations],  
    station_list=ox_stations,
    property_name="mdot",
    ax=axes2[2],
    title=False,
    labels=False,
    ylim=[13.9, 14.1],
    ylabel=r"$\mathrm{\dot{m}}$ (kg/s)"
)

subset = ['fu_tank', 'stage1_pump_in', 'pump_interstage2']
fig, ax = psf.common.plot.plot_PT_diagram(
    station_dicts=[stations],
    station_list=fu_stations,
    fluid_name='Hydrogen',
    title='Fuel‐side P–T path',
    sat_points=300, 
    scale='log'
)

plt.style.use('ggplot')
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Latin Modern Roman'],
    'font.size': 10,
    'axes.titleweight': 'bold',
    'axes.labelsize': 11,
    'axes.labelweight': 'normal',
    'legend.frameon': True,
    'legend.edgecolor': 'black',
    'grid.linestyle': '--',
    'grid.alpha': 0.7, 
    'savefig.dpi': 600
})

# -------------------------------
# Conversion Functions
# -------------------------------

def in_to_m(x):
    return x * 0.0254

def rankine_to_kelvin(T):
    return T * (5/9)

def psia_to_pa(p):
    return p * 6894.75729

def heat_flux_btu_to_w_m2(q):
    return q * (1055.06 / 0.00064516)

# -------------------------------
# Load and Convert Paper Data
# -------------------------------

def load_and_convert_paper_data(path):
    with open(path, 'r') as f:
        data = json.load(f)
    for model, ds in data.items():
        xi = ds['data']['x']
        ds['x_SI'] = [in_to_m(x) for x in xi]
        yi = ds['data']['y']
        unit = ds['units']['y']
        if unit == 'Rankine':
            ds['y_SI'] = [rankine_to_kelvin(y) for y in yi]
        elif unit == 'psia':
            ds['y_SI'] = [psia_to_pa(y) for y in yi]
        elif 'Btu' in unit:
            ds['y_SI'] = [heat_flux_btu_to_w_m2(y) for y in yi]
        else:
            ds['y_SI'] = yi
    return data

# -------------------------------
# File Paths and Data
# -------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, 'RL10_data')
files = {
    'wall_temp': 'wall_temperature.json',
    'coolant_pressure': 'coolant_static_pressure.json',
    'coolant_temp': 'coolant_static_temperature.json',
    'heat_flux': 'heat_flux.json',
}
paper_data = {k: load_and_convert_paper_data(os.path.join(data_dir, fn))
              for k, fn in files.items()}




def percent_diff(sim, ref):
    """Return (sim – ref)/ref × 100 [%]."""
    return 100.0 * (sim - ref) / ref

# -------------------------------
# Plotting Helper
# -------------------------------

def plot_comparison(x_list, y_list, labels,
                    paper_datasets, xlabel, ylabel,
                    title=None, xlim=None):
    """
    Overlay reference (paper) datasets and simulation results.

    Parameters
    ----------
    x_list : list of array-like
        X-values for simulation curves.
    y_list : list of array-like
        Y-values for simulation curves.
    labels : list of str
        Legend labels for simulation curves.
    paper_datasets : dict
        Mapping of reference dataset names to dicts with 'x_SI' and 'y_SI'.
    xlabel : str
    ylabel : str
    title : str, optional
        Plot title. If None, no title is set.
    xlim : tuple, optional
        X-axis limits.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    # Overlay paper data first
    ref_markers = ['^', 's', 'D', 'v', '>', '<', 'h', 'p', '*']
    for (model, ds), m in zip(paper_datasets.items(), ref_markers):
        ax.plot(
            ds['x_SI'], ds['y_SI'],
            marker=m, linestyle='-', color='black',
            markerfacecolor='white', markeredgecolor='black', markersize=6,
            label=f'{model.replace("_", " ").title()}'
        )

    # Plot simulation curves on top
    for x, y, lbl in zip(x_list, y_list, labels):
        ax.plot(x, y, linestyle='--', linewidth=2, label=lbl, color="tab:red")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if xlim:
        ax.set_xlim(xlim)
    ax.grid(True, which='both')

    # Combine legend entries and remove duplicates
    handles, lbls = ax.get_legend_handles_labels()
    by_label = dict(zip(lbls, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='best')

    plt.tight_layout()
    

# -------------------------------
# Extract and Plot Data (cooling_data only)
# -------------------------------

# 1) Wall Temperature Comparison (hot wall only)
T1 = np.array(cooling_data['T'])        # shape (n_x, n_layers)
hot_wall = T1[:, -1]                    # last column is hot wall
plot_comparison(
    [cooling_data['x']], [hot_wall], ['Sim: Hot Wall Temperature'],
    paper_data['wall_temp'],
    xlabel='Axial Position (m)',
    ylabel='Temperature (K)'
)

# 2) Coolant Static Pressure (Pa)
p = np.array(cooling_data['p_static'])
plot_comparison(
    [cooling_data['x']], [p], ['Pyskyfire'],
    paper_data['coolant_pressure'],
    xlabel='Axial Position (m)',
    ylabel='Pressure (Pa)'
)

# 3) Coolant Velocity
v = cooling_data['velocity']
plot_comparison(
    [cooling_data['x']], [v], ['Pyskyfire'],
    {},
    xlabel='Axial Position (m)',
    ylabel='Velocity (m/s)'
)

# 4) Coolant Temperature Comparison
x_list = [cooling_data['x'], cooling_data_2['x']]
y_list = [cooling_data['T_static'], cooling_data_2['T_static']]
labels = ['Pyskyfire', 'Pyskyfire']
plot_comparison(
    x_list, y_list, labels,
    paper_data['coolant_temp'],
    xlabel='Axial Position (m)',
    ylabel='Temperature (K)'
)

# 5) Heat Flux Comparison
q = cooling_data['dQ_dA']
plot_comparison(
    [cooling_data['x']], [q], ['Pyskyfire'],
    paper_data['heat_flux'],
    xlabel='Axial Position (m)',
    ylabel='Heat Flux (W/m²)'
)

psf.regen.plot.plot_precomputed_properties(thrust_chamber)

def plot_residuals(cooling_data, p=2):
    """
    Plot solver residuals in a paper-friendly style:
    1) Global nonlinear-residual history (semilog).
    2) Final per-cell residual distribution (semilog scatter).

    Parameters
    ----------
    cooling_data : dict
        Output from steady_heating_analysis_2 / solve_heat_exchanger_newton,
        must contain key "residuals" = (global_R, final_R).
    p : {2, np.inf}
        Order of the global norm (2 for RMS, np.inf for max).
    """
    import numpy as _np
    import matplotlib.pyplot as _plt

    global_R, final_R = cooling_data.get("residuals", (None, None))
    if global_R is None or final_R is None:
        print("No residual data to plot. Run with log_residuals=True.")
        return

    # 1) Global residual history
    fig, ax = _plt.subplots(figsize=(6, 4))
    ax.semilogy(global_R, marker='o', linestyle='-')
    ax.set_xlabel("Newton Iteration")
    norm_symbol = r"\infty" if p == _np.inf else str(p)
    ax.set_ylabel(r"Global Residual $\|R^{(k)}\|_{" + norm_symbol + r"}$")
    ax.set_title("Nonlinear Solver Convergence")
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    fig.tight_layout()

    # 2) Final per-cell residual distribution
    fig, ax = _plt.subplots(figsize=(6, 4))
    cell_idx = _np.arange(len(final_R))
    ax.semilogy(cell_idx, final_R, marker='o', linestyle='None')
    ax.set_xlabel("Cell Index")
    ax.set_ylabel(r"Final Residual $\|R_i\|_2$")
    ax.set_title("Per-Cell Final Residuals")
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    fig.tight_layout()

    
plot_residuals(cooling_data)

psf.regen.plot.plot_temperature_profile(results=cooling_data, thrust_chamber=thrust_chamber, circuit_index=0, x_query=-0.15)

plt.show()

network_name = "engine_network_reversed.html"
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

