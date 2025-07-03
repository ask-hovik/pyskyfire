import pyskyfire as psf
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import time


# Engine inputs 
p_c = 32.7501e5
p_e = 0.0377e5 
MR = 5.0
L_star = 0.95
r_c = 0.123
fu = {"cea": "LH2", 
      "cantera": "H2",
      "coolprop": "hydrogen"}
ox = {"cea": "LOX", 
      "cantera": "O2",
      "coolprop": "oxygen"}
theta_conv = 25
R_1f = 1.5
R_2f = 3
R_3f = 0.5
length_fraction = 0.713
F = 73.4e3
T_coolant_in = 33 
p_coolant_in = 69e5 

# Find the optimal values for a given set of inputs
optimals = psf.skycea.OptimalValues(ox["cea"], fu["cea"], F, MR, p_c, p_e, L_star)

V_c = optimals.V_c_opt
r_t = optimals.r_t_opt
eps = optimals.eps_opt

# Generate the contour coordinates
RL10_xs, RL10_rs = psf.regen.contour_2.get_contour_2(V_c = V_c, 
                                            r_t = r_t, 
                                            area_ratio = eps, 
                                            r_c = r_c, 
                                            theta_conv = theta_conv,
                                            nozzle = "rao", 
                                            R_1f=R_1f,
                                            R_2f=R_2f,
                                            R_3f=R_3f,
                                            length_fraction = length_fraction
                                            )

# Make a contour object using the 
RL10_contour = psf.regen.Contour(RL10_xs, RL10_rs, name = "Replicate Contour")

#psf.regen.contour.plot_theta_vs_epsilon()

wall = psf.regen.Wall(material = psf.common.material.StainlessSteel304, thickness = 0.3e-3) 

wall_group = psf.regen.WallGroup(walls=[wall])

channel_height_fn = psf.regen.make_channel_height_fn(
    contour=RL10_contour, 
    region_fractions=[-1.0, 0.25, 1.0], 
    flat_heights= [0.0032, 0.00134], 
    pinch_factors= [0.7, -5.0], 
    transition_widths=[0.1]
) # total coolant volume should be ca 0.015831543m3

cross_section = psf.regen.CrossSectionRounded()
LH2_transport = psf.skycea.CoolantTransport(fu["coolprop"])

half_pass = psf.regen.CoolingCircuit(name="Half Pass", 
                                     contour=RL10_contour, 
                                     coolant_transport=LH2_transport, 
                                     cross_section=cross_section, 
                                     span = [0.25, 1.0], 
                                     placement=psf.regen.SurfacePlacement(n_channel_positions=180), 
                                     channel_height=channel_height_fn)

full_pass = psf.regen.CoolingCircuit(name= "Full Pass",
                                     contour=RL10_contour, 
                                     coolant_transport=LH2_transport, 
                                     cross_section=cross_section, 
                                     span = [1.0, -1.0], 
                                     placement=psf.regen.SurfacePlacement(n_channel_positions=180), 
                                     channel_height=channel_height_fn)

cooling_circuit_group = psf.regen.CoolingCircuitGroup(circuit_list=[full_pass, half_pass])

combustion_transport = psf.skycea.CombustionTransport(ox=ox["cantera"], 
                                                       fu=fu["cantera"],
                                                       p_c=p_c,
                                                       p_e=p_e,
                                                       MR=MR,
                                                       T_fuel_in=111, 
                                                       T_ox_in=91, 
                                                       mode="hybrid")

combustion_transport.set_essentials(mdot=optimals.mdot, mdot_fu=optimals.mdot_fu, mdot_ox=optimals.mdot_ox)




thrust_chamber = psf.regen.ThrustChamber(contour=RL10_contour, 
                                         wall_group=wall_group,
                                         combustion_transport=combustion_transport,  
                                         cooling_circuit_group=cooling_circuit_group,
                                         roughness=0.015e-3)


"""psf.regen.plot.visualize_cooling_channels_vispy_6(
    thrust_chamber,
    draw_section_line=True,
    clip_plane=((-0.15,0,0),(0,1,0)),
    line_color=(0,0,0,1),
    draw_only_section_line=True, 
    show_axis=False
)""" # this is the one


canvas, view = psf.regen.plot.plot_engine_3D(thrust_chamber, show_axis=False)
#psf.regen.plot.save_canvas_png(canvas)

#psf.regen.plot.plot_cross_section(thrust_chamber, circuit_index = 0)
#psf.regen.plot.plot_precomputed_properties(thrust_chamber)

mdot_fu = combustion_transport.mdot_fu
boundary_conditions = psf.regen.BoundaryConditions(T_coolant_in = T_coolant_in, p_coolant_in = p_coolant_in, mdot_coolant = mdot_fu)

# direct running of project
start_time = time.time()
cooling_data_2 = psf.regen.steady_heating_analysis(thrust_chamber, n_nodes = 30, circuit_index=1, boundary_conditions=boundary_conditions, solver="newton", output=True)
T_out = cooling_data_2['T_static'][-1] # note that the endpoint of the analysis is always at the end index, no matter the direction of coolant flow
p_out = cooling_data_2['p_static'][-1]
boundary_conditions_2 = psf.regen.BoundaryConditions(T_coolant_in = T_out, p_coolant_in = p_out, mdot_coolant = mdot_fu)
cooling_data = psf.regen.steady_heating_analysis(thrust_chamber, n_nodes = 100, circuit_index=0, boundary_conditions=boundary_conditions_2, solver="newton", output=True)
end_time = time.time()
elapsed = end_time - start_time
print(f"Process took {elapsed:.4f} seconds")

print(f"Max heat flux: {max(cooling_data["dQ_dA"])}")




# -------------------------------
# Global Style Configuration
# -------------------------------
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


# -------------------------------------------------
# ❷ Gather simulation scalars
# -------------------------------------------------
sim_peak_q     = float(np.max(cooling_data["dQ_dA"]))   # W m-2
sim_T_out      = float(cooling_data["T_static"][-1])    # K
sim_p_out      = float(cooling_data["p_static"][-1])    # Pa

# -------------------------------------------------
# ❸ Extract reference scalars
#     (one value per curve in the JSON file)
# -------------------------------------------------
ref_peaks_q  = {m: max(ds["y_SI"])      for m, ds in paper_data["heat_flux"].items()}
ref_T_outs   = {m: ds["y_SI"][-1]       for m, ds in paper_data["coolant_temp"].items()}
ref_p_outs   = {m: ds["y_SI"][-1]       for m, ds in paper_data["coolant_pressure"].items()}

# -------------------------------------------------
# ❹ Compute percent differences
# -------------------------------------------------
pct_peak_q =  {m: percent_diff(sim_peak_q, r)  for m, r in ref_peaks_q.items()}
pct_T_out  =  {m: percent_diff(sim_T_out,  r)  for m, r in ref_T_outs.items()}
pct_p_out  =  {m: percent_diff(sim_p_out,  r)  for m, r in ref_p_outs.items()}

# -------------------------------------------------
# ❺ Pretty-print summary
# -------------------------------------------------
print("─" * 60)
print(f"Simulation peak heat-flux : {sim_peak_q:12.3e}  W/m²")
print(f"Simulation outlet T       : {sim_T_out:12.3f}    K")
print(f"Simulation outlet p       : {sim_p_out:12.3f}    Pa")
print("─" * 60)
print("Percent difference versus each reference curve:")
for model in sorted(set(pct_peak_q) | set(pct_T_out) | set(pct_p_out)):
    dq = pct_peak_q.get(model, None)
    dT = pct_T_out.get(model, None)
    dp = pct_p_out.get(model, None)
    print(f"{model:25s}  "
          f"Δq_peak = {dq}%   "
          f"ΔT_out = {dT}%   "
          f"Δp_out = {dp}%")
print("─" * 60)
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

    

# Example usage:
#   data = steady_heating_analysis_2(..., log_residuals=True)
#   plot_residuals(data)


# Example usage:
#   data = steady_heating_analysis_2(..., log_residuals=True)\#   plot_residuals(data)

plot_residuals(cooling_data)

psf.regen.plot.plot_temperature_profile(results=cooling_data, thrust_chamber=thrust_chamber, circuit_index=0, x_query=-0.15)

plt.show()