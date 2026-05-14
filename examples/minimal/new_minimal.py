# ============================================================
# Welcome to the minimal example showing how to use Pyskyfire!
# ============================================================

import os
import pyskyfire as psf

# Input parameters for your engine
params = dict(
    # Chamber parameters
    p_c    = 50e5,
    F      = 5e3,
    eps    = 10,
    L_star = 1.2,
    MR     = 2.8,
    AR_c   = 1.8,

    # Propellants
    cea_fu = psf.common.Fluid(type="fuel", propellants=["C2H5OH"], fractions=[1.0]),
    cea_ox = psf.common.Fluid(type="oxidizer", propellants=["N2O"], fractions=[1.0]),
    coolprop_fu = psf.common.Fluid(type="fuel", propellants=["ethanol"], fractions=[1.0]),

    # Coolant inlet conditions
    T_coolant_in = 298.15,
    p_coolant_in = 23e5,

    # Cooling system properties
    material = psf.common.solids.StainlessSteel304,
    wall_thickness = 0.5e-3,
    n_channels = 60,
    blockage_ratio = 0.0,
    roughness_height = 10e-6,

    # Film cooling inputs
    film_x_fraction = -0.85,
    film_mdot = 0.02,
    film_injection_perimeter = 0.10,
    film_liquid_absorptivity = 0.2,
    film_turbulence_intensity = 0.02,
    film_mole_fraction_CO2 = 0.1,
    film_mole_fraction_H2O = 0.4,
)

# Combustion / aerothermodynamics
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

# Coolant transport
coolant_transport = psf.skycea.CoolantTransport(params["coolprop_fu"])

# Contour
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

# Wall
wall = psf.regen.Wall(
    material=params["material"],
    thickness=params["wall_thickness"],
)

# Channel cross-section
cross_section_squared = psf.regen.CrossSectionSquared(
    blockage_ratio=params["blockage_ratio"]
)

def channel_height_function(x):
    return 2e-3

def helix_fn(x):
    return 45 * 3.14 / 180

helical_placement = psf.regen.SurfacePlacement(
    n_channel_positions=59,
    helix_angle=helix_fn,
)

cooling_circuit = psf.regen.CoolingCircuit(
    name="Cooling Pass",
    contour=contour,
    coolant_transport=coolant_transport,
    cross_section=cross_section_squared,
    span=[1.0, -1.0],
    placement=helical_placement,
    walls=[wall],
    roughness=params["roughness_height"],
    channel_height=channel_height_function,
)

film_cooling = psf.regen.FilmCooling(
    coolant_transport=coolant_transport,
    x_fraction=params["film_x_fraction"],
    coolant_mass_flow_rate=params["film_mdot"],
    film_injection_perimeter=params["film_injection_perimeter"],
    liquid_absorptivity=params["film_liquid_absorptivity"],
    turbulence_intensity=params["film_turbulence_intensity"],
    mole_fraction_CO2=params["film_mole_fraction_CO2"],
    mole_fraction_H2O=params["film_mole_fraction_H2O"],
)

thrust_chamber = psf.regen.ThrustChamber(
    contour=contour,
    combustion_transport=aerothermodynamics,
    cooling_circuits=[cooling_circuit],
    n_nodes=150,
    film_cooling=film_cooling,
)

mdot_fu = aerothermodynamics.mdot_fu

boundary_conditions = psf.regen.BoundaryConditions(
    T_coolant_in=params["T_coolant_in"],
    p_coolant_in=params["p_coolant_in"],
    mdot_coolant=mdot_fu,
)

# ============================================================
# Film cooling simulation (run first)
# ============================================================
film_model = psf.regen.GrissonFilmCoolingModel(
    thrust_chamber=thrust_chamber,
    boundary_conditions=boundary_conditions
)

film_x = thrust_chamber.contour.xs
liquid_film_results, gaseous_film_results = film_model.solve(film_x)

print("\n=== Film cooling results ===")
print(f"Injection x: {thrust_chamber.film_cooling.x:.6f} m")
print(f"Initial film mdot: {thrust_chamber.film_cooling.coolant_mass_flow_rate:.6f} kg/s")

if liquid_film_results.x_dryout is None:
    print("Liquid film survives to the end of the domain.")
else:
    print(f"Liquid film dryout at x = {liquid_film_results.x_dryout:.6f} m")
    print(f"Gaseous film nodes: {len(gaseous_film_results.x)}")

# ============================================================
# Regenerative cooling simulation (run second)
# ============================================================
cooling_data = psf.regen.steady_heating_analysis(
    thrust_chamber,
    n_nodes=100,
    circuit_index=0,
    boundary_conditions=boundary_conditions,
    output=True,
)

# ============================================================
# Plotting / report
# ============================================================

plot_3d, viewer = psf.viz.make_engine_3d(thrust_chamber)
plot_3d.show()
del plot_3d 

script_dir = os.path.dirname(os.path.abspath(__file__))

print("Started generating report")
report = psf.viz.Report("Minimal Engine")

tab_params = report.add_tab("Parameters")
optimal_values = thrust_chamber.combustion_transport.optimum
tab_params.add_table(params, caption="Input Parameters", key_title="Parameter", value_title="Value", precision=3)
tab_params.add_table(optimal_values, caption="Optimal Values", key_title="Parameters", value_title="Value", precision=3)

tab_overview = report.add_tab("Engine Overview")
fig_engine_contour = psf.viz.PlotContour(thrust_chamber.contour)
tab_overview.add_figure(fig_engine_contour)

tab_film = report.add_tab("Film Cooling")
T_g_ref = thrust_chamber.combustion_transport.get_T(thrust_chamber.film_cooling.x)
p_ref = thrust_chamber.combustion_transport.get_p(thrust_chamber.film_cooling.x)
T_sat_ref = cooling_circuit.coolant_transport.get_T_sat(p_ref) if hasattr(cooling_circuit.coolant_transport, "get_T_sat") else None

fig_film_wall_temperature = psf.viz.PlotFilmWallTemperature(
    liquid_film_results,
    gaseous_film_results,
    T_g_ref=T_g_ref,
    T_sat_ref=T_sat_ref,
)
fig_film_dryout = psf.viz.PlotFilmDryout(liquid_film_results)
fig_film_heat_flux = psf.viz.PlotFilmHeatFlux(liquid_film_results)
fig_film_evaporation_rate = psf.viz.PlotFilmEvaporationRate(liquid_film_results)
fig_film_hconv = psf.viz.PlotFilmHeatTransferCoefficient(liquid_film_results, gaseous_film_results)
fig_film_radfrac = psf.viz.PlotFilmRadiativeFraction(liquid_film_results)
fig_film_boundary = psf.viz.PlotFilmBoundaryLayer(gaseous_film_results)
fig_film_gas_rad = psf.viz.PlotGaseousFilmRadiation(gaseous_film_results)
fig_film_gas_gap = psf.viz.PlotGaseousFilmTemperatureGap(gaseous_film_results)
#fig_film_wall_htc = psf.viz.PlotFilmEffectiveWallHTC(liquid_film_results, gaseous_film_results)

tab_film.add_figure(fig_film_wall_temperature)
tab_film.add_figure(fig_film_dryout)
tab_film.add_figure(fig_film_heat_flux)
tab_film.add_figure(fig_film_evaporation_rate)
tab_film.add_figure(fig_film_hconv)
tab_film.add_figure(fig_film_radfrac)
tab_film.add_figure(fig_film_boundary)
tab_film.add_figure(fig_film_gas_rad)
tab_film.add_figure(fig_film_gas_gap)
#tab_film.add_figure(fig_film_wall_htc)

tab_cooling_data = report.add_tab("Cooling Data")
fig_wall_temperature = psf.viz.PlotWallTemperature(cooling_data, plot_hot=True, plot_coolant_wall=True)
fig_coolant_temperature = psf.viz.PlotCoolantTemperature(cooling_data)
fig_coolant_pressure = psf.viz.PlotCoolantPressure(cooling_data)
fig_heat_flux = psf.viz.PlotHeatFlux(cooling_data)
fig_coolant_velocity = psf.viz.PlotVelocity(cooling_data)
fig_heat_transfer_coefficient = psf.viz.PlotHeatTransferCoefficient(cooling_data)

tab_cooling_data.add_figure(fig_wall_temperature)
tab_cooling_data.add_figure(fig_coolant_temperature)
tab_cooling_data.add_figure(fig_coolant_pressure)
tab_cooling_data.add_figure(fig_heat_flux)
tab_cooling_data.add_figure(fig_coolant_velocity)
tab_cooling_data.add_figure(fig_heat_transfer_coefficient)

# After all the content has been added, the report can be written to file
out_path = os.path.join(script_dir, "film_report.html")
report.save_html(out_path)
print(f"Report saved to {out_path}")