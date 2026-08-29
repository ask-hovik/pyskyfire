# ==================================================================
# Coupled film + regenerative cooling
#
# The film is injected part-way down the chamber rather than at the
# injector face, so the barrel upstream of it is cooled regeneratively
# only and the two regimes can be compared on one engine.
#
# The same chamber is run three ways:
#   1. film cooling on its own
#   2. regenerative cooling on its own
#   3. the two coupled, with the film setting the hot-side heat flux
#      wherever it covers the wall
# ==================================================================

import argparse
from pathlib import Path

import numpy as np

import pyskyfire as psf


parser = argparse.ArgumentParser()
parser.add_argument(
    "--output-dir",
    type=Path,
    help="Directory for generated HTML outputs.",
)
args = parser.parse_args()

output_dir = args.output_dir or Path(__file__).resolve().parent
output_dir.mkdir(parents=True, exist_ok=True)

# Single-figure HTML files, written next to the report so the documentation can
# embed one plot without carrying a second copy of the case definition.
standalone_dir = output_dir / "standalone"
standalone_dir.mkdir(parents=True, exist_ok=True)

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

    # Film cooling inputs.
    # x_fraction is signed: -1 = injector face, 0 = throat, +1 = nozzle exit.
    # -0.7 puts the film ring ~30% of the chamber length downstream of the
    # injector, leaving that first stretch regeneratively cooled only.
    film_x_fraction = -0.70,
    film_mdot = 0.02,
    film_injection_perimeter = 0.10,
    film_liquid_absorptivity = 0.2,
    film_turbulence_intensity = 0.02,
    film_mole_fraction_CO2 = 0.1,
    film_mole_fraction_H2O = 0.4,
)

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

wall = psf.regen.Wall(material=params["material"], thickness=params["wall_thickness"])
cross_section_squared = psf.regen.CrossSectionSquared(blockage_ratio=params["blockage_ratio"])

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

# The film is bled off the fuel side, so the channels only see what is left.
mdot_regen = aerothermodynamics.mdot_fu - params["film_mdot"]

boundary_conditions = psf.regen.BoundaryConditions(
    T_coolant_in=params["T_coolant_in"],
    p_coolant_in=params["p_coolant_in"],
    mdot_coolant=mdot_regen,
)

x_injection = thrust_chamber.film_cooling.x

# ------------------------------------------------------------------
# 1. Film cooling alone
# ------------------------------------------------------------------
liquid, gaseous = psf.regen.solve_film_cooling(thrust_chamber, boundary_conditions)

print("\n=== Film cooling only ===")
print(f"Chamber runs from x = {contour.xs[0]:.4f} m to the throat at x = {contour.x_t:.4f} m")
print(f"Film injected at   x = {x_injection:.4f} m")
if liquid.x_dryout is None:
    print("Liquid film survives to the end of the domain")
else:
    print(f"Dryout at x = {liquid.x_dryout:.4f} m, {len(gaseous.x)} gaseous-film nodes")

# ------------------------------------------------------------------
# 2. Regenerative cooling alone
# ------------------------------------------------------------------
regen_only = psf.regen.coupled_steady_heating_analysis(
    thrust_chamber,
    boundary_conditions,
    nodes=100,
    circuit_index=0,
    film=False,
    heat_curvature_correction=False,
)

# ------------------------------------------------------------------
# 3. Coupled
# ------------------------------------------------------------------
coupled = psf.regen.coupled_steady_heating_analysis(
    thrust_chamber,
    boundary_conditions,
    nodes=100,
    circuit_index=0,
    film=True,
    h_liquid_wall=psf.regen.DEFAULT_H_LIQUID_WALL,
    heat_curvature_correction=False,
)

# The regen plotters label each trace with ``result.circuit_name``, so naming
# the two runs apart is all it takes to overlay them on one set of axes.
regen_only.circuit_name = "Regen only"
coupled.circuit_name = "Regen + film"

# ------------------------------------------------------------------
# Comparison
# ------------------------------------------------------------------
# Upstream of the film ring the hot side is identical between the two runs; any
# difference there comes back through the coolant, which passes the film-cooled
# stretch first because the circuit flows nozzle-to-injector.
upstream = regen_only.x < x_injection

summary = {
    "Peak hot-wall temperature, regen only [K]": regen_only.T[:, -1].max(),
    "Peak hot-wall temperature, regen + film [K]": coupled.T[:, -1].max(),
    "Peak heat flux, regen only [MW/m^2]": regen_only.dQ_dA.max() / 1e6,
    "Peak heat flux, regen + film [MW/m^2]": coupled.dQ_dA.max() / 1e6,
    "Coolant outlet temperature, regen only [K]": regen_only.T_static[-1],
    "Coolant outlet temperature, regen + film [K]": coupled.T_static[-1],
    "Coolant pressure drop, regen only [bar]":
        (regen_only.p_stagnation[0] - regen_only.p_stagnation[-1]) / 1e5,
    "Coolant pressure drop, regen + film [bar]":
        (coupled.p_stagnation[0] - coupled.p_stagnation[-1]) / 1e5,
    "Wall cooling upstream of the film ring, max [K]":
        (regen_only.T[upstream, -1] - coupled.T[upstream, -1]).max(),
    "Film injection station [m]": x_injection,
    "Liquid film dryout station [m]":
        liquid.x_dryout if liquid.x_dryout is not None else float("nan"),
    "Film mass flow [kg/s]": params["film_mdot"],
    "Regenerative coolant mass flow [kg/s]": mdot_regen,
}

print("\n=== Regen only vs. coupled ===")
print(f"{'quantity':<28}{'regen only':>14}{'coupled':>14}")
rows = [
    ("peak hot-wall T [K]", regen_only.T[:, -1].max(), coupled.T[:, -1].max()),
    ("coolant outlet T [K]", regen_only.T_static[-1], coupled.T_static[-1]),
    ("peak heat flux [MW/m2]", regen_only.dQ_dA.max() / 1e6, coupled.dQ_dA.max() / 1e6),
    ("coolant dp [bar]", (regen_only.p_stagnation[0] - regen_only.p_stagnation[-1]) / 1e5,
     (coupled.p_stagnation[0] - coupled.p_stagnation[-1]) / 1e5),
]
for name, a, b in rows:
    print(f"{name:<28}{a:>14.3f}{b:>14.3f}")

regime_rows = {}
for regime in ("gas", "liquid_film", "gaseous_film"):
    n = int(np.sum(coupled.film_regime == regime))
    if n:
        x_reg = coupled.x[coupled.film_regime == regime]
        print(f"{regime:<14} {n:>4} nodes  x = [{x_reg.min():.4f}, {x_reg.max():.4f}] m")
        regime_rows[f"{regime} coverage [m]"] = f"{x_reg.min():.4f} to {x_reg.max():.4f}"
        regime_rows[f"{regime} nodes"] = n

# ==================================================================
# Report
# ==================================================================
print("\nStarted generating report")
report = psf.viz.Report("Coupled Film + Regenerative Cooling")

# --- Parameters ---------------------------------------------------
tab_params = report.add_tab("Parameters")
tab_params.add_table(params, caption="Input Parameters", key_title="Parameter",
                     value_title="Value", precision=3)
tab_params.add_table(thrust_chamber.combustion_transport.optimum, caption="Optimal Values",
                     key_title="Parameter", value_title="Value", precision=3)

# --- Summary ------------------------------------------------------
tab_summary = report.add_tab("Summary")
tab_summary.add_markdown(
    "The film ring sits downstream of the injector, so the upstream stretch of "
    "the barrel is cooled regeneratively only. Both runs use the same chamber, "
    "the same channel geometry and the same regenerative coolant flow, so the "
    "difference between them is the film alone."
)
tab_summary.add_table(summary, caption="Headline Results", key_title="Quantity",
                      value_title="Value", precision=4)
tab_summary.add_table(regime_rows, caption="Hot-Side Regime Coverage",
                      key_title="Regime", value_title="Value", precision=4)

# --- Engine overview ----------------------------------------------
tab_overview = report.add_tab("Engine Overview")
fig_contour = psf.viz.PlotContour(thrust_chamber.contour)
fig_contour.add.vline(x=x_injection, line_dash="dot", line_color="firebrick",
                      annotation_text="film injection", annotation_position="top")
tab_overview.add_figure(fig_contour, caption="Chamber contour with the film injection station.")

# --- Film cooling -------------------------------------------------
tab_film = report.add_tab("Film Cooling")
T_g_ref = thrust_chamber.combustion_transport.get_T(x_injection)
p_ref = thrust_chamber.combustion_transport.get_p(x_injection)
T_sat_ref = (coolant_transport.get_T_sat(params["p_coolant_in"])
             if hasattr(coolant_transport, "get_T_sat") else None)

tab_film.add_markdown(
    "Grisson film model, solved on its own. The liquid film absorbs the gas-side "
    "heat load as evaporation until it dries out, after which its vapour "
    "continues to shield the wall at a reduced recovery temperature."
)
for fig, caption in [
    (psf.viz.PlotFilmWallTemperature(liquid, gaseous, T_g_ref=T_g_ref, T_sat_ref=T_sat_ref),
     "Film and wall temperature through the liquid and gaseous regions."),
    (psf.viz.PlotFilmDryout(liquid), "Film flow rate per unit width down to dryout."),
    (psf.viz.PlotFilmHeatFlux(liquid), "Convective and radiative load on the liquid film."),
    (psf.viz.PlotFilmEvaporationRate(liquid), "Local evaporation rate."),
    (psf.viz.PlotFilmHeatTransferCoefficient(liquid, gaseous),
     "Film-side heat transfer coefficient. Note the gaseous branch reports an "
     "entrainment conductance, not a wall coefficient."),
    (psf.viz.PlotFilmRadiativeFraction(liquid), "Radiative share of the liquid-film heat load."),
    (psf.viz.PlotFilmBoundaryLayer(gaseous), "Boundary-layer mass flow after dryout."),
    (psf.viz.PlotGaseousFilmRadiation(gaseous), "Radiative flux through the gaseous film."),
    (psf.viz.PlotGaseousFilmTemperatureGap(gaseous),
     "Gap between the film-cooled adiabatic wall temperature and the wall."),
]:
    tab_film.add_figure(fig, caption=caption)

# --- Cooling data, with and without film --------------------------
tab_cooling = report.add_tab("Cooling Data")
tab_cooling.add_markdown(
    "Both runs overlaid. Note that the circuit flows nozzle-to-injector, so the "
    "coolant crosses the film-cooled stretch before it reaches the bare barrel "
    "upstream of the injection station. The hot-side boundary condition there is "
    "identical between the two runs, but the wall still comes out cooler with the "
    "film, because the coolant arrives having picked up less heat on the way. The "
    "curves therefore do not merge upstream of injection."
)

fig_wall = psf.viz.PlotWallTemperature(regen_only, coupled,
                                       plot_hot=True, plot_coolant_wall=False)
fig_wall.add.vline(x=x_injection, line_dash="dot", line_color="firebrick",
                   annotation_text="film injection", annotation_position="top left")
if liquid.x_dryout is not None:
    fig_wall.add.vline(x=liquid.x_dryout, line_dash="dot", line_color="darkorange",
                       annotation_text="dryout", annotation_position="top right")
fig_wall.save_html(standalone_dir / "film-regen-wall-temperature.html")
tab_cooling.add_figure(
    fig_wall,
    caption="Hot-side wall temperature with and without film cooling.",
)

fig_flux = psf.viz.PlotHeatFlux(regen_only, coupled)
fig_flux.add.vline(x=x_injection, line_dash="dot", line_color="firebrick",
                   annotation_text="film injection", annotation_position="top left")
fig_flux.save_html(standalone_dir / "film-regen-heat-flux.html")
tab_cooling.add_figure(fig_flux, caption="Wall heat flux with and without film cooling.")

for fig, caption in [
    (psf.viz.PlotCoolantTemperature(regen_only, coupled), "Bulk coolant temperature."),
    (psf.viz.PlotCoolantPressure(regen_only, coupled), "Coolant pressure."),
    (psf.viz.PlotHeatTransferCoefficient(regen_only, coupled),
     "Hot- and coolant-side heat transfer coefficients."),
    (psf.viz.PlotVelocity(regen_only, coupled), "Coolant velocity."),
]:
    tab_cooling.add_figure(fig, caption=caption)

out_path = output_dir / "coupled_film_report.html"
report.save_html(out_path)
print(f"Report saved to {out_path}")
