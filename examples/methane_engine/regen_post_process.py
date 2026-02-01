import os
import pyskyfire as psf

script_dir = os.path.dirname(os.path.abspath(__file__)) 
in_path = os.path.join(script_dir, "regen_results.pkl")
out_path = os.path.join(script_dir, "regen_report.html")

res = psf.common.Results.load(in_path)
params = res.params
thrust_chamber = res.thrust_chamber
cooling_data_a = res.cooling_data[0]
cooling_data_b = res.cooling_data[1]
cooling_data_c = res.cooling_data[2]
cooling_data_d = res.cooling_data[3]



# ========
# Plotting
# ========

# Generate report
print(f"Started generating report")
report = psf.viz.Report("Methane Engine")

# Engine Parameters
tab_params = report.add_tab("Parameters")
optimal_values = thrust_chamber.combustion_transport.optimum # A dictionary with the calculated optimums
tab_params.add_table(params, caption="Input Parameters", key_title="Parameter", value_title="Value", precision = 3) # input parameters we started with at the top of this file
tab_params.add_table(optimal_values, caption="Optimal Values", key_title="Parameters", value_title="Value", precision = 3)

# Engine Overview
tab_overview = report.add_tab("Engine Overview")

# 3D engine visualisation
plotter, viewer_src = psf.viz.make_engine_3d(thrust_chamber, show=False)
tab_overview.add_iframe(
    viewer_src,
    caption="Engine 3D",
) # deferred iframe lets you load the heavy engine model more intelligently

# Engine 2D contour
fig_engine_contour = psf.viz.PlotContour(thrust_chamber.contour) # Plots the engine contour
tab_overview.add_figure(fig_engine_contour)

# Cooling Data 
tab_cooling_data = report.add_tab("Cooling Data")
fig_wall_temperature = psf.viz.PlotWallTemperature(*res.cooling_data, plot_hot=True, plot_coolant_wall=True, plot_interfaces=True) # Plot both hot and cold side wall temperature
fig_coolant_temperature = psf.viz.PlotCoolantTemperature(*res.cooling_data) # plots coolant temperature through the engine
fig_coolant_pressure = psf.viz.PlotCoolantPressure(*res.cooling_data) # plots coolant pressure through the engine
fig_heat_flux = psf.viz.PlotHeatFlux(*res.cooling_data) # plots heat flux through wall
fig_coolant_velocity = psf.viz.PlotVelocity(*res.cooling_data)
tab_cooling_data.add_figure(fig_wall_temperature)
tab_cooling_data.add_figure(fig_coolant_temperature)
tab_cooling_data.add_figure(fig_coolant_pressure)
tab_cooling_data.add_figure(fig_heat_flux)
tab_cooling_data.add_figure(fig_coolant_velocity)

# Thrust Chamber Properties
tab_thrust_chamber_properties = report.add_tab("Thrust Chamber Properties")
fig_coolant_area = psf.viz.PlotCoolantArea(thrust_chamber)
fig_hydraulic_diameter = psf.viz.PlotHydraulicDiameter(thrust_chamber)
fig_radius_of_curvature = psf.viz.PlotRadiusOfCurvature(thrust_chamber)
fig_dAdx_thermal_hot_gas = psf.viz.PlotdAdxThermalHotGas(thrust_chamber)
fig_dAdx_thermal_coolant = psf.viz.PlotdAdxThermalCoolant(thrust_chamber)
fig_dAdx_coolant_area = psf.viz.PlotdAdxCoolantArea(thrust_chamber)
tab_thrust_chamber_properties.add_figure(fig_coolant_area)
tab_thrust_chamber_properties.add_figure(fig_hydraulic_diameter)
tab_thrust_chamber_properties.add_figure(fig_radius_of_curvature, caption="Currently some issues with the radius of curvature computation")
tab_thrust_chamber_properties.add_figure(fig_dAdx_thermal_hot_gas)
tab_thrust_chamber_properties.add_figure(fig_dAdx_thermal_coolant)
tab_thrust_chamber_properties.add_figure(fig_dAdx_coolant_area)

# Combustion
tab_combustion_reaction = report.add_tab("Combustion")
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
tab_combustion_reaction.add_figure(M)
tab_combustion_reaction.add_figure(gamma)
tab_combustion_reaction.add_figure(T)
tab_combustion_reaction.add_figure(p)
tab_combustion_reaction.add_figure(h)
tab_combustion_reaction.add_figure(cp)
tab_combustion_reaction.add_figure(k)
tab_combustion_reaction.add_figure(mu)
tab_combustion_reaction.add_figure(Pr)
tab_combustion_reaction.add_figure(rho)
tab_combustion_reaction.add_figure(a)

# Thermal Gradient
tab_thermal_gradient = report.add_tab("Thermal Gradient") 
chamber_gradient = psf.viz.PlotTemperatureProfile(cooling_data_c, thrust_chamber, 0, -0.1) #this functionality is still a bit experimental
throat_gradient = psf.viz.PlotTemperatureProfile(cooling_data_a, thrust_chamber, 0, 0)
exit_gradient = psf.viz.PlotTemperatureProfile(cooling_data_b, thrust_chamber, 0, 0.05)
tab_thermal_gradient.add_figure(chamber_gradient)
tab_thermal_gradient.add_figure(throat_gradient)
tab_thermal_gradient.add_figure(exit_gradient)

# After all the content has been added, the report can be written to file
report.save_html(out_path)
print(f"Report saved to {out_path}")


