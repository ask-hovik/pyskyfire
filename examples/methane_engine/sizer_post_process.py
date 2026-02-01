import os
import pyskyfire as psf

script_dir = os.path.dirname(os.path.abspath(__file__)) 
in_path = os.path.join(script_dir, "sizer_results.pkl")

res = psf.common.Results.load(in_path)

net = res.net
stations = res.stations
signals = res.signals
residuals = res.residuals
params = res.input_params
cooling_data = [res.block_results["regen_throat_pass"], res.block_results["ox_regen"], ]
thrust_chamber = params["thrust_chamber"]



# ========
# Plotting
# ========

script_dir = os.path.dirname(os.path.abspath(__file__)) # Find current folder

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
#stl_path = os.path.join(script_dir, "engine_channels.stl") # save path and export filetype (all gmsh output filetypes available, including step)
#psf.viz.make_engine_gmsh(thrust_chamber, filename=stl_path) # Generates an STL of the engine and saves it in the script folder
#fig_cooling_channel_stl = psf.viz.EmbedSTL(stl_path) # Makes an embedable object out of an stl
fig_engine_contour = psf.viz.PlotContour(thrust_chamber.contour) # Plots the engine contour
#tab_overview.add_figure(fig_cooling_channel_stl)
tab_overview.add_figure(fig_engine_contour)

# Cooling Data 
tab_cooling_data = report.add_tab("Cooling Data")
fig_wall_temperature = psf.viz.PlotWallTemperature(*cooling_data, plot_hot=True, plot_coolant_wall=True, plot_interfaces=True) # Plot both hot and cold side wall temperature
fig_coolant_temperature = psf.viz.PlotCoolantTemperature(*cooling_data) # plots coolant temperature through the engine
fig_coolant_pressure = psf.viz.PlotCoolantPressure(*cooling_data) # plots coolant pressure through the engine
fig_heat_flux = psf.viz.PlotHeatFlux(*cooling_data) # plots heat flux through wall
fig_coolant_velocity = psf.viz.PlotVelocity(*cooling_data)
tab_cooling_data.add_figure(fig_wall_temperature)
tab_cooling_data.add_figure(fig_coolant_temperature)
tab_cooling_data.add_figure(fig_coolant_pressure)
tab_cooling_data.add_figure(fig_heat_flux)
tab_cooling_data.add_figure(fig_coolant_velocity)

# -------------- Engine Sizer Results ----------------
tab_sizer = report.add_tab("Engine Sizer")
fig_residual_history = psf.viz.PlotResidualHistory(residuals)
tab_sizer.add_figure(fig_residual_history, caption="Residual development over iterations")

# station list now available both from reference data and model
fu_stations = [
    "fu_engine_in", "fu_pump_in", "fu_pump_out",
    "fu_regen_in", "fu_regen_out", "fu_turbine_in",
    "fu_turbine_out", "fu_injector_plenum", "fu_chamber_in"
]

ox_stations = [
    "ox_engine_in", "ox_pump_in", "ox_pump_out",
    "ox_regen_in", "ox_regen_out", "ox_turbine_in",
    "ox_turbine_out", "ox_injector_plenum", "ox_chamber_in"
]

fig_fuel_side_pressure = psf.viz.PlotStationProperty(station_dicts=stations,
                            station_list=fu_stations,
                            property_name="p",
                            #labels=["Pyskyfire", "NASA, Binder"]
                            )

fig_fuel_side_temperature = psf.viz.PlotStationProperty(station_dicts=stations,
                            station_list=fu_stations,
                            property_name="T",
                            #labels=["Pyskyfire", "NASA, Binder"]
                            )

fig_fuel_side_mdot = psf.viz.PlotStationProperty(station_dicts=stations,
                            station_list=fu_stations,
                            property_name="mdot",
                            #labels=["Pyskyfire", "NASA, Binder"]
                            )

fig_ox_side_pressure = psf.viz.PlotStationProperty(station_dicts=stations,
                            station_list=ox_stations,
                            property_name="p",
                            #labels=["Pyskyfire", "NASA, Binder"]
                            )

fig_ox_side_temperature = psf.viz.PlotStationProperty(station_dicts=stations,
                            station_list=ox_stations,
                            property_name="T",
                            #labels=["Pyskyfire", "NASA, Binder"]
                            )

fig_ox_side_mdot = psf.viz.PlotStationProperty(station_dicts=stations,
                            station_list=ox_stations,
                            property_name="mdot",
                            #labels=["Pyskyfire", "NASA, Binder"]
                            )

tab_sizer.add_figure(fig_fuel_side_pressure, caption="The discrepancy here is a sort of compounding effect of the pressure drop rhtough the cooling channels being underestimated. When pressure drop rises over cooling channels, higher pressure rise is needed across pump, which again means higher pressure drop across turbine. I decided to leave it in because I do not actually know exactly why the pressure drop in the channels is lower than expected." )
tab_sizer.add_figure(fig_fuel_side_temperature, caption="Unknown why the reference data has a lower temperature rise than Pyskyfire. If you look at the temperature rise comparison in the regenerative cooling results they are not this excessive. ")
tab_sizer.add_figure(fig_fuel_side_mdot, caption="The discrepancies all over the board is exacerbated by the fact that I am not correcting for ideal vs delivered Isp. This lowers mass flow which has an effect on the whole system. When better Isp estimates are made in the future, these predictions will improve. ")
tab_sizer.add_figure(fig_ox_side_pressure)
tab_sizer.add_figure(fig_ox_side_temperature)
tab_sizer.add_figure(fig_ox_side_mdot)

fig_fuel_side_PT_diagram = psf.viz.PlotPTDiagram(station_dicts=[stations],
                                                 station_list=fu_stations,
                                                 fluid_name="methane",
                                                 title='Fuel-side P-T path',
                                                 scale="linear") #TODO: annotations in log mode don't work properly

fig_ox_side_PT_diagram = psf.viz.PlotPTDiagram(station_dicts=[stations],
                                                 station_list=ox_stations,
                                                 fluid_name="oxygen",
                                                 title='Ox-side P-T path',
                                                 scale="linear") 

tab_sizer.add_figure(fig_fuel_side_PT_diagram)
tab_sizer.add_figure(fig_ox_side_PT_diagram)

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
chamber_gradient = psf.viz.PlotTemperatureProfile(cooling_data[0], thrust_chamber, 0, -0.1) #this functionality is still a bit experimental
throat_gradient = psf.viz.PlotTemperatureProfile(cooling_data[0], thrust_chamber, 0, 0)
exit_gradient = psf.viz.PlotTemperatureProfile(cooling_data[0], thrust_chamber, 0, 0.05)
tab_thermal_gradient.add_figure(chamber_gradient)
tab_thermal_gradient.add_figure(throat_gradient)
tab_thermal_gradient.add_figure(exit_gradient)

# After all the content has been added, the report can be written to file
out_path = os.path.join(script_dir, "methane_engine_report.html")
report.save_html(out_path)
print(f"Report saved to {out_path}")