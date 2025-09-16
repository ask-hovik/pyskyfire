# ============================================================
# Welcome to the minimal example showing how to use Pyskyfire!
# ============================================================

# This example will guide you through the creation of an engine, 
# and the running of a regenerative cooling simulation, before
# showing you how to plot the results using Pyskyfire's built-in
# plotting functionality. 

# The example is meant to showcase the capabilities of the Pyskyfire
# library in the simplest way possible, not being a good engine 
# design. Please follow proper design principles when designing
# your engine. 

import os
import pyskyfire as psf

# Input parameters for your engine
params = dict(
    # Chamber parameters
    p_c    = 20e5, # 20 bar chamber pressure
    F      = 5e3,  # 5kN thrust
    eps    = 4,    # Nozzle area ratio
    L_star = 1.2,  # Ratio of chamber volume to throat area
    MR     = 2.8,  # Mixture ratio of the combustion reaction
    AR_c   = 1.8,  # Chamber aspect ratio. Engines typically have a value between 1 and 2

    # Propellants for the two underlying codes, NASA CEA and coolprop, has to be loaded
    cea_fu = psf.common.Fluid(type="fuel", propellants=["C2H5OH"], fractions=[1.0]),
    cea_ox = psf.common.Fluid(type="oxidizer", propellants=["N2O"], fractions=[1.0]),
    coolprop_fu = "ethanol",
    T_coolant_in = 298.15, # The temperature of the fuel as it enters the cooling channels
    p_coolant_in = 23e5, # The pressure of the fuel as it enters the cooling channels. 

    # Cooling system properties
    material = psf.common.solids.StainlessSteel304, # There is built in support for a few different materials
    wall_thickness = 0.5e-3, # Half a millimeter wall thickness between the coolant and the combustion reaction
    n_channels = 80, # Number of cooling channels distributed around the engine
    blockage_ratio = 0.2, # This is a fraction that describes how many percent of the chamber cross section that are used up by ribs (20% in this case)
    roughness_height = 10e-6 # This is the average roughness height of the surface on the inside of the cooling channel. Here 10 micrometers
)  

# Any sim requires an description of the combustion reaction. The code behind this in Pyskyfire is NASA CEA.
# This object will also be able to tell you something about the optimal values of your engine
aerothermodynamics = psf.skycea.Aerothermodynamics.from_F_eps_Lstar(fu=params["cea_fu"], 
                                                                    ox=params["cea_ox"], 
                                                                    MR=params["MR"], 
                                                                    p_c=params["p_c"], 
                                                                    F=params["F"], 
                                                                    eps=params["eps"], 
                                                                    L_star=params["L_star"], 
                                                                    ) 

# The aerothermodynamics object has calculated the optimal chamber volume based on your input:
params["V_c"] = aerothermodynamics.V_c # Add to params dictionary for convencience
params["r_t"] = aerothermodynamics.r_t

# At the same time we might as well create an object to describe the coolant properties
# In pyskyfire the coolant properties code is based on CoolProp. In this example, the fuel is used as 
# the cooling fluid. 
coolant_transport = psf.skycea.CoolantTransport(params["coolprop_fu"])

# Helper function to generate a suitable contour. There are a lot more optional parameters to shape the chamber contour to your liking
xs, rs = psf.regen.contour.get_contour(V_c=params["V_c"],
                                       AR_c=params["AR_c"],
                                       r_t=params["r_t"],
                                       area_ratio=params["eps"], 
                                       nozzle = "conical", 
                                       )

# The raw contour data is encapsulated in a Contour object, which calculates a lot of the properties of the engine contour created
contour = psf.regen.Contour(xs, rs, name = "Minimal Contour")

# A wall object has to be created to represent the barrier between the coolant and the combustion reaction
wall = psf.regen.Wall(material  = params["material"],
                      thickness = params["wall_thickness"]) 

# Support for multiple layered walls is implemented, therefore a wall group object has to be created to calculate the properties of the 
# entire wall stack. For this simulationi a single wall with no coatings or multiple materials is used. 
wall_group = psf.regen.WallGroup(walls=[wall])

# The cooling channels can have multiple different cross sectional shapes. Examples of this are rounded cooling 
# channels like the RL10 engine, or squared/rectangular cooling channels like on the Vulcain engine. 
# For this simulation, a squared cooling channel will be used.
cross_section = psf.regen.CrossSectionSquared()

# A channel height function has to be supplied, to let the geometry creation code know how "tall" 
# the cooling channel is at any point. For this sim we will just create a constant function
def channel_height_function(x):
    return 2e-3 # constant 2mm tall cooling channel

# The SurfacePlacement class is a class that defines how the cooling channels are wrapped around the thrust chamber. 
# The class is in need of a bit of an update, SurfacePlacement in particular describes a vertical set of cooling channels. 
# In principle, cooling channels could be placed in a helical configuration, or even a sine pattern or anything you can imagine really. 
# Channels could even be placed inside the chamber or outside the chamber for various reasons. This is encapsulated in this class. 
# For this simulation, we are using a simple vertical channel configuration. The number of channel positions just means how many channels
# are distributed around the engine. 
surface_placement = psf.regen.SurfacePlacement(n_channel_positions=params["n_channels"])

# Now a cooling channel circuit has to be made. The cooling circuit object takes responsibility to supply information about
# the geometry and contents of the cooling channel. Pyskyfire supports placing multiple cooling circuits around a single engine.
# This means different circuits can cover different parts of the engine, and even interlace with each other. 
# For this example, a single cooling circuit running from the bottom of the nozzle to the top of the chamber is created. 
cooling_circuit = psf.regen.CoolingCircuit(name="Cooling Pass", 
                                     contour=contour, # The contour we created earlier
                                     coolant_transport=coolant_transport, # Coolant properties we establised earlier
                                     cross_section=cross_section, # The type of cross section we are using (rounded/squared)
                                     span = [1.0, -1.0], # This is over what span of the contour this particular circuit is placed. [1.0, -1.0] means the cooling circuit goes from the bottom of the nozzle to the top of the chamber.
                                     placement=surface_placement, # Surface placement is a class that describe how the thrust chamber is wrapped in channels
                                     channel_height=channel_height_function, # Channel height as a function of x
                                     blockage_ratio=params["blockage_ratio"]) # How many percent of the circumference of the engine is taken up by ribs, only needed for certain cross sections, like the squared cross section. 

# In the case multiple cooling circuits are created, they are grouped together in this group object. 
cooling_circuit_group = psf.regen.CoolingCircuitGroup(circuit_list=[cooling_circuit])

# The thrust chamber object is an umbrella class that takes in and consolidates all the objects we have created into a complete thrust chamber
# It calculates values that are only completely defined when all pieces of the puzzle are present, and then exposes those properties to you
thrust_chamber = psf.regen.ThrustChamber(contour=contour, 
                                         wall_group=wall_group,
                                         combustion_transport=aerothermodynamics,  
                                         cooling_circuit_group=cooling_circuit_group,
                                         roughness=params["roughness_height"])


# We now have a thrust chamber. At this point it might be useful to plot the contour and export a 3d-model of the 
# thrust chamber. You can experiment with that as you design your own engine. Plotting will however be a major 
# theme later in this example. The most interesting thing that can be done with this thrust chamber is to simulate
# the regenerative cooling of this chamber at steady state during a firing. 

# If we imagine that all the fuel that are flowing through the cooling channels end up in the chamber
# then the mass flow is the same between the two, and we can pull the fuel mass flow from the combustion object. 
mdot_fu = aerothermodynamics.mdot_fu

# These are the inlet conditions of the regenerative cooling loop. We just chose these in the input parameters, but in more 
# sophisticated runs this can be part of an iterative loop such that this value is solved for. 
boundary_conditions = psf.regen.BoundaryConditions(T_coolant_in = params["T_coolant_in"], p_coolant_in = params["p_coolant_in"], mdot_coolant = mdot_fu)

# This is the function running the simulation. 
cooling_data = psf.regen.steady_heating_analysis(thrust_chamber, 
                                                 n_nodes = 100, # Number of simulation nodes
                                                 circuit_index=0, # If there were more circuits, this index would indicate which one is being simulated
                                                 boundary_conditions=boundary_conditions, 
                                                 solver="newton", # The newton solver is the only one implemented at the moment
                                                 output=True) # Prints simulation progres


# At this point we have obtained the regenerative cooling result for the engine we just designed. We have created
# a lot of data. You could of course take this raw data and process it yourself, that's great! But Pyskyfire has a 
# built in plotting library that helps you plot and tabulate the data you are most likely interested in. 
# The built in plotting consists of a number of classes that plot the data with sensible presets. 
# The instances of those classes can then be added to a report, which is an HTML file that is organized
# easily readable, and portable, so can be shared easily. 

# ========
# Plotting
# ========

script_dir = os.path.dirname(os.path.abspath(__file__)) # Find current folder

# Generate report
print(f"Started generating report")
report = psf.viz.Report("Minimal Engine")

# Engine Parameters
tab_params = report.add_tab("Parameters")
optimal_values = thrust_chamber.combustion_transport.optimum # A dictionary with the calculated optimums
tab_params.add_table(params, caption="Input Parameters", key_title="Parameter", value_title="Value", precision = 3) # input parameters we started with at the top of this file
tab_params.add_table(optimal_values, caption="Optimal Values", key_title="Parameters", value_title="Value", precision = 3)

# Engine Overview
tab_overview = report.add_tab("Engine Overview")
stl_path = os.path.join(script_dir, "engine_channels.stl") # save path and export filetype (all gmsh output filetypes available, including step)
psf.viz.make_engine_gmsh(thrust_chamber, filename=stl_path) # Generates an STL of the engine and saves it in the script folder
fig_cooling_channel_stl = psf.viz.EmbedSTL(stl_path) # Makes an embedable object out of an stl
fig_engine_contour = psf.viz.PlotContour(thrust_chamber.contour) # Plots the engine contour
tab_overview.add_figure(fig_cooling_channel_stl)
tab_overview.add_figure(fig_engine_contour)

# Cooling Data 
tab_cooling_data = report.add_tab("Cooling Data")
fig_wall_temperature = psf.viz.PlotWallTemperature(cooling_data, plot_hot=True, plot_coolant_wall=True) # Plot both hot and cold side wall temperature
fig_coolant_temperature = psf.viz.PlotCoolantTemperature(cooling_data) # plots coolant temperature through the engine
fig_coolant_pressure = psf.viz.PlotCoolantPressure(cooling_data) # plots coolant pressure through the engine
fig_heat_flux = psf.viz.PlotHeatFlux(cooling_data) # plots heat flux through wall
fig_coolant_velocity = psf.viz.PlotVelocity(cooling_data)
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
chamber_gradient = psf.viz.PlotTemperatureProfile(cooling_data, thrust_chamber, 0, -0.1) #this functionality is still a bit experimental
throat_gradient = psf.viz.PlotTemperatureProfile(cooling_data, thrust_chamber, 0, 0)
exit_gradient = psf.viz.PlotTemperatureProfile(cooling_data, thrust_chamber, 0, 0.05)
tab_thermal_gradient.add_figure(chamber_gradient)
tab_thermal_gradient.add_figure(throat_gradient)
tab_thermal_gradient.add_figure(exit_gradient)

# After all the content has been added, the report can be written to file
out_path = os.path.join(script_dir, "minimal_report.html")
report.save_html(out_path)
print(f"Report saved to {out_path}")