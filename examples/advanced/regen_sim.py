import os
import pyskyfire as psf


# =============
# Engine inputs
# =============

params = dict(
    # All values from Binder et. al. unless stated otherwise
    p_c = 100e5,
    p_e = 0.8e5,
    MR = 2.8,
    AR_c = 1.8,
    F = 300e3,

    # Fuel/oxidizer parameters
    cea_fu = psf.common.Fluid(type="fuel", propellants=["CH4"], fractions=[1.0]),
    cea_ox = psf.common.Fluid(type="oxidizer", propellants=["O2"], fractions=[1.0]),
    coolant_fu =  psf.common.Fluid(type="fuel", propellants=["methane"], fractions=[1.0]),
    coolant_ox = psf.common.Fluid(type="oxidizer", propellants=["oxygen"], fractions=[1.0]),
    T_gas_fu_in = 250,  # Adjusted by a few kelvin to fall within NASA CEA tables
    T_gas_ox_in = 250,  # Also adjusted

    T_coolant_fu_in = 111,
    p_coolant_fu_in = 150e5, 

    T_coolant_ox_in = 93,
    p_coolant_ox_in = 150e5,

    # Chamber-nozzle parameters, contour from Binder et. al.
    theta_conv = 35,
    R_1f = 1.5,
    R_2f = 2,
    R_3f = 0.5,
    length_fraction = 0.8, # 80% nozzle
    L_star = 1.1,  # I can't remember where I read this (not Binder)

    # Cooling channel
    copper_roughness_height = 0.0030e-3,  # found in Binder et al 4.6e-5 inches. Conventional SS tubing is often cited as 0.0015e-3
    inconel_roughtness_height = 0.0015e-3,

    # thicknesses
    oxygen_copper_thickness = 0.6e-3,
    fuel_copper_thickness = 0.6e-3,
    barrier_thickness = 0.1e-3,
    inconel_thickness = 0.4e-3,

    # Materials
    copper = psf.common.solids.GRCop42,
    zirconia = psf.common.solids.ZirconiumOxide,
    inconel = psf.common.solids.Inconel625

)




# =====
# Setup
# =====

aerothermodynamics = psf.skycea.Aerothermodynamics.from_F_pe_Lstar(fu=params["cea_fu"], 
                                                                    ox=params["cea_ox"], 
                                                                    T_fu_in=params["T_gas_fu_in"], 
                                                                    T_ox_in=params["T_gas_ox_in"], 
                                                                    MR=params["MR"], 
                                                                    p_c=params["p_c"], 
                                                                    F=params["F"], 
                                                                    p_e=params["p_e"], 
                                                                    L_star=params["L_star"], 
                                                                    p_amb=1e5, # ca atmospheric pressure, optional
                                                                    npts=15 # good compromise between speed and precision
                                                                    ) 

V_c = aerothermodynamics.V_c 
r_t = aerothermodynamics.r_t
eps = aerothermodynamics.eps


# Generate the contour coordinates
xs, rs = psf.regen.contour.get_contour(V_c = V_c, 
                                       r_t = r_t, 
                                       area_ratio = eps, 
                                       AR_c = params["AR_c"], 
                                       theta_conv = params["theta_conv"],
                                       nozzle = "rao", 
                                       R_1f=params["R_1f"],
                                       R_2f=params["R_2f"],
                                       R_3f=params["R_3f"],
                                       length_fraction = params["length_fraction"]
                                       )

# Make a contour object using the coordinates
contour = psf.regen.Contour(xs, rs, name = "Methane Engine")
#fig = psf.viz.PlotContour(contour)
#fig.show()
#input()

oxygen_copper_wall = psf.regen.Wall(material = params["copper"], thickness = params["oxygen_copper_thickness"]) 
fuel_copper_wall = psf.regen.Wall(material = params["copper"], thickness = params["fuel_copper_thickness"])
barrier_wall = psf.regen.Wall(material = params["zirconia"], thickness = params["barrier_thickness"])

inconel_wall = psf.regen.Wall(material=params["inconel"], thickness=params["inconel_thickness"])
#wall_group = psf.regen.WallGroup(walls=[barrier_wall, copper_wall])

channel_height_fn = psf.regen.make_channel_height_fn(
    contour=contour, 
    region_fractions=[-1.0, 0.3, 1.0], 
    flat_heights= [0.00032, 0.00134], 
    pinch_factors= [0.8, -5.0], 
    transition_widths=[0.1]
) # total coolant volume should be ca 0.015831543m3

def simple_height(x):
    return 1.5e-3

def internal_width(x):
    return 20e-3

def internal_height(x):
    return 20e-3


cross_section = psf.regen.CrossSectionRounded()
fuel_transport = psf.skycea.CoolantTransport(params["coolant_fu"])
ox_transport = psf.skycea.CoolantTransport(params["coolant_ox"])
channel_placement = psf.regen.SurfacePlacement(n_channel_positions=120)

fu_copper_pass = psf.regen.CoolingCircuit(name="Fuel Copper Pass", 
                                     contour=contour, 
                                     coolant_transport=fuel_transport, 
                                     cross_section=cross_section, 
                                     span = [-0.3, 0.29], 
                                     placement=channel_placement,
                                     walls = [barrier_wall, fuel_copper_wall],
                                     roughness = params["copper_roughness_height"],
                                     channel_height=channel_height_fn)

fu_inco_pass_cocurrent = psf.regen.CoolingCircuit(name="Fuel Inconel Pass Cocurrent", 
                                     contour=contour, 
                                     coolant_transport=fuel_transport, 
                                     cross_section=cross_section, 
                                     span = [0.3, 1.0], 
                                     placement=channel_placement,
                                     walls = [inconel_wall],
                                     roughness = params["inconel_roughtness_height"],
                                     channel_height=channel_height_fn)

fu_inco_pass_countercurrent = psf.regen.CoolingCircuit(name="Fuel Inconel Pass Countercurrent", 
                                          contour=contour,
                                          coolant_transport=fuel_transport,
                                          cross_section=cross_section,
                                          span = [1.0, 0.3],
                                          placement=channel_placement,
                                          walls = [inconel_wall],
                                          roughness = params["inconel_roughtness_height"],
                                          channel_height=channel_height_fn)

ox_copper_pass = psf.regen.CoolingCircuit(name= "Oxidizer Copper Pass",
                                     contour=contour, 
                                     coolant_transport=ox_transport, 
                                     cross_section=cross_section, 
                                     span = [-0.31, -1.0], 
                                     placement=psf.regen.SurfacePlacement(n_channel_positions=240),
                                     walls = [barrier_wall, oxygen_copper_wall], 
                                     roughness = params["copper_roughness_height"],
                                     channel_height=simple_height)

"""ox_internal_pass = psf.regen.CoolingCircuit(name = "Oxidizer Internal Pass", 
                                            contour=contour,
                                            coolant_transport = ox_transport, 
                                            cross_section = cross_section, 
                                            span = [-0.31, -1.0], 
                                            placement = psf.regen.InternalPlacement(channel_width=internal_width, n_channel_positions=16, n_channels_per_leaf=1),
                                            walls = [barrier_wall, oxygen_copper_wall],
                                            roughness = params["copper_roughness_height"],
                                            channel_height=internal_height
                                            )"""




thrust_chamber = psf.regen.ThrustChamber(contour=contour, 
                                         combustion_transport=aerothermodynamics,  
                                         cooling_circuits=[fu_copper_pass, fu_inco_pass_cocurrent, fu_inco_pass_countercurrent, ox_copper_pass],
                                         h_gas_corr=1.0, # No correction applied
                                         h_cold_corr=1.0, # No correction applied
                                         ) 
#psf.viz.make_engine_3d(thrust_chamber)
#input()

script_dir = os.path.dirname(os.path.abspath(__file__)) 
stl_path = os.path.join(script_dir, "internal_channels.stl")
#psf.viz.make_engine_gmsh(thrust_chamber, filename=stl_path, display_channels=8) # Generates an STL of the engine and saves it in the script folder
#input()

#eng_3d = psf.viz.make_engine_3d_pyvista(thrust_chamber)
#eng_3d.show()
#del eng_3d

# Fuel Cooling Pass
mdot_fu = aerothermodynamics.mdot_fu
boundary_conditions_a = psf.regen.BoundaryConditions(T_coolant_in = params["T_coolant_fu_in"], 
                                                     p_coolant_in = params["p_coolant_fu_in"], 
                                                     mdot_coolant = mdot_fu)

cooling_data_a = psf.regen.steady_heating_analysis(thrust_chamber, 
                                                   n_nodes = 25, 
                                                   circuit_index=0, 
                                                   boundary_conditions=boundary_conditions_a, 
                                                   solver="newton", 
                                                   output=True)

T_out = cooling_data_a['T_stagnation'][-1]
p_out = cooling_data_a['p_stagnation'][-1]

boundary_conditions_b = psf.regen.BoundaryConditions(T_coolant_in = T_out, 
                                                     p_coolant_in = p_out, 
                                                     mdot_coolant = mdot_fu)

cooling_data_b = psf.regen.steady_heating_analysis(thrust_chamber, 
                                                   n_nodes = 30, 
                                                   circuit_index=1, 
                                                   boundary_conditions=boundary_conditions_b, 
                                                   solver="newton", 
                                                   output=True)

T_out = cooling_data_b['T_stagnation'][-1]
p_out = cooling_data_b['p_stagnation'][-1]

boundary_conditions_c = psf.regen.BoundaryConditions(T_coolant_in = T_out, 
                                                     p_coolant_in = p_out, 
                                                     mdot_coolant = mdot_fu)

cooling_data_c = psf.regen.steady_heating_analysis(thrust_chamber, 
                                                   n_nodes = 30, 
                                                   circuit_index=2, 
                                                   boundary_conditions=boundary_conditions_c)

# Ox pass
mdot_ox_1 = aerothermodynamics.mdot_ox*0.5
boundary_conditions_d = psf.regen.BoundaryConditions(T_coolant_in=params["T_coolant_ox_in"],
                                                     p_coolant_in=params["p_coolant_ox_in"],
                                                     mdot_coolant=mdot_ox_1)

cooling_data_d = psf.regen.steady_heating_analysis(thrust_chamber, 
                                                   n_nodes=50,
                                                   circuit_index=3,
                                                   boundary_conditions=boundary_conditions_d,
                                                   solver="newton",
                                                   output=True)

"""mdot_ox_2 = aerothermodynamics.mdot_ox*0.1
boundary_conditions_e = psf.regen.BoundaryConditions(T_coolant_in=params["T_coolant_ox_in"],
                                                     p_coolant_in=params["p_coolant_ox_in"],
                                                     mdot_coolant=mdot_ox_2)

cooling_data_e = psf.regen.steady_heating_analysis(thrust_chamber, 
                                                   n_nodes=50,
                                                   circuit_index=4,
                                                   boundary_conditions=boundary_conditions_d,
                                                   solver="newton",
                                                   output=True)"""


# Save results
script_dir = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(script_dir, "regen_results.pkl")

res = psf.common.Results()
res.add(name="params", obj=params)
res.add(name="thrust_chamber", obj=thrust_chamber)
res.add(name="cooling_data", obj=[cooling_data_a, cooling_data_b, cooling_data_c, cooling_data_d])
res.save(out_path)