import os
import pyskyfire as psf


# =============
# Engine inputs
# =============

params = dict(
    # All values from Binder et. al. unless stated otherwise
    p_c = 32.7501e5,
    eps = 56.1,  # Without silver insert
    MR = 5.0,
    r_c = 0.123,  # From graph
    F = 73.4e3,

    # Fuel/oxidizer parameters
    cea_fu = psf.common.Fluid(type="fuel", propellants=["H2"], fractions=[1.0]),
    cea_ox = psf.common.Fluid(type="oxidizer", propellants=["O2(L)"], fractions=[1.0]),
    coolprop_fu = "hydrogen",
    coolprop_ox = "oxygen",
    T_gas_fu_in = 200,  # Adjusted by a few kelvin to fall within NASA CEA tables
    T_gas_ox_in = 100,  # Also adjusted
    T_coolant_in = 33,
    p_coolant_in = 69e5, 

    # Chamber-nozzle parameters, contour from Binder et. al.
    theta_conv = 25,
    R_1f = 1.5,
    R_2f = 3,
    R_3f = 0.5,
    length_fraction = 0.713,
    L_star = 0.95,  # I can't remember where I read this (not Binder)

    # Cooling channel
    roughness_height = 1.1684e-6,  # found in Binder et al 4.6e-5 inches. Conventional SS tubing is often cited as 0.0015e-3
    wall_thickness = 0.31e-3,
)


# =====
# Setup
# =====

aerothermodynamics = psf.skycea.Aerothermodynamics.from_F_eps_Lstar(fu=params["cea_fu"], 
                                                                    ox=params["cea_ox"], 
                                                                    T_fu_in=params["T_gas_fu_in"], 
                                                                    T_ox_in=params["T_gas_ox_in"], 
                                                                    MR=params["MR"], 
                                                                    p_c=params["p_c"], 
                                                                    F=params["F"], 
                                                                    eps=params["eps"], 
                                                                    L_star=params["L_star"], 
                                                                    p_amb=1e5, # ca atmospheric pressure, optional
                                                                    npts=15 # good compromise between speed and precision
                                                                    ) 

V_c = aerothermodynamics.V_c 
r_t = aerothermodynamics.r_t

# Generate the contour coordinates
xs, rs = psf.regen.contour.get_contour(V_c = V_c, 
                                       r_t = r_t, 
                                       area_ratio = params["eps"], 
                                       r_c = params["r_c"], 
                                       theta_conv = params["theta_conv"],
                                       nozzle = "rao", 
                                       R_1f=params["R_1f"],
                                       R_2f=params["R_2f"],
                                       R_3f=params["R_3f"],
                                       length_fraction = params["length_fraction"]
                                       )

# Make a contour object using the coordinates
contour = psf.regen.Contour(xs, rs, name = "RL10 contour")

wall = psf.regen.Wall(material = psf.common.solids.StainlessSteel304, thickness = params["wall_thickness"]) 
wall_group = psf.regen.WallGroup(walls=[wall])

channel_height_fn = psf.regen.make_channel_height_fn(
    contour=contour, 
    region_fractions=[-1.0, 0.25, 1.0], 
    flat_heights= [0.0032, 0.00134], 
    pinch_factors= [0.7, -5.0], 
    transition_widths=[0.1]
) # total coolant volume should be ca 0.015831543m3


cross_section = psf.regen.CrossSectionRounded()
LH2_transport = psf.skycea.CoolantTransport(params["coolprop_fu"])

half_pass = psf.regen.CoolingCircuit(name="Half Pass", 
                                     contour=contour, 
                                     coolant_transport=LH2_transport, 
                                     cross_section=cross_section, 
                                     span = [0.25, 1.0], 
                                     placement=psf.regen.SurfacePlacement(n_channel_positions=180), 
                                     channel_height=channel_height_fn)

full_pass = psf.regen.CoolingCircuit(name= "Full Pass",
                                     contour=contour, 
                                     coolant_transport=LH2_transport, 
                                     cross_section=cross_section, 
                                     span = [1.0, -1.0], 
                                     placement=psf.regen.SurfacePlacement(n_channel_positions=180), 
                                     channel_height=channel_height_fn)

cooling_circuit_group = psf.regen.CoolingCircuitGroup(circuit_list=[half_pass, full_pass])

thrust_chamber = psf.regen.ThrustChamber(contour=contour, 
                                         wall_group=wall_group,
                                         combustion_transport=aerothermodynamics,  
                                         cooling_circuit_group=cooling_circuit_group,
                                         roughness=params["roughness_height"],
                                         h_gas_corr=1.0, # No correction applied
                                         h_cold_corr=1.0) # No correction applied

mdot_fu = aerothermodynamics.mdot_fu
boundary_conditions_a = psf.regen.BoundaryConditions(T_coolant_in = params["T_coolant_in"], 
                                                     p_coolant_in = params["p_coolant_in"], 
                                                     mdot_coolant = mdot_fu)

cooling_data_a = psf.regen.steady_heating_analysis(thrust_chamber, 
                                                   n_nodes = 30, 
                                                   circuit_index=0, 
                                                   boundary_conditions=boundary_conditions_a, 
                                                   solver="newton", 
                                                   output=True)

T_out = cooling_data_a['T_static'][-1] # TODO: check if this should actually be total values, i forget
p_out = cooling_data_a['p_static'][-1]

boundary_conditions_b = psf.regen.BoundaryConditions(T_coolant_in = T_out, 
                                                     p_coolant_in = p_out, 
                                                     mdot_coolant = mdot_fu)

cooling_data_b = psf.regen.steady_heating_analysis(thrust_chamber, 
                                                   n_nodes = 80, 
                                                   circuit_index=1, 
                                                   boundary_conditions=boundary_conditions_b, 
                                                   solver="newton", 
                                                   output=True)

# Save results
script_dir = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(script_dir, "regen_results.pkl")

print(cooling_data_b)

res = psf.common.Results()
res.add(name="params", obj=params)
res.add(name="thrust_chamber", obj=thrust_chamber)
res.add(name="cooling_data", obj=[cooling_data_a, cooling_data_b])
res.save(out_path)