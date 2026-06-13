# Welcome to the Heimdall rocket engine simulation


import numpy as np
import os
import json
import time
from typing import Callable, Iterable, Union
import warnings

start = time.time()
import pyskyfire as psf
end = time.time()
# The actual added 3% TEOS contains enough Si to create approx 1% SiO2, which is added as a substitute for TEOS
# The remaining organic part of the TEOS is just approximated as ethanol

#print(f"Import took {end - start} s")

# Engine inputs 

params = dict(
    p_c = 21e5, # Measured static pressure during hotfire
    r_t = 30.888e-3, # Measured from Catia
    r_c = 64.972e-3, # Catia
    r_e = 56.876e-3, # Catia
    L_star = 1.927, # calculated from geometry
    MR = 3.0, # Best guess from Gaute
    L_c = 0.210,
    theta_conv = 30,
    R_1f = 1.5,
    R_2f = 3,
    R_3f = 0.5,
    length_fraction = 0.8,
    F = 8.8e3,
    n_channels = 15,
    T_coolant_in = 298.15,
    p_coolant_in = 28e5,
    T_ox_chamber_in = 300,
    T_fu_chamber_in = 410,
    tube_roughness_height = 20e-6, # This is a guess
    fu_coolant  = psf.common.Fluid(type="fuel", propellants=["Ethanol"], fractions=[1.0]),
    fu          = psf.common.Fluid(type="fuel", propellants=["C2H5OH"], fractions=[1.0]),
    ox          = psf.common.Fluid(type="oxidizer", propellants=["N2O"], fractions=[1.0]),
    wall_teos   = psf.common.solids.TEOS,
    teos_thickness = 10e-6, # guesswork
    wall_inco   = psf.common.solids.Inconel625,
    inco_thickness = 0.7e-3,
)

# Do some additional calcs
params["A_t"] = np.pi*params["r_t"]**2
params["A_c"] = np.pi*params["r_c"]**2
params["A_e"] = np.pi*params["r_e"]**2
params["eps"] = params["A_e"]/params["A_t"]


# Generate the contour coordinates
start = time.time() 
xs, rs = psf.regen.contour.get_contour(L_c = params["L_c"], 
                                            r_t = params["r_t"], 
                                            area_ratio = params["eps"], 
                                            r_c = params["r_c"], 
                                            theta_conv = params["theta_conv"],
                                            nozzle = "rao", 
                                            R_1f=params["R_1f"],
                                            R_2f=params["R_2f"],
                                            R_3f=params["R_3f"],
                                            length_fraction = params["length_fraction"]
                                            )

# Make a contour object
contour = psf.regen.Contour(xs, rs, name = "Heimdall Contour")

wall_teos = psf.regen.Wall(material = params["wall_teos"], thickness=params["teos_thickness"]) # Guesswork as to thickness and conductivity
wall_inco = psf.regen.Wall(material = params["wall_inco"], thickness = params["inco_thickness"]) 


def channel_height(x): 
    return 2e-3

#cross_section = psf.regen.CrossSectionSquared() 
cross_section = psf.regen.CrossSectionRounded()


coolant_transport = psf.skycea.CoolantTransport(params["fu_coolant"])

pass1 = psf.regen.CoolingCircuit(name="Half Pass",
                                    contour=contour, 
                                     coolant_transport=coolant_transport, 
                                     cross_section=cross_section, 
                                     span = [-1.0, 1.0], 
                                     placement=psf.regen.SurfacePlacement(n_channel_positions=20), 
                                     channel_height=channel_height, 
                                     walls=[wall_teos, wall_inco],
                                     roughness=params["tube_roughness_height"],
                                     blockage_ratio=0.1)

"""pass2 = psf.regen.CoolingCircuit(name="Half Pass",
                                    contour=contour, 
                                     coolant_transport=coolant_transport, 
                                     cross_section=cross_section, 
                                     span = [-0.51, -1.0], 
                                     placement=psf.regen.SurfacePlacement(n_channel_positions=20), 
                                     channel_height=channel_height, 
                                     walls=[wall_teos, wall_inco],
                                     roughness=params["tube_roughness_height"],
                                     blockage_ratio=0.2)"""

pass3 = psf.regen.CoolingCircuit(name="Full Pass", 
                                     contour=contour, 
                                     coolant_transport=coolant_transport, 
                                     cross_section=cross_section, 
                                     span = [1.0, -0.5], 
                                     placement=psf.regen.SurfacePlacement(n_channel_positions=20), 
                                     channel_height=channel_height, 
                                     walls=[wall_teos, wall_inco],
                                     roughness=params["tube_roughness_height"],
                                     blockage_ratio=0.1)



cooling_circuit_group = psf.regen.CoolingCircuitGroup(circuit_list=[pass1, pass3])


aerothermodynamics = psf.skycea.Aerothermodynamics.from_F_eps_Lstar(fu=params["fu"], ox=params["ox"], T_fu_in=params["T_fu_chamber_in"], T_ox_in=params["T_ox_chamber_in"], 
                                                                    MR=params["MR"], p_c=params["p_c"], F=params["F"], eps=params["eps"], L_star=params["L_star"], npts=15)

thrust_chamber = psf.regen.ThrustChamber(contour=contour, 
                                         combustion_transport=aerothermodynamics,  
                                         cooling_circuit_group=cooling_circuit_group, 
                                         h_gas_corr=1.0)


mdot_fu = aerothermodynamics.mdot_fu
boundary_conditions_a = psf.regen.BoundaryConditions(T_coolant_in = params["T_coolant_in"], p_coolant_in = params["p_coolant_in"], mdot_coolant = mdot_fu)

cooling_data_a = psf.regen.steady_heating_analysis(thrust_chamber, n_nodes = 50, circuit_index=0, boundary_conditions=boundary_conditions_a, solver="newton", output=True)

T_out = cooling_data_a["T_stagnation"][-1]
p_out = cooling_data_a["p_stagnation"][-1]

boundary_conditions_b = psf.regen.BoundaryConditions(T_coolant_in=T_out, p_coolant_in=p_out, mdot_coolant=mdot_fu)

cooling_data_b = psf.regen.steady_heating_analysis(thrust_chamber, n_nodes = 50, circuit_index=1, boundary_conditions=boundary_conditions_b, solver="newton", output=True)


# Save results
script_dir = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(script_dir, "regen.pkl")
res = psf.common.Results()
res.add(name="params", obj=params)
res.add(name="cooling_data", obj=[cooling_data_a, cooling_data_b])
res.add(name="thrust_chamber", obj=thrust_chamber)
res.save(out_path)


