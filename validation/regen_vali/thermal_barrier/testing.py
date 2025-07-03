import pyskyfire as psf
import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import subprocess
import time


# The following are the inputs that generate a contour remarkably similar to the actual RL10 contour. 
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
T_coolant_in = 33 # Kelvin
p_coolant_in = 69e5 # 68 bar

optimals = psf.skycea.OptimalValues(ox["cea"], fu["cea"], F, MR, p_c, p_e, L_star)

V_c = optimals.V_c_opt
r_t = optimals.r_t_opt
eps = optimals.eps_opt

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

RL10_contour = psf.regen.Contour(RL10_xs, RL10_rs, name = "Replicate Contour")

#psf.regen.contour.plot_theta_vs_epsilon()

wall = psf.regen.Wall(material = psf.common.material.StainlessSteel304, thickness = 0.3e-3) 

wall_group = psf.regen.WallGroup(walls=[wall])

channel_height_fn = psf.regen.make_channel_height_fn(
    contour=RL10_contour, 
    region_fractions=[-1.0, 0.25, 1.0], 
    flat_heights= [0.03, 0.0134], 
    pinch_factors= [0.7, -5.0], 
    transition_widths=[0.1]
) # total coolant volume should be ca 0.015831543m3

cross_section = psf.regen.CrossSectionSquared()
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
                                         cooling_circuit_group=cooling_circuit_group)


"""psf.regen.plot.visualize_cooling_channels_vispy_6(
    thrust_chamber,
    draw_section_line=True,
    clip_plane=((-0.15,0,0),(0,1,0)),
    line_color=(0,0,0,1),
    draw_only_section_line=True, 
    show_axis=False
)""" # this is the one


#psf.regen.plot.plot_engine_3D(thrust_chamber, show_axis=False)
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

psf.regen.plot.plot_contour([thrust_chamber.contour])
psf.regen.plot.plot_coolant_pressure(cooling_data, cooling_data_2, static=True, stagnation=True)
psf.regen.plot.plot_coolant_temperature(cooling_data, cooling_data_2)
psf.regen.plot.plot_wall_temperature(cooling_data, cooling_data_2, plot_hot=True, plot_interfaces=True, plot_coolant_wall=True, circuits=[thrust_chamber.cooling_circuit_group.circuits[0], thrust_chamber.cooling_circuit_group.circuits[1]])
psf.regen.plot.plot_velocity(cooling_data, cooling_data_2)
psf.regen.plot.plot_heat_flux(cooling_data, cooling_data_2)
plt.show()

mdot_ox = thrust_chamber.combustion_transport.mdot_ox
mdot_fu = thrust_chamber.combustion_transport.mdot_fu
print(f"mdot_ox: {mdot_ox}")
print(f"mdot_fu: {mdot_fu}")