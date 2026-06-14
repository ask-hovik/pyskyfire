# -------------------------------------------------------------
# RL10_engine_sizer.py
# -------------------------------------------------------------
import numpy as np
import CoolProp.CoolProp as CP
import pyskyfire as psf
import os
import time


def setup_initial_stations(params):
    
    # Create the stations used in the engine representation. Seed with initial values
    st = {}

    # ==== CH4 side ====
    mdot_fu_est = params["thrust_chamber"].combustion_transport.mdot_fu
    st["fu_engine_in"]              = psf.common.Station(params["p_tank_fu"],   params["T_tank_fu"],        mdot_fu_est)
    st["fu_pump_in"]                = psf.common.Station(st["fu_engine_in"].p,  st["fu_engine_in"].T,       mdot_fu_est)
    st["fu_pump_out"]               = psf.common.Station(params["p_c"]*1.4,     params["T_tank_fu"]+10,     mdot_fu_est)
    st["fu_shaft_recirc"]           = psf.common.Station(params["p_c"]*1.4,     params["T_tank_fu"]+10,     mdot_fu_est)
    st["fu_regen_duct_in"]          = psf.common.Station(params["p_c"]*1.4,     params["T_tank_fu"]+10,     mdot_fu_est)
    st["fu_regen_in"]               = psf.common.Station(params["p_c"]*1.4,     params["T_tank_fu"]+10,     mdot_fu_est)
    st["fu_regen_interstage_1"]     = psf.common.Station(params["p_c"]*1.4,     params["T_tank_fu"]+100,     mdot_fu_est)
    st["fu_regen_interstage_2"]     = psf.common.Station(params["p_c"]*1.4,     params["T_tank_fu"]+180,     mdot_fu_est)
    st["fu_regen_out"]              = psf.common.Station(params["p_c"]*1.4,     params["T_tank_fu"]+250,     mdot_fu_est)
    st["fu_tubrine_inlet_split"]    = psf.common.Station(params["p_c"]*1.3,     params["T_tank_fu"]+250,     mdot_fu_est)
    st["fu_bypass_valve"]           = psf.common.Station(params["p_c"]*1.3,    params["T_tank_fu"]+250,     mdot_fu_est)
    st["fu_turbine_in"]             = psf.common.Station(params["p_c"]*1.2,     params["T_tank_fu"]+250,     mdot_fu_est)
    st["fu_turbine_out"]            = psf.common.Station(params["p_c"]*1.1,     params["T_tank_fu"]+200,     mdot_fu_est)
    st["fu_turbine_outlet_merge"]   = psf.common.Station(params["p_c"]*1.1,     params["T_tank_fu"]+200,     mdot_fu_est)
    st["fu_injector_plenum"]        = psf.common.Station(params["p_c"]*1.1,     params["T_tank_fu"]+200,     mdot_fu_est)
    st["fu_chamber_in"]             = psf.common.Station(params["p_c"],         params["T_tank_fu"]+200,     mdot_fu_est)

    

    # ==== LOX side ====
    mdot_ox_est = params["thrust_chamber"].combustion_transport.mdot_ox
    st["ox_engine_in"]              = psf.common.Station(params["p_tank_ox"],   params["T_tank_ox"],      mdot_ox_est)
    st["ox_pump_in"]                = psf.common.Station(st["ox_engine_in"].p,   st["ox_engine_in"].T,    mdot_ox_est)
    st["ox_pump_out"]               = psf.common.Station(params["p_c"]*1.6,     params["T_tank_ox"]+10,     mdot_ox_est)
    st["ox_shaft_recirc"]           = psf.common.Station(params["p_c"]*1.6,     params["T_tank_ox"]+10,     mdot_ox_est)
    st["ox_regen_duct_in"]          = psf.common.Station(params["p_c"]*1.6,     params["T_tank_ox"]+10,     mdot_ox_est)
    st["ox_regen_in"]               = psf.common.Station(params["p_c"]*1.6,     params["T_tank_ox"]+10,     mdot_ox_est)
    st["ox_regen_1_in"]             = psf.common.Station(params["p_c"]*1.6,     params["T_tank_ox"]+250,     mdot_ox_est/2)
    st["ox_regen_2_in"]             = psf.common.Station(params["p_c"]*1.6,     params["T_tank_ox"]+250,     mdot_ox_est/2)
    st["ox_regen_1_out"]            = psf.common.Station(params["p_c"]*1.6,     params["T_tank_ox"]+250,     mdot_ox_est/2)
    st["ox_regen_2_out"]            = psf.common.Station(params["p_c"]*1.6,     params["T_tank_ox"]+250,     mdot_ox_est/2)
    st["ox_regen_out"]              = psf.common.Station(params["p_c"]*1.6,     params["T_tank_ox"]+250,     mdot_ox_est/2)
    st["ox_turbine_inlet_split"]    = psf.common.Station(params["p_c"]*1.5,     params["T_tank_ox"]+250,    mdot_ox_est)
    st["ox_bypass_valve"]           = psf.common.Station(params["p_c"]*1.3,     params["T_tank_ox"]+250,     mdot_ox_est)
    st["ox_turbine_in"]             = psf.common.Station(params["p_c"]*1.3,     params["T_tank_ox"]+250,     mdot_ox_est)
    st["ox_turbine_out"]            = psf.common.Station(params["p_c"]*1.2,     params["T_tank_ox"]+200,     mdot_ox_est)
    st["ox_turbine_outlet_merge"]   = psf.common.Station(params["p_c"]*1.1,     params["T_tank_ox"]+200,    mdot_ox_est)
    st["ox_injector_plenum"]        = psf.common.Station(params["p_c"]*1.1,     params["T_tank_ox"]+200,    mdot_ox_est)
    st["ox_chamber_in"]             = psf.common.Station(params["p_c"],         params["T_tank_ox"]+200,    mdot_ox_est)

    return st

# ------- 2. helper: initialise scalar signals ------------------------------
def setup_initial_signals(params):
    sg = {}

    # boundary conditions
    sg["p_c"]            = params["p_c"] # TODO: seems unnecessary.

    # bootstrap guesses
    sg["P_fuel_pump"] = 2.8e5
    sg["P_ox_pump"] = 2.0e5
    sg["P_fuel_turbine_required"]        = sg["P_fuel_pump"]
    sg["P_ox_turbine_required"]         = sg["P_ox_pump"]

    return sg
def setup_thrust_chamber(params):
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
        return 2e-3

    def internal_width(x):
        return 20e-3

    def internal_height(x):
        return 20e-3


    cross_section = psf.regen.CrossSectionRounded()
    fuel_transport = psf.skycea.CoolantTransport(params["coolprop_fu"])
    ox_transport = psf.skycea.CoolantTransport(params["coolprop_ox"])
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
                                        placement=psf.regen.SurfacePlacement(n_channel_positions=300),
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
                                            cooling_circuits = [fu_copper_pass, fu_inco_pass_cocurrent, fu_inco_pass_countercurrent, ox_copper_pass],
                                            h_gas_corr=1.0, # No correction applied
                                            h_cold_corr=1.0) # No correction applied

    """plot_3d, viewer = psf.viz.make_engine_3d(thrust_chamber)
    plot_3d.show()
    del plot_3d 
    input()"""

    
    return thrust_chamber


# ------- 4. build network ---------------------------------------------------
def engine_sizer(params):

    # Make thrust chamber and add to the params dictionary. Initialise stations, signals and blocks
    thrust_chamber = setup_thrust_chamber(params)
    params["thrust_chamber"] = thrust_chamber
    stations = setup_initial_stations(params)
    signals  = setup_initial_signals(params)
    blocks   = []

    
    # Create engine components. Connect engine network using constructor
    # ======= Fuel side blocks ========
    blocks.append(psf.common.MassFlowMergerBlock(name="fuel_inlet_merge", 
                                                 st_in=["fu_engine_in", "fu_shaft_recirc"],
                                                 st_out="fu_pump_in",
                                                 medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.PumpBlock(name="fuel_pump", 
                                       st_in="fu_pump_in",
                                       st_out="fu_pump_out",
                                       overcome=["duct_pump_regen_fuel", "regen_throat_pass",  "duct_regen_turbine_fuel", "fuel_turbine", "duct_turbine_injector_fuel", "fu_injector"], # "regen_cocurrent_nozzle_pass", "regen_countercurrent_nozzle_pass",
                                       load_fraction= 1.0,
                                       p_base = params["p_c"],
                                       input_p = params["p_tank_fu"],
                                       eta=params["eta_pump_fu"], 
                                       n = params["n_fu"], 
                                       medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="fuel_recirc_split", 
                                                   st_in="fu_pump_out",
                                                   st_out=["fu_regen_duct_in", "fu_shaft_recirc"],
                                                   fractions=[1-params["zeta_fu_recirc"], params["zeta_fu_recirc"]],
                                                   medium=params["coolprop_fu"].propellants[0]))    
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_pump_regen_fuel",  
                                             st_in="fu_regen_duct_in", 
                                             st_out="fu_regen_in", 
                                             pressure_ratio=params["eta_pump_regen_fu"],
                                             medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.RegenBlock(name="regen_throat_pass", 
                                        st_in="fu_regen_in", 
                                        st_out="fu_regen_interstage_1", 
                                        circuit_index=0, 
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.RegenBlock(name="regen_cocurrent_nozzle_pass", 
                                        st_in="fu_regen_interstage_1", 
                                        st_out="fu_regen_interstage_2", 
                                        circuit_index=1, 
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.RegenBlock(name="regen_countercurrent_nozzle_pass", 
                                        st_in="fu_regen_interstage_2", 
                                        st_out="fu_regen_out", 
                                        circuit_index=2, 
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_regen_turbine_fuel", 
                                             st_in="fu_regen_out", 
                                             st_out="fu_turbine_inlet_split", 
                                             pressure_ratio=params["eta_regen_turbine_fu"],
                                             medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="fuel_split_turbine_bypass", 
                                                   st_in="fu_turbine_inlet_split",
                                                   st_out=["fu_turbine_in", "fu_bypass_valve"],
                                                   fractions=[1-params["zeta_fu_turbine_bypass"], params["zeta_fu_turbine_bypass"]],
                                                   medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.TurbineBlock(name="fuel_turbine",  
                                          st_in="fu_turbine_in", 
                                          st_out="fu_turbine_out", 
                                          P_req_key="P_fuel_turbine_required",
                                          eta=params["eta_turbine_fu"], 
                                          medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.MassFlowMergerBlock(name="merge_turbine_bypass", 
                                                 st_in=["fu_turbine_out", "fu_bypass_valve"],
                                                 st_out="fu_turbine_outlet_merge",
                                                 medium=params["coolprop_fu"].propellants[0],))
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_turbine_injector_fuel", 
                                             st_in="fu_turbine_outlet_merge", 
                                             st_out="fu_injector_plenum", 
                                             pressure_ratio=params["eta_turbine_injector_fu"],
                                             medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="fu_injector", 
                                             st_in="fu_injector_plenum", 
                                             st_out="fu_chamber_in", 
                                             pressure_ratio=params["eta_fu_injector"],
                                             medium=params["coolprop_fu"].propellants[0]))
    
    blocks.append(psf.common.TransmissionBlock(name="fuel_shaft", sink_keys=["P_fuel_pump"], source_keys=["P_fuel_turbine_required"]))


    # ========= Ox side blocks =========
    blocks.append(psf.common.MassFlowMergerBlock(name="merge_ox_recirc", 
                                                 st_in=["ox_engine_in", "ox_shaft_recirc"],
                                                 st_out="ox_pump_in",
                                                 medium=params["coolprop_ox"].propellants[0]))
    
    blocks.append(psf.common.PumpBlock(name="ox_pump", 
                                       st_in="ox_pump_in", 
                                       st_out="ox_pump_out",
                                       overcome=["duct_pump_regen_ox", "ox_regen", "duct_regen_turbine_ox", "ox_turbine", "duct_turbine_injector_ox", "ox_injector"],
                                       load_fraction=1.0,
                                       p_base=params["p_c"],
                                       input_p = params["p_tank_ox"],
                                       n=params["n_ox"], eta=params["eta_pump_ox"],
                                       medium=params["coolprop_ox"].propellants[0]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="split_ox_recirc", 
                                                   st_in="ox_pump_out",
                                                   st_out=["ox_regen_duct_in", "ox_shaft_recirc"],
                                                   fractions=[1-params["zeta_ox_recirc"], params["zeta_ox_recirc"]],
                                                   medium=params["coolprop_ox"].propellants[0]))
    

    blocks.append(psf.common.SimpleDuctBlock(name="duct_pump_regen_ox", 
                                             st_in="ox_regen_duct_in", 
                                             st_out="ox_regen_in",
                                             pressure_ratio=params["eta_pump_regen_ox"],
                                             medium=params["coolprop_ox"].propellants[0]))

    blocks.append(psf.common.MassFlowSplitterBlock(name="ox_regen_split", 
                                                st_in="ox_regen_in",
                                                st_out=["ox_regen_1_in", "ox_regen_2_in"],
                                                fractions=[0.5, 0.5],
                                                medium=params["coolprop_ox"].propellants[0]))
    
    blocks.append(psf.common.RegenBlock(name="ox_regen",
                                        st_in="ox_regen_1_in",
                                        st_out="ox_regen_1_out",
                                        circuit_index=3, # Just pretending there is more area than there is
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["coolprop_ox"].propellants[0]))

    blocks.append(psf.common.RegenBlock(name="ox_regen",
                                        st_in="ox_regen_2_in",
                                        st_out="ox_regen_2_out",
                                        circuit_index=3, # here, same circuit index is deliberate in this sim
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["coolprop_ox"].propellants[0]))

    blocks.append(psf.common.MassFlowMergerBlock(name="ox_regen_merge", 
                                                 st_in=["ox_regen_1_out", "ox_regen_2_out"],
                                                 st_out="ox_regen_out",
                                                 medium=params["coolprop_ox"].propellants[0]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_regen_turbine_ox", 
                                             st_in="ox_regen_out", 
                                             st_out="ox_turbine_inlet_split",
                                             pressure_ratio=params["eta_pump_regen_ox"],
                                             medium=params["coolprop_ox"].propellants[0]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="ox_split_turbine_bypass", 
                                                st_in="ox_turbine_inlet_split",
                                                st_out=["ox_turbine_in", "ox_bypass_valve"],
                                                fractions=[1-params["zeta_ox_turbine_bypass"], params["zeta_ox_turbine_bypass"]],
                                                medium=params["coolprop_ox"].propellants[0]))
    
    blocks.append(psf.common.TurbineBlock(name="ox_turbine",  
                                          st_in="ox_turbine_in", 
                                          st_out="ox_turbine_out", 
                                          P_req_key="P_ox_turbine_required",
                                          eta=params["eta_turbine_ox"], 
                                          medium=params["coolprop_ox"].propellants[0]))
    
    blocks.append(psf.common.MassFlowMergerBlock(name="merge_ox_bypass", 
                                                 st_in=["ox_turbine_out", "ox_bypass_valve"],
                                                 st_out="ox_turbine_outlet_merge",
                                                 medium=params["coolprop_ox"].propellants[0]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_turbine_injector_ox", 
                                             st_in="ox_turbine_outlet_merge", 
                                             st_out="ox_injector_plenum",
                                             pressure_ratio=params["eta_pump_regen_ox"],
                                             medium=params["coolprop_ox"].propellants[0]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="ox_injector", 
                                             st_in="ox_injector_plenum", 
                                             st_out="ox_chamber_in",
                                             pressure_ratio=params["eta_ox_injector"],
                                             medium=params["coolprop_ox"].propellants[0]))
    

    blocks.append(psf.common.TransmissionBlock(name="ox_shaft", sink_keys=["P_ox_pump"], source_keys=["P_ox_turbine_required"]))

    # Initialise signals for pressure drop across each component
    for blk in blocks:
        if hasattr(blk, "dp_key"):
            st_in  = blk.station_inputs[0]  if blk.station_inputs  else None
            st_out = blk.station_outputs[0] if blk.station_outputs else None

            if st_in in stations and st_out in stations:
                dp0 = max(stations[st_in].p - stations[st_out].p, 0.0)
            else:
                dp0 = 0.0 # fall-back

            signals.setdefault(blk.dp_key, dp0) # TODO: internalise this or something. Very hacky. 

    # Instantiate simulation object and run simulation
    net = psf.common.EngineNetwork(stations, signals, blocks)
    net.run_fixed_point(tol=1e-5, max_iter=100)

    return net.stations, net.signals, net.residuals, net.block_results, params, net


def main():
    """params = {
        # Core engine parameters
        'p_c': 32.7501e5,
        'MR': 5.055,
        'AR_c': 1.5,
        'F': 73.182e3,
        'p_e': 0.0377e5,
        'A_e': 0.708,
        'p_tank_fu': 1.908e5,
        'p_tank_ox': 2.92e5,
        'L_star': 0.95,
        "eps": 56.1,  # Without silver insert

        # Propellants
        "cea_fu": psf.common.Fluid(type="fuel", propellants=["H2"], fractions=[1.0]),
        "cea_ox": psf.common.Fluid(type="oxidizer", propellants=["O2(L)"], fractions=[1.0]),
        "coolant_fu": psf.common.Fluid(type="fuel", propellants=["Hydrogen"], fractions=[1.0]),
        "coolprop_fu": "hydrogen",
        "coolprop_ox": "oxygen",

        "T_gas_fu_in": 200,  # Adjusted by a few kelvin to fall within NASA CEA tables
        "T_gas_ox_in": 100,  # Also adjusted

        # Chamber geometry parameters
        'theta_conv': 25,
        'r_c': 0.123, # measured
        'R_1f': 1.5,
        'R_2f': 3.0,
        'R_3f': 0.5,
        'length_fraction': 0.713,
        'n_fu': 31537,
        'n_ox': 12615,

        # Cooling channel parameters
        "wall_thickness": 0.31e-3,
        "roughness_height": 1.1684e-6,  # found in Binder et al 4.6e-5 inches. Conventional SS tubing is often cited as 0.0015e-3
        'ht_up': 0.0032,
        'ht_lo': 0.00134,
        'pn_up': 0.7,
        'pn_lo': -5.0,

        # Pump efficiencies
        'eta_pump_fu': 0.5810,
        'eta_pump_ox': 0.6422,

        # Pump load fraction
        'stage1_fraction': 0.5,
        
    }"""

    params = dict(
        p_c = 100e5,
        p_e = 0.8e5,
        MR = 2.8,
        AR_c = 1.8,
        F = 300e3,

        # Fuel/oxidizer parameters
        cea_fu = psf.common.Fluid(type="fuel", propellants=["CH4"], fractions=[1.0]),
        cea_ox = psf.common.Fluid(type="oxidizer", propellants=["O2"], fractions=[1.0]),
        coolprop_fu =  psf.common.Fluid(type="fuel", propellants=["methane"], fractions=[1.0]),
        coolprop_ox = psf.common.Fluid(type="oxidizer", propellants=["oxygen"], fractions=[1.0]),
        T_gas_fu_in = 300,  # Adjusted by a few kelvin to fall within NASA CEA tables
        T_gas_ox_in = 300,  # Also adjusted

        # Propellant tanks
        p_tank_ox = 5e5,
        p_tank_fu = 5e5,

        # Chamber-nozzle parameters
        theta_conv = 35,
        R_1f = 1.5,
        R_2f = 2,
        R_3f = 0.5,
        length_fraction = 0.8, # 80% nozzle
        L_star = 1.1,  

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
        inconel = psf.common.solids.Inconel625,

        # Pump paramters
        eta_pump_fu = 0.6,
        n_fu = 50000, #rpm
        eta_pump_ox = 0.6,
        n_ox = 25000, #rpm

        # Turbine efficiencies
        eta_turbine_fu = 0.73,
        eta_turbine_ox = 0.73,

        # Duct efficiencies fuel side
        eta_pump_regen_fu = 0.95,
        eta_regen_turbine_fu = 0.98,
        eta_turbine_injector_fu = 0.94,
        eta_fu_injector = 0.88,

        # Duct efficiencies ox side
        eta_pump_regen_ox = 0.95,
        eta_regen_turbine_ox = 0.95,
        eta_pump_injector_ox = 0.88,
        eta_ox_injector = 0.88,

        # Mass flow leakage fractions
        zeta_fu_recirc = 0.005,
        zeta_ox_recirc = 0.005, 
        zeta_fu_turbine_bypass = 0.005,
        zeta_ox_turbine_bypass = 0.005 ,

    )

    # Additional derived parameters
    params['T_tank_ox'] = CP.PropsSI('T', 'P', params['p_tank_fu'], 'Q', 0, params["coolprop_ox"].propellants[0]) 
    params['T_tank_fu'] = CP.PropsSI('T', 'P', params['p_tank_fu'], 'Q', 0, params["coolprop_fu"].propellants[0]) 

    
    params['rho_ox_tank'] = CP.PropsSI('D', 'P', params['p_tank_ox'], 'Q', 0, params["coolprop_ox"].propellants[0])
    params['rho_fu_tank'] = CP.PropsSI('D', 'P', params['p_tank_fu'], 'Q', 0, params["coolprop_fu"].propellants[0])

    # Start clock
    start_time = time.time()

    # Run engine sizer
    stations, signals, residuals, block_results, input_params, net = engine_sizer(params)

    # Stop clock and print
    end_time = time.time()
    duration = end_time - start_time
    print(f"Sim took {duration} s")

    # Save results to a file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(script_dir, "sizer_results.pkl")

    res = psf.common.Results()
    res.add(name="net", obj=net)
    res.add(name="stations", obj=stations)
    res.add(name="signals", obj=signals)
    res.add(name="residuals", obj=residuals)
    res.add(name="block_results", obj=block_results)
    res.add(name="input_params", obj=input_params)
    res.save(out_path)

if __name__ == "__main__":
    main()