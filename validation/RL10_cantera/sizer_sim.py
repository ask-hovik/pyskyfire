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

    # ==== LH2 side ====
    st["fu_engine_in"]              = psf.common.Station(params["p_tank_fu"], params["T_tank_fu"], params["thrust_chamber"].combustion_transport.mdot_fu)
    st["stage1_pump_in"]            = psf.common.Station(st["fu_engine_in"].p, st["fu_engine_in"].T, 2.8)
    st["stage1_shaft_recirc"]       = psf.common.Station(params["p_c"]*1.20, st["stage1_pump_in"].T, 2.8)
    st["stage2_shaft_recirc"]       = psf.common.Station(params["p_c"]*1.20, st["stage1_pump_in"].T, 2.8)
    st["pump_interstage1"]          = psf.common.Station(params["p_c"]*1.20, st["stage1_pump_in"].T, 2.8)
    st["pump_interstage2"]          = psf.common.Station(params["p_c"]*1.20, st["stage1_pump_in"].T, 2.8)
    st["stage2_pump_out"]           = psf.common.Station(params["p_c"]*2.30, st["pump_interstage2"].T, 2.8)
    st["regen_duct_in"]             = psf.common.Station(params["p_c"]*2.30, st["pump_interstage2"].T, 2.8)
    st["regen_in"]                  = psf.common.Station(st["stage2_pump_out"].p, st["stage2_pump_out"].T, 2.8)
    st["regen_interstage"]          = psf.common.Station(st["regen_in"].p, st["regen_in"].T+30.0, 2.8)
    st["regen_out"]                 = psf.common.Station(params["p_c"]*1.80, st["regen_interstage"].T+150.0, 2.8)
    st["bypass_in"]                 = psf.common.Station(params["p_c"]*1.70, st["regen_interstage"].T+150.0, 2.8)
    st["bypass_valve"]              = psf.common.Station(params["p_c"]*1.60, st["regen_interstage"].T+150.0, 2.8)
    st["turbine_in"]                = psf.common.Station(st["regen_out"].p, st["regen_out"].T, 2.8)
    st["turbine_out"]               = psf.common.Station(params["p_c"],       st["turbine_in"].T/1.10, 2.8)
    st["bypass_out"]                = psf.common.Station(params["p_c"],       st["turbine_in"].T/1.10, 2.8)
    st["fu_injector_plenum"]        = psf.common.Station(params["p_c"],       st["turbine_in"].T/1.10, 2.8)
    st["fu_chamber_in"]             = psf.common.Station(params["p_c"],       st["turbine_out"].T, 2.8)

    # ==== LOX side ====
    st["ox_engine_in"]              = psf.common.Station(params["p_tank_ox"], params["T_tank_ox"], params["thrust_chamber"].combustion_transport.mdot_ox)
    st["ox_pump_in"]                = psf.common.Station(st["ox_engine_in"].p, st["ox_engine_in"].T, 14.0)
    st["ox_pump_out"]               = psf.common.Station(params["p_c"]*1.30, st["ox_pump_in"].T, 14.0)
    st["ox_shaft_recirc"]           = psf.common.Station(params["p_c"]*1.30, st["ox_pump_in"].T, 14.0)
    st["ox_duct_in"]                = psf.common.Station(params["p_c"]*1.30, st["ox_pump_in"].T, 14.0)
    st["ox_injector_plenum"]        = psf.common.Station(params["p_c"]*1.30, st["ox_pump_in"].T, 14.0)
    st["ox_chamber_in"]             = psf.common.Station(params["p_c"],       st["ox_pump_out"].T, 14.0)

    return st

# ------- 2. helper: initialise scalar signals ------------------------------
def setup_initial_signals(params):
    sg = {}

    # boundary conditions
    sg["p_c"]            = params["p_c"]

    # bootstrap guesses
    sg["P_stage1_fuel_pump"] = 2.8e5
    sg["P_stage2_fuel_pump"] = 2.0e5
    sg["P_ox_pump"]         = 1.0e5
    sg["P_required"]        = sg["P_stage1_fuel_pump"] + sg["P_stage2_fuel_pump"] + sg["P_ox_pump"]

    return sg
def setup_thrust_chamber(params):
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

    # Define a channel height function using the built in maker
    channel_height_fn = psf.regen.make_channel_height_fn(
        contour=contour, 
        region_fractions=[-1.0, 0.25, 1.0], 
        flat_heights= [params["ht_up"], params["ht_lo"]], 
        pinch_factors= [params["pn_up"], params["pn_lo"]], 
        transition_widths=[0.1]
    ) # total coolant volume should be ca 0.015831543m3


    cross_section = psf.regen.CrossSectionRounded()
    LH2_transport = psf.skycea.CoolantTransport(params["coolant_fu"])

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
    blocks.append(psf.common.MassFlowMergerBlock(name="merge_stage1_recirc", 
                                                 st_in=["fu_engine_in", "stage1_shaft_recirc", "stage2_shaft_recirc"],
                                                 st_out="stage1_pump_in",
                                                 medium=params["coolprop_fu"]))
    
    blocks.append(psf.common.PumpBlock(name="stage1_fuel_pump", 
                                       st_in="stage1_pump_in",
                                       st_out="pump_interstage1",
                                       overcome=["duct_pump_regen", "regen_half_pass", "regen_full_pass", "duct_regen_turbine", "turbine", "duct_turbine_injector", "fu_injector"],
                                       load_fraction= 0.5,
                                       p_base = params["p_c"],
                                       input_p = params["p_tank_fu"],
                                       eta=params["eta_pump_fu"], 
                                       n = params["n_fu"], 
                                       medium=params["coolprop_fu"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="split_stage1_recirc", 
                                                   st_in="pump_interstage1",
                                                   st_out=["pump_interstage2", "stage1_shaft_recirc"],
                                                   fractions=[1-params["zeta_stage1_recirc"], params["zeta_stage1_recirc"]],
                                                   medium=params["coolprop_fu"]))
    

    blocks.append(psf.common.PumpBlock(name="stage2_fuel_pump", 
                                       st_in="pump_interstage2", 
                                       st_out="stage2_pump_out", 
                                       overcome=["duct_pump_regen", "regen_half_pass", "regen_full_pass", "duct_regen_turbine", "turbine", "duct_turbine_injector", "fu_injector"],
                                       load_fraction= 0.5,
                                       p_base= params["p_c"],
                                       input_p = params["p_tank_fu"],
                                       eta=params["eta_pump_fu"], 
                                       n = params["n_fu"], 
                                       medium=params["coolprop_fu"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="split_stage2_recirc", 
                                                   st_in="stage2_pump_out",
                                                   st_out=["regen_duct_in", "stage2_shaft_recirc", "gearbox_dump"],
                                                   fractions=[1-params["zeta_stage2_recirc"]-params["zeta_stage2_gearbox"], params["zeta_stage2_recirc"], params["zeta_stage2_gearbox"]],
                                                   medium=params["coolprop_fu"]))
    
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_pump_regen",  
                                             st_in="regen_duct_in", 
                                             st_out="regen_in", 
                                             eta=params["eta_pump_regen_fu"],
                                             medium=params["coolprop_fu"]))
    
    blocks.append(psf.common.RegenBlock(name="regen_half_pass", 
                                        st_in="regen_in", 
                                        st_out="regen_interstage", 
                                        circuit_index=0, 
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["coolprop_fu"]))
    
    blocks.append(psf.common.RegenBlock(name="regen_full_pass", 
                                        st_in="regen_interstage", 
                                        st_out="regen_out", 
                                        circuit_index=1, 
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["coolprop_fu"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_regen_turbine", 
                                             st_in="regen_out", 
                                             st_out="bypass_in", 
                                             eta=params["eta_regen_turbine_fu"],
                                             medium=params["coolprop_fu"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="split_turbine_bypass", 
                                                   st_in="bypass_in",
                                                   st_out=["turbine_in", "bypass_valve"],
                                                   fractions=[1-params["zeta_turbine_bypass"], params["zeta_turbine_bypass"]],
                                                   medium=params["coolprop_fu"]))
    
    
    blocks.append(psf.common.TurbineBlock(name="turbine",  
                                          st_in="turbine_in", 
                                          st_out="bypass_out", 
                                          P_req_key="P_required",
                                          eta=params["eta_turbine_fu"], 
                                          medium=params["coolprop_fu"]))
    
    blocks.append(psf.common.MassFlowMergerBlock(name="merge_turbine_bypass", 
                                                 st_in=["bypass_out", "bypass_valve"],
                                                 st_out="turbine_out",
                                                 medium=params["coolprop_fu"],))
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_turbine_injector", 
                                             st_in="turbine_out", 
                                             st_out="fu_injector_plenum", 
                                             eta=params["eta_turbine_injector_fu"],
                                             medium=params["coolprop_fu"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="fu_injector", 
                                             st_in="fu_injector_plenum", 
                                             st_out="fu_chamber_in", 
                                             eta=params["eta_fu_injector"],
                                             medium=params["coolprop_fu"]))
    


    # ========= Ox side blocks =========
    blocks.append(psf.common.MassFlowMergerBlock(name="merge_ox_recirc", 
                                                 st_in=["ox_engine_in", "ox_shaft_recirc"],
                                                 st_out="ox_pump_in",
                                                 medium=params["coolprop_ox"]))
    
    blocks.append(psf.common.PumpBlock(name="ox_pump", 
                                       st_in="ox_pump_in", 
                                       st_out="ox_pump_out",
                                       overcome=["duct_pump_chamber_ox", "ox_injector"],
                                       load_fraction=1.0,
                                       p_base=params["p_c"],
                                       input_p = params["p_tank_ox"],
                                       n=params["n_ox"], eta=params["eta_pump_ox"],
                                       medium=params["coolprop_ox"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="split_ox_recirc", 
                                                   st_in="ox_pump_out",
                                                   st_out=["ox_duct_in", "ox_shaft_recirc"],
                                                   fractions=[1-params["zeta_pump_ox_recirc"], params["zeta_pump_ox_recirc"]],
                                                   medium=params["coolprop_ox"]))
    

    blocks.append(psf.common.SimpleDuctBlock(name="duct_pump_chamber_ox", 
                                             st_in="ox_duct_in", 
                                             st_out="ox_injector_plenum",
                                             eta=params["eta_pump_injector_ox"],
                                             medium=params["coolprop_ox"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="ox_injector", 
                                             st_in="ox_injector_plenum", 
                                             st_out="ox_chamber_in",
                                             eta=params["eta_ox_injector"],
                                             medium=params["coolprop_ox"]))
    

    blocks.append(psf.common.TransmissionBlock(name="transmission", sink_keys=["P_stage1_fuel_pump", "P_stage2_fuel_pump", "P_ox_pump"]))

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
    params = {
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

        # Turbine efficiencies
        'eta_turbine_fu': 0.7353,

        # Duct efficiencies fuel side
        'eta_stage1_stage2_fu': 0.95,
        'eta_pump_regen_fu': 0.95,
        'eta_regen_turbine_fu': 0.98,
        'eta_turbine_injector_fu': 0.94,
        'eta_fu_injector': 0.88,

        # Duct efficiencies ox side
        'eta_pump_regen_ox': 0.95,
        'eta_regen_turbine_ox': 0.95,
        'eta_pump_injector_ox': 0.88,
        'eta_ox_injector': 0.88,

        # Mass flow leakage fractions
        'zeta_stage1_recirc': 0.0070,
        'zeta_stage2_recirc': 0.00294, 
        'zeta_stage2_gearbox': 0.00294,
        'zeta_turbine_bypass': 0.00427,
        'zeta_pump_ox_recirc': 0.0005 ,

        
    }

    # Additional derived parameters
    params['T_tank_ox'] = 97.05 # Taken from Binder et al.
    params['T_tank_fu'] = CP.PropsSI('T', 'P', params['p_tank_fu'], 'Q', 0, params["coolprop_fu"]) 

    
    params['rho_ox_tank'] = CP.PropsSI('D', 'P', params['p_tank_ox'], 'Q', 0, params["coolprop_ox"])
    params['rho_fu_tank'] = CP.PropsSI('D', 'P', params['p_tank_fu'], 'Q', 0, params["coolprop_fu"])

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