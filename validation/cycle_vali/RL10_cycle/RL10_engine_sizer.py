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
    st["fu_engine_in"]              = psf.common.Station(params["p_tank_fu"], params["T_tank_fu"], params["thrust_chamber"].optimal_values.mdot_fu)
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
    st["ox_engine_in"]              = psf.common.Station(params["p_tank_ox"], params["T_tank_ox"], params["thrust_chamber"].optimal_values.mdot_ox)
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

    # Get the optimal values for the input combination
    optimals = psf.skycea.OptimalValues(ox = params["ox"]["cea"], 
                                    fu = params["fu"]["cea"], 
                                    F = params["F"], 
                                    MR = params["MR"], 
                                    p_c = params["p_c"], 
                                    p_e = params["p_e"], 
                                    L_star = params["L_star"], 
                                    Isp=445.6) # Isp is taken from Binder et al. 
    
    V_c = optimals.V_c_opt
    r_t = optimals.r_t_opt
    eps = optimals.eps_opt
    p_e = optimals.p_e

    # Get contour coordinates
    xs, rs = psf.regen.contour_2.get_contour_2(V_c = V_c, 
                                        r_t = r_t, 
                                        area_ratio = eps, 
                                        AR_c=params["AR_c"], 
                                        theta_conv = params["theta_conv"],
                                        nozzle = "rao", 
                                        R_1f = params["R_1f"],
                                        R_2f = params["R_2f"],
                                        R_3f = params["R_3f"],
                                        length_fraction = params["length_fraction"]
                                        )

    # Create contour class
    contour = psf.regen.Contour(xs, rs, name = "Replicate Contour")

    # Create wall arrangement
    wall = psf.regen.Wall(material = psf.common.material.StainlessSteel304, thickness = params["wall_thickness"]) 
    wall_group = psf.regen.WallGroup(walls=[wall])

    # Define a channel height function using the built in maker
    channel_height_fn = psf.regen.make_channel_height_fn(
        contour=contour, 
        region_fractions=[-1.0, 0.25, 1.0], 
        flat_heights= [params["ht_up"], params["ht_lo"]], 
        pinch_factors= [params["pn_up"], params["pn_lo"]], 
        transition_widths=[0.1]
    )

    # Choose a cross section 
    cross_section = psf.regen.CrossSectionRounded()

    # Choose a coolant fluid 
    LH2_transport = psf.skycea.CoolantTransport(params["fu"]["coolprop"])

    # Create cooling circuits
    half_pass = psf.regen.CoolingCircuit(name="Half Pass",
                                        contour=contour, 
                                        coolant_transport=LH2_transport, 
                                        cross_section=cross_section, 
                                        span = [0.25, 1.0], 
                                        placement=psf.regen.SurfacePlacement(n_channel_positions=180), 
                                        channel_height=channel_height_fn)

    full_pass = psf.regen.CoolingCircuit(name="Full Pass",
                                        contour=contour, 
                                        coolant_transport=LH2_transport, 
                                        cross_section=cross_section, 
                                        span = [1.0, -1.0], 
                                        placement=psf.regen.SurfacePlacement(n_channel_positions=180), 
                                        channel_height=channel_height_fn)

    cooling_circuit_group = psf.regen.CoolingCircuitGroup(circuit_list=[half_pass, full_pass])

    # Create chamber combustion reaction
    combustion_transport = psf.skycea.CombustionTransport(ox = params["ox"]["cantera"], 
                                                        fu = params["fu"]["cantera"],
                                                        p_c = params["p_c"],
                                                        p_e = p_e,
                                                        MR = params["MR"],
                                                        T_fuel_in = 111, 
                                                        T_ox_in = 91)

    combustion_transport.set_essentials(mdot=optimals.mdot, mdot_fu=optimals.mdot_fu, mdot_ox=optimals.mdot_ox)

    # Create the thrust chamber object
    thrust_chamber = psf.regen.ThrustChamber(contour=contour, 
                                            wall_group=wall_group,
                                            combustion_transport=combustion_transport, 
                                            optimal_values=optimals, 
                                            cooling_circuit_group=cooling_circuit_group)
    
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
                                                 medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.PumpBlock(name="stage1_fuel_pump", 
                                       st_in="stage1_pump_in",
                                       st_out="pump_interstage1",
                                       overcome=["duct_pump_regen", "regen_half_pass", "regen_full_pass", "duct_regen_turbine", "turbine", "duct_turbine_injector", "fu_injector"],
                                       load_fraction= 0.5,
                                       p_base = params["p_c"],
                                       input_p = params["p_tank_fu"],
                                       eta=params["eta_pump_fu"], 
                                       n = params["n_fu"], 
                                       medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="split_stage1_recirc", 
                                                   st_in="pump_interstage1",
                                                   st_out=["pump_interstage2", "stage1_shaft_recirc"],
                                                   fractions=[1-params["zeta_stage1_recirc"], params["zeta_stage1_recirc"]],
                                                   medium=params["fu"]["coolprop"]))
    

    blocks.append(psf.common.PumpBlock(name="stage2_fuel_pump", 
                                       st_in="pump_interstage2", 
                                       st_out="stage2_pump_out", 
                                       overcome=["duct_pump_regen", "regen_half_pass", "regen_full_pass", "duct_regen_turbine", "turbine", "duct_turbine_injector", "fu_injector"],
                                       load_fraction= 0.5,
                                       p_base= params["p_c"],
                                       input_p = params["p_tank_fu"],
                                       eta=params["eta_pump_fu"], 
                                       n = params["n_fu"], 
                                       medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="split_stage2_recirc", 
                                                   st_in="stage2_pump_out",
                                                   st_out=["regen_duct_in", "stage2_shaft_recirc", "gearbox_dump"],
                                                   fractions=[1-params["zeta_stage2_recirc"]-params["zeta_stage2_gearbox"], params["zeta_stage2_recirc"], params["zeta_stage2_gearbox"]],
                                                   medium=params["fu"]["coolprop"]))
    
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_pump_regen",  
                                             st_in="regen_duct_in", 
                                             st_out="regen_in", 
                                             eta=params["eta_pump_regen_fu"],
                                             medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.RegenBlock(name="regen_half_pass", 
                                        st_in="regen_in", 
                                        st_out="regen_interstage", 
                                        circuit_index=0, 
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.RegenBlock(name="regen_full_pass", 
                                        st_in="regen_interstage", 
                                        st_out="regen_out", 
                                        circuit_index=1, 
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_regen_turbine", 
                                             st_in="regen_out", 
                                             st_out="bypass_in", 
                                             eta=params["eta_regen_turbine_fu"],
                                             medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="split_turbine_bypass", 
                                                   st_in="bypass_in",
                                                   st_out=["turbine_in", "bypass_valve"],
                                                   fractions=[1-params["zeta_turbine_bypass"], params["zeta_turbine_bypass"]],
                                                   medium=params["fu"]["coolprop"]))
    
    
    blocks.append(psf.common.TurbineBlock(name="turbine",  
                                          st_in="turbine_in", 
                                          st_out="bypass_out", 
                                          P_req_key="P_required",
                                          eta=params["eta_turbine_fu"], 
                                          medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.MassFlowMergerBlock(name="merge_turbine_bypass", 
                                                 st_in=["bypass_out", "bypass_valve"],
                                                 st_out="turbine_out",
                                                 medium=params["fu"]["coolprop"],))
    
    blocks.append(psf.common.SimpleDuctBlock(name="duct_turbine_injector", 
                                             st_in="turbine_out", 
                                             st_out="fu_injector_plenum", 
                                             eta=params["eta_turbine_injector_fu"],
                                             medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="fu_injector", 
                                             st_in="fu_injector_plenum", 
                                             st_out="fu_chamber_in", 
                                             eta=params["eta_fu_injector"],
                                             medium=params["fu"]["coolprop"]))
    


    # ========= Ox side blocks =========
    blocks.append(psf.common.MassFlowMergerBlock(name="merge_ox_recirc", 
                                                 st_in=["ox_engine_in", "ox_shaft_recirc"],
                                                 st_out="ox_pump_in",
                                                 medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.PumpBlock(name="ox_pump", 
                                       st_in="ox_pump_in", 
                                       st_out="ox_pump_out",
                                       overcome=["duct_pump_chamber_ox", "ox_injector"],
                                       load_fraction=1.0,
                                       p_base=params["p_c"],
                                       input_p = params["p_tank_ox"],
                                       n=params["n_ox"], eta=params["eta_pump_ox"],
                                       medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="split_ox_recirc", 
                                                   st_in="ox_pump_out",
                                                   st_out=["ox_duct_in", "ox_shaft_recirc"],
                                                   fractions=[1-params["zeta_pump_ox_recirc"], params["zeta_pump_ox_recirc"]],
                                                   medium=params["ox"]["coolprop"]))
    

    blocks.append(psf.common.SimpleDuctBlock(name="duct_pump_chamber_ox", 
                                             st_in="ox_duct_in", 
                                             st_out="ox_injector_plenum",
                                             eta=params["eta_pump_injector_ox"],
                                             medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="ox_injector", 
                                             st_in="ox_injector_plenum", 
                                             st_out="ox_chamber_in",
                                             eta=params["eta_ox_injector"],
                                             medium=params["ox"]["coolprop"]))
    

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

        # Propellants
        'fu': {"cea": "H2",
               "cantera": "H2", 
               "coolprop": "hydrogen"},
        
        "ox": {"cea": "LOX",
               "cantera": "O2", 
               "coolprop": "oxygen"},

        # Chamber geometry parameters
        'theta_conv': 25,
        'R_1f': 1.5,
        'R_2f': 3.0,
        'R_3f': 0.5,
        'length_fraction': 0.713,
        'n_fu': 31537,
        'n_ox': 12615,

        # Cooling channel parameters
        "wall_thickness": 0.3e-3,
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
    params['T_tank_fu'] = CP.PropsSI('T', 'P', params['p_tank_fu'], 'Q', 0, params["fu"]["coolprop"]) 

    
    params['rho_ox_tank'] = CP.PropsSI('D', 'P', params['p_tank_ox'], 'Q', 0, params["ox"]['coolprop'])
    params['rho_fu_tank'] = CP.PropsSI('D', 'P', params['p_tank_fu'], 'Q', 0, params["fu"]['coolprop'])

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
    out_path = os.path.join(script_dir, "engine_sizer_res.pkl")

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