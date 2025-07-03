# -------------------------------------------------------------
# RL10_engine_sizer.py
# -------------------------------------------------------------
import numpy as np
import CoolProp.CoolProp as CP
import pyskyfire as psf
import os
import time

def setup_initial_stations(params):
    st = {}

    # ==== fuel side ====
    st["fu_engine_in"]              = psf.common.Station(p=params["p_fu_tank"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_pump_in"]                = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_pump_out"]               = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_seal_shaft"]             = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_duct_to_regen_in"]       = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_regen_in"]               = psf.common.Station(p=params["p_c"]*2, T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_regen_out"]              = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_duct_to_turbine_out"]    = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_turbine_in"]             = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_turbine_out"]            = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_bypass_valve"]           = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_turbine_housing_out"]    = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_injector_plenum"]        = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)
    st["fu_chamber_in"]             = psf.common.Station(p=params["p_c"], T=params["T_fu_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_fu)


    # ==== LOX side ====
    st["ox_engine_in"]              = psf.common.Station(p=params["p_ox_tank"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_pump_in"]                = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_pump_out"]               = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_seal_shaft"]             = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_duct_to_regen_in"]       = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_regen_in"]               = psf.common.Station(p=params["p_c"]*2, T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_regen_interstage"]       = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_regen_out"]              = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_duct_to_turbine_out"]    = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_turbine_in"]             = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_turbine_out"]            = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_bypass_valve"]           = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_turbine_housing_out"]    = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_injector_plenum"]        = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)
    st["ox_chamber_in"]             = psf.common.Station(p=params["p_c"], T=params["T_ox_tank"], mdot=params["thrust_chamber"].optimal_values.mdot_ox)

    return st

# ------- 2. helper: initialise scalar signals ------------------------------
def setup_initial_signals(params):
    sg = {}

    # boundary conditions
    sg["p_c"]               = params["p_c"]

    # bootstrap guesses
    sg["P_fu_pump"]         = 2.8e5
    sg["P_ox_pump"]         = 1.0e5
    sg["P_ox_required"]     = sg["P_ox_pump"]
    sg["P_fu_required"]      = sg["P_fu_pump"]

    return sg

def setup_thrust_chamber(params): 

    
    optimals = psf.skycea.OptimalValues(params["ox"]["cea"], params["fu"]["cea"], params["F"], params["MR"], params["p_c"], params["p_e"], params["L_star"])

    V_c = optimals.V_c_opt
    r_t = optimals.r_t_opt
    eps = optimals.eps_opt

    xs, rs = psf.regen.contour_2.get_contour_2(V_c = V_c, 
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

    contour = psf.regen.Contour(xs, rs, name = "Methane Sea Level Engine Contour")

    #psf.regen.contour.plot_theta_vs_epsilon()

    wall = psf.regen.Wall(material = psf.common.material.CopperC106, thickness = params["copper_thickness"]) 
    coating = psf.regen.Wall(material = psf.common.material.ZirconiumOxide, thickness = params["zirconia_thickness"])
    wall_group = psf.regen.WallGroup(walls=[coating, wall])

    fuel_channel_height_fn = psf.regen.make_channel_height_fn(
        contour=contour, 
        region_fractions=[-1.0, 1.0], 
        flat_heights= [params["fu_flat"]], 
        pinch_factors= [params["fu_pinch"]], 
        transition_widths=[]
    )

    ox_1_channel_height_fn = psf.regen.make_channel_height_fn(
        contour=contour, 
        region_fractions=[-1.0, 1.0], 
        flat_heights= [params["ox_1_flat"]], 
        pinch_factors= [params["ox_1_pinch"]], 
        transition_widths=[]
    )

    ox_2_channel_height_fn = psf.regen.make_channel_height_fn(
        contour=contour, 
        region_fractions=[-1.0, 1.0], 
        flat_heights= [params["ox_2_flat"]], 
        pinch_factors= [params["ox_2_pinch"]], 
        transition_widths=[]
    )


    cross_section = psf.regen.CrossSectionRounded()
    fuel_transport = psf.skycea.CoolantTransport(params["fu"]["coolprop"])
    ox_transport = psf.skycea.CoolantTransport(params["ox"]["coolprop"])

    fuel_pass = psf.regen.CoolingCircuit(name="Fuel Pass", 
                                        contour=contour, 
                                        coolant_transport=fuel_transport, 
                                        cross_section=cross_section, 
                                        span = [0.25, -1.0], 
                                        placement=psf.regen.SurfacePlacement(n_channel_positions=200), 
                                        channel_height=fuel_channel_height_fn)

    ox_pass_1 = psf.regen.CoolingCircuit(name="Ox Pass 1", 
                                        contour=contour, 
                                        coolant_transport=ox_transport, 
                                        cross_section=cross_section, 
                                        span = [0.26, 1.0], 
                                        placement=psf.regen.SurfacePlacement(n_channel_positions=100), 
                                        channel_height=ox_1_channel_height_fn)

    ox_pass_2 = psf.regen.CoolingCircuit(name="Ox Pass 2", 
                                        contour=contour, 
                                        coolant_transport=ox_transport, 
                                        cross_section=cross_section, 
                                        span = [1.0, 0.26], 
                                        placement=psf.regen.SurfacePlacement(n_channel_positions=100), 
                                        channel_height=ox_2_channel_height_fn)

    cooling_circuit_group = psf.regen.CoolingCircuitGroup(circuit_list=[fuel_pass, ox_pass_1, ox_pass_2])

    combustion_transport = psf.skycea.CombustionTransport(ox=params["ox"]["cantera"], 
                                                        fu=params["fu"]["cantera"],
                                                        p_c=params["p_c"],
                                                        p_e=params["p_e"],
                                                        MR=params["MR"],
                                                        T_fuel_in=111, 
                                                        T_ox_in=91, 
                                                        mode="hybrid")

    combustion_transport.set_essentials(mdot=optimals.mdot, mdot_fu=optimals.mdot_fu, mdot_ox=optimals.mdot_ox)


    thrust_chamber = psf.regen.ThrustChamber(contour=contour, 
                                            wall_group=wall_group,
                                            combustion_transport=combustion_transport,  
                                            cooling_circuit_group=cooling_circuit_group,
                                            optimal_values = optimals)
    
    return thrust_chamber


# ------- 4. build network ---------------------------------------------------
def engine_sizer(params):

    thrust_chamber = setup_thrust_chamber(params)
    params["thrust_chamber"] = thrust_chamber
    stations = setup_initial_stations(params)
    signals  = setup_initial_signals(params)
    blocks   = []

    
    # ======== fuel side blocks ========
    blocks.append(psf.common.MassFlowMergerBlock(name="fu_merge_tank_recirc", 
                                                 st_in=["fu_engine_in", "fu_seal_shaft"],
                                                 st_out="fu_pump_in",
                                                 medium=params["fu"]["coolprop"]))

    blocks.append(psf.common.PumpBlock(name="fu_pump", 
                                       st_in="fu_pump_in",
                                       st_out="fu_pump_out",
                                       overcome=["fu_duct_to_regen", "fu_regen", "fu_duct_to_turbine", "fu_turbine", "fu_duct_to_injector", "fu_injector"],
                                       load_fraction= 1.0,
                                       p_base = params["p_c"],
                                       input_p = params["p_fu_tank"],
                                       eta=params["eta_fu_pump"], 
                                       n = params["n_fu"], 
                                       medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="fu_split_regen_recirc", 
                                                   st_in="fu_pump_out",
                                                   st_out=["fu_duct_to_regen_in", "fu_seal_shaft"],
                                                   fractions=[1-params["zeta_fu_pump_recirc"], params["zeta_fu_pump_recirc"]],
                                                   medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="fu_duct_to_regen",  
                                             st_in="fu_duct_to_regen_in", 
                                             st_out="fu_regen_in", 
                                             eta=params["eta_fu_duct_to_regen"],
                                             medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.RegenBlock(name="fu_regen", 
                                        st_in="fu_regen_in", 
                                        st_out="fu_regen_out", 
                                        circuit_index=0, 
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="fu_duct_to_turbine", 
                                             st_in="fu_regen_out", 
                                             st_out="fu_duct_to_turbine_out", 
                                             eta=params["eta_fu_duct_to_turbine"],
                                             medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="fu_split_turbine_bypass", 
                                                   st_in="fu_duct_to_turbine_out",
                                                   st_out=["fu_turbine_in", "fu_bypass_valve"],
                                                   fractions=[1-params["zeta_fu_turbine_bypass"], params["zeta_fu_turbine_bypass"]],
                                                   medium=params["fu"]["coolprop"]))
    
    
    blocks.append(psf.common.TurbineBlock(name="fu_turbine",  
                                          st_in="fu_turbine_in", 
                                          st_out="fu_turbine_out", 
                                          P_req_key="P_fu_required",
                                          eta=params["eta_fu_turbine"], 
                                          medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.MassFlowMergerBlock(name="fu_merge_turbine_bypass", 
                                                 st_in=["fu_turbine_out", "fu_bypass_valve"],
                                                 st_out="fu_turbine_housing_out",
                                                 medium=params["fu"]["coolprop"]))

    blocks.append(psf.common.SimpleDuctBlock(name="fu_duct_to_injector", 
                                             st_in="fu_turbine_housing_out", 
                                             st_out="fu_injector_plenum", 
                                             eta=params["eta_fu_duct_to_injector"],
                                             medium=params["fu"]["coolprop"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="fu_injector", 
                                             st_in="fu_injector_plenum", 
                                             st_out="fu_chamber_in", 
                                             eta=params["eta_fu_injector"],
                                             medium=params["fu"]["coolprop"]))
    


    # ======== ox side blocks =========
    blocks.append(psf.common.MassFlowMergerBlock(name="ox_merge_tank_recirc", 
                                                 st_in=["ox_engine_in", "ox_seal_shaft"],
                                                 st_out="ox_pump_in",
                                                 medium=params["ox"]["coolprop"]))

    blocks.append(psf.common.PumpBlock(name="ox_pump", 
                                       st_in="ox_pump_in",
                                       st_out="ox_pump_out",
                                       overcome=["ox_duct_to_regen", "ox_regen_1", "ox_regen_2", "ox_duct_to_turbine", "ox_turbine", "ox_duct_to_injector", "ox_injector"],
                                       load_fraction= 1.0,
                                       p_base = params["p_c"],
                                       input_p = params["p_ox_tank"],
                                       eta=params["eta_ox_pump"], 
                                       n = params["n_ox"], 
                                       medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="ox_split_regen_recirc", 
                                                   st_in="ox_pump_out",
                                                   st_out=["ox_duct_to_regen_in", "ox_seal_shaft"],
                                                   fractions=[1-params["zeta_ox_pump_recirc"], params["zeta_ox_pump_recirc"]],
                                                   medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="ox_duct_to_regen",  
                                             st_in="ox_duct_to_regen_in", 
                                             st_out="ox_regen_in", 
                                             eta=params["eta_ox_duct_to_regen"],
                                             medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.RegenBlock(name="ox_regen_1", 
                                        st_in="ox_regen_in", 
                                        st_out="ox_regen_interstage", 
                                        circuit_index=1, 
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.RegenBlock(name="ox_regen_2", 
                                        st_in="ox_regen_interstage", 
                                        st_out="ox_regen_out", 
                                        circuit_index=2, 
                                        thrust_chamber=params["thrust_chamber"],
                                        medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="ox_duct_to_turbine", 
                                             st_in="ox_regen_out", 
                                             st_out="ox_duct_to_turbine_out", 
                                             eta=params["eta_ox_duct_to_turbine"],
                                             medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.MassFlowSplitterBlock(name="ox_split_turbine_bypass", 
                                                   st_in="ox_duct_to_turbine_out",
                                                   st_out=["ox_turbine_in", "ox_bypass_valve"],
                                                   fractions=[1-params["zeta_ox_turbine_bypass"], params["zeta_ox_turbine_bypass"]],
                                                   medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.TurbineBlock(name="ox_turbine",  
                                          st_in="ox_turbine_in", 
                                          st_out="ox_turbine_out", 
                                          P_req_key="P_ox_required",
                                          eta=params["eta_ox_turbine"], 
                                          medium=params["ox"]["coolprop"],
                                          ))
    
    blocks.append(psf.common.MassFlowMergerBlock(name="ox_merge_turbine_bypass", 
                                                 st_in=["ox_turbine_out", "ox_bypass_valve"],
                                                 st_out="ox_turbine_housing_out",
                                                 medium=params["ox"]["coolprop"]))

    blocks.append(psf.common.SimpleDuctBlock(name="ox_duct_to_injector", 
                                             st_in="ox_turbine_housing_out", 
                                             st_out="ox_injector_plenum", 
                                             eta=params["eta_ox_duct_to_injector"],
                                             medium=params["ox"]["coolprop"]))
    
    blocks.append(psf.common.SimpleDuctBlock(name="ox_injector", 
                                             st_in="ox_injector_plenum", 
                                             st_out="ox_chamber_in", 
                                             eta=params["eta_ox_injector"],
                                             medium=params["ox"]["coolprop"]))
    

    blocks.append(psf.common.TransmissionBlock(name="fu_shaft", sink_keys=["P_fu_pump"], out_key="P_fu_required"))
    blocks.append(psf.common.TransmissionBlock(name="ox_shaft", sink_keys=["P_ox_pump"], out_key="P_ox_required"))


    for blk in blocks:
        if hasattr(blk, "dp_key"):
            st_in  = blk.station_inputs[0]  if blk.station_inputs  else None
            st_out = blk.station_outputs[0] if blk.station_outputs else None

            if st_in in stations and st_out in stations:
                dp0 = max(stations[st_in].p - stations[st_out].p, 0.0)
            else:
                
                dp0 = 0.0                      # fall-back if we lack both stations

            signals.setdefault(blk.dp_key, dp0) # TODO: internalise this or something. Very hacky. 

    # Instantiate simulation object and run simulation
    net = psf.common.EngineNetwork(stations, signals, blocks)
    net.run_fixed_point(tol=1e-5, max_iter=100)

    return net, params


def main():
    params = {
        # Core engine parameters
        'p_c': 100e5,
        'MR': 3.4,
        'AR_c': 1.8,
        'F': 100e3,
        'p_e': 0.8e5,
        'p_fu_tank': 3e5,
        'p_ox_tank': 4e5,
        'L_star': 1.2,

        # Propellants
        'fu': {"cea": "CH4",
               "cantera": "CH4", 
               "coolprop": "methane"},
        
        "ox": {"cea": "LOX",
               "cantera": "O2", 
               "coolprop": "oxygen"},

        # Chamber geometry parameters
        'theta_conv': 35,
        'R_1f': 1.5,
        'R_2f': 2.0,
        'R_3f': 0.5,
        'length_fraction': 0.8,
        'n_fu': 31537,
        'n_ox': 12615,

        # Cooling channel parameters
        "copper_thickness": 0.5e-3,
        'zirconia_thickness': 0.05e-3,
        'fu_flat': 0.003,
        'fu_pinch': -0.5,
        'ox_1_flat': 0.0001,
        'ox_1_pinch': 0.0,
        'ox_2_flat' :0.0001,
        'ox_2_pinch':-5.0,    

        # Pump efficiencies
        'eta_fu_pump': 0.5810,
        'eta_ox_pump': 0.6422,

        # Turbine efficiencies
        'eta_fu_turbine': 0.7353,
        'eta_ox_turbine': 0.7353,

        # Duct efficiencies fuel side
        'eta_fu_duct_to_regen': 0.95,
        'eta_fu_duct_to_turbine': 0.95,
        'eta_fu_duct_to_injector': 0.95,
        'eta_fu_injector': 0.94,

        # Duct efficiencies ox side
        'eta_ox_duct_to_regen': 0.95,
        'eta_ox_duct_to_turbine': 0.95,
        'eta_ox_duct_to_injector': 0.95,
        'eta_ox_injector': 0.94,

        # Mass flow leakage fractions
        'zeta_fu_pump_recirc': 0.0070,
        'zeta_fu_turbine_bypass': 0.00427, 
        'zeta_ox_pump_recirc': 0.0005,
        'zeta_ox_turbine_bypass': 0.00427, 
        
    }

    # Additional derived parameters
    params['T_ox_tank'] = CP.PropsSI('T', 'P', params['p_ox_tank'], 'Q', 0, params["ox"]['coolprop']) -2
    params['T_fu_tank'] = CP.PropsSI('T', 'P', params['p_fu_tank'], 'Q', 0, params["fu"]["coolprop"]) -2

    
    params['rho_ox_tank'] = CP.PropsSI('D', 'P', params['p_ox_tank'], 'Q', 0, params["ox"]['coolprop'])
    params['rho_fu_tank'] = CP.PropsSI('D', 'P', params['p_fu_tank'], 'Q', 0, params["fu"]['coolprop'])

    # time and run simulation
    start_time = time.time()
    net, input_params = engine_sizer(params)
    end_time = time.time()
    duration = end_time - start_time
    print(f"Sim took {duration} s")

    # Save results
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(script_dir, "methane_engine_sizer_res.pkl")
    res = psf.common.Results()
    res.add(name="net", obj=net)
    res.add(name="input_params", obj=input_params)
    res.save(out_path)

if __name__ == "__main__":
    # Code here only runs when the file is executed directly
    main()