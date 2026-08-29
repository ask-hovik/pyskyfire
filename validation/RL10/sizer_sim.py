"""RL10A-3-3A engine-cycle balance.

This script is the second half of the RL10 validation case. It reads the thrust
chamber written by ``regen_sim.py`` and closes the expander cycle around it:
two-stage fuel pump, cooling jacket, turbine, gearbox-driven oxidizer pump, and
the injectors.

The cycle is *not* the regenerative-cooling validation, and it does not run at
the same operating point. ``regen_sim.py`` sizes the chamber from Binder's
Tables 2.4.1/2.5.1 and drives the jacket with the measured inlet state (33 K,
69 bar), so its wall temperatures, heat fluxes, and pressures can be compared
directly with the published curves. The cycle here runs at Binder's *engine
station* balance instead, which is the one the paper's cycle section uses, and
solves for its own jacket inlet state. Its jacket results are still saved (the
``block_results`` entries for the two regen blocks are full ``RegenResult``
objects) but ``post_process.py`` deliberately does not plot them; the regen tabs
of the report show the validation run instead.

Run ``regen_sim.py`` first, then this script, then ``post_process.py``.
"""

import json
import os
import time

import pyskyfire as psf


script_dir = os.path.dirname(os.path.abspath(__file__))
reference_dir = os.path.join(script_dir, "reference_data")
REGEN_RESULTS_FILENAME = "regen_results.pkl"
SIZER_RESULTS_FILENAME = "sizer_results.pkl"
STATION_REFERENCE_FILENAME = "reference_station_data.json"


# =============
# Cycle inputs
# =============

def load_station_reference():
    """Return the Binder et al. engine station table."""
    path = os.path.join(reference_dir, STATION_REFERENCE_FILENAME)
    with open(path, encoding="utf-8") as reference_file:
        return json.load(reference_file)


def make_params(thrust_chamber):
    """Return the cycle inputs.

    Geometry and the gas-side heat load are inherited from ``thrust_chamber``.
    The mass flows and inlet states come from the engine station table in
    Binder et al., because that is the balance the paper's *cycle* section is
    written at, and it is the one the report compares against.

    Note that this is deliberately not the balance the thrust chamber was
    built on. ``regen_sim.py`` sizes the chamber from Tables 2.4.1 and 2.5.1
    (2.7093 kg/s fuel, 14.2369 kg/s ox, MR 5.255), which is the operating point
    the jacket reference curves are quoted at. The engine station data is a
    different balance (2.7587 and 13.9480 kg/s at the injector, MR 5.056). The
    cycle therefore transports the station-table flows through a chamber whose
    combustion state was solved 1.8% lean on fuel and at 3.8% higher MR. That
    inconsistency is the price of reusing one chamber for both studies; it
    shows up as a jacket heat pickup slightly different from the paper's.
    """
    reference = load_station_reference()
    fuel_reference = reference["engine_condition_fuel"]
    ox_reference = reference["engine_condition_oxidizer"]

    params = dict(
        # Chamber pressure is inherited; the chamber was sized for it.
        p_c=thrust_chamber.combustion_transport.p_c,

        # Injector-face mass flows, Binder et al. engine station data.
        mdot_fu=fuel_reference["fuel_injector_plenum"]["mass_flow"],
        mdot_ox=ox_reference["ox_injector_plenum"]["mass_flow"],

        # Coolant/propellant media for the network blocks.
        coolprop_fu="hydrogen",
        coolprop_ox="oxygen",

        # Engine inlet conditions, Binder et al. engine station data. The fuel
        # tank is subcooled by about 1.3 K, so its temperature is taken from
        # the table rather than assumed to sit on the saturation line.
        p_tank_fu=fuel_reference["fuel_tank"]["total_pressure"],
        T_tank_fu=fuel_reference["fuel_tank"]["total_temperature"],
        p_tank_ox=ox_reference["oxygen_tank"]["total_pressure"],
        T_tank_ox=ox_reference["oxygen_tank"]["total_temperature"],

        # Turbomachinery
        eta_pump_fu=0.5810,
        eta_pump_ox=0.6422,
        eta_turbine_fu=0.7353,
        n_fu=31537,   # rpm
        n_ox=12615,   # rpm
        stage1_load_fraction=0.5,

        # Duct pressure ratios, fuel side
        eta_pump_regen_fu=0.95,
        eta_regen_turbine_fu=0.98,
        eta_turbine_injector_fu=0.94,
        eta_fu_injector=0.88,

        # Duct pressure ratios, oxidizer side
        eta_pump_injector_ox=0.88,
        eta_ox_injector=0.88,

        # Mass-flow split fractions, read off the station table's successive
        # mass-flow steps: -0.698% across stage 1, -0.591% across stage 2
        # (split evenly between the shaft seal and the gearbox), -0.424% into
        # the turbine bypass, and -0.032% across the oxidizer pump.
        zeta_stage1_recirc=0.0070,
        zeta_stage2_recirc=0.00294,
        zeta_stage2_gearbox=0.00294,
        zeta_turbine_bypass=0.00427,
        zeta_pump_ox_recirc=0.000323,
    )

    # Engine inlet flows chosen so the *injector* sees the station-table flows.
    # The recirculation streams return to the pump inlet and cancel out, but the
    # gearbox stream is dumped overboard, so the fuel side must be fed slightly
    # rich of the injector requirement. Solving the merge/split balance:
    #   pump_in = engine_in + z1*pump_in + z2*(1 - z1)*pump_in
    #   jacket  = (1 - z1)*(1 - z2 - zg)*pump_in
    #
    # The station table cannot distinguish recirculation from overboard loss:
    # both give the same pump-inlet-to-injector ratio, 1/((1-z1)(1-z2-zg)) =
    # 1.01301, and differ only at the engine inlet, which the table does not
    # list. Feeding the injector 2.7587 kg/s therefore lands the pump inlet on
    # 2.7946 kg/s either way, matching the table to 0.004%.
    z1 = params["zeta_stage1_recirc"]
    z2 = params["zeta_stage2_recirc"]
    zg = params["zeta_stage2_gearbox"]
    params["mdot_fu_engine_in"] = params["mdot_fu"] * (
        (1.0 - z1 - z2 * (1.0 - z1)) / ((1.0 - z1) * (1.0 - z2 - zg))
    )
    # The oxidizer recirculation is fully internal, so no correction is needed.
    params["mdot_ox_engine_in"] = params["mdot_ox"]

    return params


# =================
# Initial guesses
# =================

def setup_initial_stations(params):
    """Seed every station in the network with a plausible state."""
    p_c = params["p_c"]
    mdot_fu = params["mdot_fu"]
    mdot_ox = params["mdot_ox"]
    T_fu = params["T_tank_fu"]
    T_ox = params["T_tank_ox"]

    st = {}

    # ==== LH2 side ====
    st["fu_engine_in"]        = psf.common.Station(params["p_tank_fu"], T_fu, params["mdot_fu_engine_in"])
    st["stage1_pump_in"]      = psf.common.Station(params["p_tank_fu"], T_fu, mdot_fu)
    st["stage1_shaft_recirc"] = psf.common.Station(p_c * 1.20, T_fu, mdot_fu * params["zeta_stage1_recirc"])
    st["stage2_shaft_recirc"] = psf.common.Station(p_c * 2.30, T_fu, mdot_fu * params["zeta_stage2_recirc"])
    st["pump_interstage1"]    = psf.common.Station(p_c * 1.20, T_fu + 5.0, mdot_fu)
    st["pump_interstage2"]    = psf.common.Station(p_c * 1.20, T_fu + 5.0, mdot_fu)
    st["stage2_pump_out"]     = psf.common.Station(p_c * 2.30, T_fu + 10.0, mdot_fu)
    st["gearbox_dump"]        = psf.common.Station(p_c * 2.30, T_fu + 10.0, mdot_fu * params["zeta_stage2_gearbox"])
    st["regen_duct_in"]       = psf.common.Station(p_c * 2.30, T_fu + 10.0, mdot_fu)
    st["regen_in"]            = psf.common.Station(p_c * 2.20, T_fu + 10.0, mdot_fu)
    st["regen_interstage"]    = psf.common.Station(p_c * 2.10, T_fu + 40.0, mdot_fu)
    st["regen_out"]           = psf.common.Station(p_c * 1.80, T_fu + 220.0, mdot_fu)
    st["bypass_in"]           = psf.common.Station(p_c * 1.70, T_fu + 220.0, mdot_fu)
    st["turbine_in"]          = psf.common.Station(p_c * 1.70, T_fu + 220.0, mdot_fu)
    st["bypass_valve"]        = psf.common.Station(p_c * 1.60, T_fu + 220.0, mdot_fu * params["zeta_turbine_bypass"])
    st["bypass_out"]          = psf.common.Station(p_c * 1.20, T_fu + 200.0, mdot_fu)
    st["turbine_out"]         = psf.common.Station(p_c * 1.20, T_fu + 200.0, mdot_fu)
    st["fu_injector_plenum"]  = psf.common.Station(p_c * 1.14, T_fu + 200.0, mdot_fu)
    st["fu_chamber_in"]       = psf.common.Station(p_c, T_fu + 200.0, mdot_fu)

    # ==== LOX side ====
    st["ox_engine_in"]        = psf.common.Station(params["p_tank_ox"], T_ox, params["mdot_ox_engine_in"])
    st["ox_pump_in"]          = psf.common.Station(params["p_tank_ox"], T_ox, mdot_ox)
    st["ox_pump_out"]         = psf.common.Station(p_c * 1.30, T_ox + 3.0, mdot_ox)
    st["ox_shaft_recirc"]     = psf.common.Station(p_c * 1.30, T_ox + 3.0, mdot_ox * params["zeta_pump_ox_recirc"])
    st["ox_duct_in"]          = psf.common.Station(p_c * 1.30, T_ox + 3.0, mdot_ox)
    st["ox_injector_plenum"]  = psf.common.Station(p_c * 1.14, T_ox + 3.0, mdot_ox)
    st["ox_chamber_in"]       = psf.common.Station(p_c, T_ox + 3.0, mdot_ox)

    return st


def setup_initial_signals(params):
    """Seed the scalar signals: chamber pressure and the shaft powers."""
    signals = {
        "p_c": params["p_c"],
        "P_stage1_fuel_pump": 4.0e5,
        "P_stage2_fuel_pump": 4.0e5,
        "P_ox_pump": 0.6e5,
    }
    signals["P_required"] = (
        signals["P_stage1_fuel_pump"]
        + signals["P_stage2_fuel_pump"]
        + signals["P_ox_pump"]
    )
    return signals


# =====================
# Network construction
# =====================

# Everything the fuel pump has to push through, from its own discharge to the
# injector face. Both stages carry half of it.
FUEL_PUMP_LOAD = [
    "duct_pump_regen",
    "regen_half_pass",
    "regen_full_pass",
    "duct_regen_turbine",
    "turbine",
    "duct_turbine_injector",
    "fu_injector",
]

OX_PUMP_LOAD = [
    "duct_pump_chamber_ox",
    "ox_injector",
]

def build_blocks(params, thrust_chamber):
    """Assemble the RL10 expander-cycle block list in solve order."""
    fu = params["coolprop_fu"]
    ox = params["coolprop_ox"]
    blocks = []

    # ======= Fuel side =======
    blocks.append(psf.common.MassFlowMergerBlock(
        name="merge_stage1_recirc",
        st_in=["fu_engine_in", "stage1_shaft_recirc", "stage2_shaft_recirc"],
        st_out="stage1_pump_in",
        medium=fu,
    ))

    blocks.append(psf.common.PumpBlock(
        name="stage1_fuel_pump",
        st_in="stage1_pump_in",
        st_out="pump_interstage1",
        overcome=FUEL_PUMP_LOAD,
        load_fraction=params["stage1_load_fraction"],
        p_base=params["p_c"],
        input_p=params["p_tank_fu"],
        eta=params["eta_pump_fu"],
        n=params["n_fu"],
        medium=fu,
    ))

    blocks.append(psf.common.MassFlowSplitterBlock(
        name="split_stage1_recirc",
        st_in="pump_interstage1",
        st_out=["pump_interstage2", "stage1_shaft_recirc"],
        fractions=[
            1.0 - params["zeta_stage1_recirc"],
            params["zeta_stage1_recirc"],
        ],
        medium=fu,
    ))

    blocks.append(psf.common.PumpBlock(
        name="stage2_fuel_pump",
        st_in="pump_interstage2",
        st_out="stage2_pump_out",
        overcome=FUEL_PUMP_LOAD,
        load_fraction=1.0 - params["stage1_load_fraction"],
        p_base=params["p_c"],
        input_p=params["p_tank_fu"],
        eta=params["eta_pump_fu"],
        n=params["n_fu"],
        medium=fu,
    ))

    blocks.append(psf.common.MassFlowSplitterBlock(
        name="split_stage2_recirc",
        st_in="stage2_pump_out",
        st_out=["regen_duct_in", "stage2_shaft_recirc", "gearbox_dump"],
        fractions=[
            1.0 - params["zeta_stage2_recirc"] - params["zeta_stage2_gearbox"],
            params["zeta_stage2_recirc"],
            params["zeta_stage2_gearbox"],
        ],
        medium=fu,
    ))

    blocks.append(psf.common.SimpleDuctBlock(
        name="duct_pump_regen",
        st_in="regen_duct_in",
        st_out="regen_in",
        pressure_ratio=params["eta_pump_regen_fu"],
        medium=fu,
    ))

    blocks.append(psf.common.RegenBlock(
        name="regen_half_pass",
        st_in="regen_in",
        st_out="regen_interstage",
        circuit_index=0,
        thrust_chamber=thrust_chamber,
        medium=fu,
    ))

    blocks.append(psf.common.RegenBlock(
        name="regen_full_pass",
        st_in="regen_interstage",
        st_out="regen_out",
        circuit_index=1,
        thrust_chamber=thrust_chamber,
        medium=fu,
    ))

    blocks.append(psf.common.SimpleDuctBlock(
        name="duct_regen_turbine",
        st_in="regen_out",
        st_out="bypass_in",
        pressure_ratio=params["eta_regen_turbine_fu"],
        medium=fu,
    ))

    blocks.append(psf.common.MassFlowSplitterBlock(
        name="split_turbine_bypass",
        st_in="bypass_in",
        st_out=["turbine_in", "bypass_valve"],
        fractions=[
            1.0 - params["zeta_turbine_bypass"],
            params["zeta_turbine_bypass"],
        ],
        medium=fu,
    ))

    blocks.append(psf.common.TurbineBlock(
        name="turbine",
        st_in="turbine_in",
        st_out="bypass_out",
        P_req_key="P_required",
        eta=params["eta_turbine_fu"],
        medium=fu,
    ))

    blocks.append(psf.common.MassFlowMergerBlock(
        name="merge_turbine_bypass",
        st_in=["bypass_out", "bypass_valve"],
        st_out="turbine_out",
        medium=fu,
    ))

    blocks.append(psf.common.SimpleDuctBlock(
        name="duct_turbine_injector",
        st_in="turbine_out",
        st_out="fu_injector_plenum",
        pressure_ratio=params["eta_turbine_injector_fu"],
        medium=fu,
    ))

    blocks.append(psf.common.SimpleDuctBlock(
        name="fu_injector",
        st_in="fu_injector_plenum",
        st_out="fu_chamber_in",
        pressure_ratio=params["eta_fu_injector"],
        medium=fu,
    ))

    # ======= Oxidizer side =======
    blocks.append(psf.common.MassFlowMergerBlock(
        name="merge_ox_recirc",
        st_in=["ox_engine_in", "ox_shaft_recirc"],
        st_out="ox_pump_in",
        medium=ox,
    ))

    blocks.append(psf.common.PumpBlock(
        name="ox_pump",
        st_in="ox_pump_in",
        st_out="ox_pump_out",
        overcome=OX_PUMP_LOAD,
        load_fraction=1.0,
        p_base=params["p_c"],
        input_p=params["p_tank_ox"],
        eta=params["eta_pump_ox"],
        n=params["n_ox"],
        medium=ox,
    ))

    blocks.append(psf.common.MassFlowSplitterBlock(
        name="split_ox_recirc",
        st_in="ox_pump_out",
        st_out=["ox_duct_in", "ox_shaft_recirc"],
        fractions=[
            1.0 - params["zeta_pump_ox_recirc"],
            params["zeta_pump_ox_recirc"],
        ],
        medium=ox,
    ))

    blocks.append(psf.common.SimpleDuctBlock(
        name="duct_pump_chamber_ox",
        st_in="ox_duct_in",
        st_out="ox_injector_plenum",
        pressure_ratio=params["eta_pump_injector_ox"],
        medium=ox,
    ))

    blocks.append(psf.common.SimpleDuctBlock(
        name="ox_injector",
        st_in="ox_injector_plenum",
        st_out="ox_chamber_in",
        pressure_ratio=params["eta_ox_injector"],
        medium=ox,
    ))

    # The single turbine drives both fuel-pump stages and, through the gearbox,
    # the oxidizer pump.
    blocks.append(psf.common.TransmissionBlock(
        name="transmission",
        sink_keys=["P_stage1_fuel_pump", "P_stage2_fuel_pump", "P_ox_pump"],
        source_keys=["P_required"],
    ))

    return blocks


def engine_sizer(params, thrust_chamber, *, tol=1e-5, max_iter=100):
    """Build and converge the RL10 cycle around the given thrust chamber."""
    stations = setup_initial_stations(params)
    signals = setup_initial_signals(params)
    blocks = build_blocks(params, thrust_chamber)

    # Seed each block's pressure-drop signal from the station guesses.
    # TODO: internalise this in EngineNetwork.
    for block in blocks:
        if not hasattr(block, "dp_key"):
            continue
        st_in = block.station_inputs[0] if block.station_inputs else None
        st_out = block.station_outputs[0] if block.station_outputs else None
        if st_in in stations and st_out in stations:
            dp0 = max(stations[st_in].p - stations[st_out].p, 0.0)
        else:
            dp0 = 0.0
        signals.setdefault(block.dp_key, dp0)

    net = psf.common.EngineNetwork(stations, signals, blocks)
    net.run_fixed_point(tol=tol, max_iter=max_iter)

    return net


# =====
# Main
# =====

def main(regen_path=None, out_path=None):
    if regen_path is None:
        regen_path = os.path.join(script_dir, REGEN_RESULTS_FILENAME)
    if out_path is None:
        out_path = os.path.join(script_dir, SIZER_RESULTS_FILENAME)

    if not os.path.exists(regen_path):
        raise FileNotFoundError(
            f"{regen_path} not found. Run regen_sim.py before sizer_sim.py."
        )

    regen_results = psf.common.Results.load(regen_path)
    thrust_chamber = regen_results.thrust_chamber

    params = make_params(thrust_chamber)

    start_time = time.time()
    net = engine_sizer(params, thrust_chamber)
    print(f"Cycle simulation completed in {time.time() - start_time:.1f} s")

    results = psf.common.Results()
    results.add(name="params", obj=params)
    results.add(name="regen_params", obj=regen_results.params)
    results.add(name="thrust_chamber", obj=thrust_chamber)
    results.add(name="net", obj=net)
    results.add(name="stations", obj=net.stations)
    results.add(name="signals", obj=net.signals)
    results.add(name="residuals", obj=net.residuals)
    results.add(name="block_results", obj=net.block_results)
    results.save(out_path)
    print(f"Results saved to {out_path}")

    return results


if __name__ == "__main__":
    main()
