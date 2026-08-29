"""
Methane/oxygen dual-expander example.

Set RUN_MODE to:
    "regen_only" - solve the thrust-chamber cooling circuits using explicit
                   standalone inlet boundary conditions.
    "full_cycle" - solve the complete engine network.

Both modes save the same core result contract:
    mode, params, thrust_chamber, cooling_data

Full-cycle results additionally include:
    net, stations, signals, residuals, block_results
"""

import argparse
from pathlib import Path
import time

import CoolProp.CoolProp as CP
import pyskyfire as psf


# tutorial:start:run-mode
# Choose: "regen_only" or "full_cycle" 
RUN_MODE = "full_cycle"
#RUN_MODE = "regen_only" 
RESULTS_FILENAME = "results.pkl"
# tutorial:end:run-mode


# tutorial:start:engine-inputs
def make_params():
    """Return the editable input set for this example."""

    params = dict(
        # Core engine parameters
        p_c=100e5,
        p_e=0.8e5,
        MR=2.8,
        AR_c=2.5,
        F=100e3,

        # Fuel/oxidizer parameters
        cea_fu=psf.common.Fluid(type="fuel", propellants=["CH4"], fractions=[1.0]),
        cea_ox=psf.common.Fluid(type="oxidizer", propellants=["O2"], fractions=[1.0]),
        coolprop_fu=psf.common.Fluid(type="fuel", propellants=["methane"], fractions=[1.0],),
        coolprop_ox=psf.common.Fluid(type="oxidizer", propellants=["oxygen"], fractions=[1.0],),
        T_gas_fu_in=300,  # Adjusted to fall within NASA CEA tables
        T_gas_ox_in=300,

        # Propellant tanks
        p_tank_ox=5e5,
        p_tank_fu=5e5,

        # Chamber/nozzle parameters
        theta_conv=35,
        R_1f=0.5,
        R_2f=1.0,
        R_3f=0.5,
        length_fraction=1.0,  # 80 % nozzle
        L_star=1.1,

        # Cooling channels
        copper_roughness_height=0.0030e-3,
        inconel_roughness_height=0.0015e-3,
        throat_channel_count=140,
        nozzle_channel_count=120,
        chamber_channel_count=220,
        ox_channel_height=6e-3,

        # Wall thicknesses
        oxygen_copper_thickness=0.6e-3,
        fuel_copper_thickness=0.6e-3,
        barrier_thickness=0.05e-3,
        inconel_thickness=0.4e-3,

        # Materials
        copper=psf.common.solids.GRCop42,
        zirconia=psf.common.solids.ZirconiumOxide,
        inconel=psf.common.solids.Inconel625,

        # Pump parameters
        eta_pump_fu=0.6,
        n_fu=50000,  # rpm
        eta_pump_ox=0.6,
        n_ox=25000,  # rpm

        # Turbine efficiencies
        eta_turbine_fu=0.73,
        eta_turbine_ox=0.73,

        # Duct pressure ratios: fuel side
        eta_pump_regen_fu=0.95,
        eta_regen_turbine_fu=0.98,
        eta_turbine_injector_fu=0.94,
        eta_fu_injector=0.88,

        # Duct pressure ratios: oxidizer side
        eta_pump_regen_ox=0.95,
        eta_regen_turbine_ox=0.95,
        eta_turbine_injector_ox=0.88,
        eta_ox_injector=0.88,

        # Mass-flow leakage fractions
        zeta_fu_recirc=0.005,
        zeta_ox_recirc=0.005,
        zeta_fu_turbine_bypass=0.005,
        zeta_ox_turbine_bypass=0.005,

        # Standalone regeneration-analysis boundary conditions, used in regen_only mode. 
        T_regen_fu_in=111.0,
        p_regen_fu_in=150e5,
        T_regen_ox_in=93.0,
        p_regen_ox_in=150e5,
    )

    # Derived tank properties used by the full-cycle initial guesses.
    params["T_tank_ox"] = CP.PropsSI("T", "P", params["p_tank_ox"], "Q", 0, params["coolprop_ox"].propellants[0],)
    params["T_tank_fu"] = CP.PropsSI("T", "P", params["p_tank_fu"], "Q", 0, params["coolprop_fu"].propellants[0],)

    return params

# tutorial:end:engine-inputs

# tutorial:start:thrust-chamber
def setup_thrust_chamber(params):
    """Build the chamber geometry, walls, channels, and gas-side transport."""

    aerothermodynamics = psf.skycea.Aerothermodynamics.from_F_pe_Lstar(
        fu=params["cea_fu"],
        ox=params["cea_ox"],
        T_fu_in=params["T_gas_fu_in"],
        T_ox_in=params["T_gas_ox_in"],
        MR=params["MR"],
        p_c=params["p_c"],
        F=params["F"],
        p_e=params["p_e"],
        L_star=params["L_star"],
        p_amb=1e5,
    )

    xs, rs = psf.regen.contour.get_contour(
        V_c=aerothermodynamics.V_c,
        r_t=aerothermodynamics.r_t,
        area_ratio=aerothermodynamics.eps,
        AR_c=params["AR_c"],
        theta_conv=params["theta_conv"],
        nozzle="rao",
        R_1f=params["R_1f"],
        R_2f=params["R_2f"],
        R_3f=params["R_3f"],
        length_fraction=params["length_fraction"],
    )
    contour = psf.regen.Contour(xs, rs, name="Methane Engine")

    oxygen_copper_wall = psf.regen.Wall(
        material=params["copper"],
        thickness=params["oxygen_copper_thickness"],
    )
    fuel_copper_wall = psf.regen.Wall(
        material=params["copper"],
        thickness=params["fuel_copper_thickness"],
    )
    barrier_wall = psf.regen.Wall(
        material=params["zirconia"],
        thickness=params["barrier_thickness"],
    )
    inconel_wall = psf.regen.Wall(
        material=params["inconel"],
        thickness=params["inconel_thickness"],
    )

    channel_height_fn_1 = psf.regen.make_channel_height_fn(
        contour=contour,
        region_fractions=[-1.0, 0.3, 1.0],
        flat_heights=[0.0025, 0.004],
        pinch_factors=[-1.0, -4.0],
        transition_widths=[0.000001],
    )

    channel_height_fn_2 = psf.regen.make_channel_height_fn(
        contour=contour,
        region_fractions=[-1.0, 0.3, 1.0],
        flat_heights=[0.0032, 0.004],
        pinch_factors=[0.8, -4.0],
        transition_widths=[0.0000001],
    )

    def ox_channel_height(_x):
        return params["ox_channel_height"]

    cross_section = psf.regen.CrossSectionRounded()
    fuel_transport = psf.skycea.CoolantTransport(params["coolprop_fu"])
    ox_transport = psf.skycea.CoolantTransport(params["coolprop_ox"])

    chamber_channel_placement = psf.regen.SurfacePlacement(n_channel_positions=params["chamber_channel_count"],)
    throat_channel_placement = psf.regen.SurfacePlacement(n_channel_positions=params["throat_channel_count"],)
    nozzle_channel_placement = psf.regen.SurfacePlacement(n_channel_positions=params["nozzle_channel_count"],)

    fu_copper_pass = psf.regen.CoolingCircuit(
        name="Fuel Copper Pass",
        contour=contour,
        coolant_transport=fuel_transport,
        cross_section=cross_section,
        span=[-0.3, 0.29],
        placement=throat_channel_placement,
        walls=[barrier_wall, fuel_copper_wall],
        roughness=params["copper_roughness_height"],
        channel_height=channel_height_fn_1,
    )

    fu_inco_pass_cocurrent = psf.regen.CoolingCircuit(
        name="Fuel Inconel Pass Cocurrent",
        contour=contour,
        coolant_transport=fuel_transport,
        cross_section=cross_section,
        span=[0.31, 1.0],
        placement=nozzle_channel_placement,
        walls=[fuel_copper_wall],
        roughness=params["inconel_roughness_height"],
        channel_height=channel_height_fn_1,
    )

    fu_inco_pass_countercurrent = psf.regen.CoolingCircuit(
        name="Fuel Inconel Pass Countercurrent",
        contour=contour,
        coolant_transport=fuel_transport,
        cross_section=cross_section,
        span=[1.0, 0.31],
        placement=nozzle_channel_placement,
        walls=[fuel_copper_wall],
        roughness=params["inconel_roughness_height"],
        channel_height=channel_height_fn_2,
    )

    ox_copper_pass = psf.regen.CoolingCircuit(
        name="Oxidizer Copper Pass",
        contour=contour,
        coolant_transport=ox_transport,
        cross_section=cross_section,
        span=[-0.31, -1.0],
        placement=chamber_channel_placement,
        walls=[oxygen_copper_wall],
        roughness=params["copper_roughness_height"],
        channel_height=ox_channel_height,
    )

    return psf.regen.ThrustChamber(
        contour=contour,
        combustion_transport=aerothermodynamics,
        cooling_circuits=[fu_copper_pass, fu_inco_pass_cocurrent, fu_inco_pass_countercurrent, ox_copper_pass,],
        h_gas_corr=1.0,
        h_cold_corr=1.0,
        n_nodes=100,
    )


# tutorial:end:thrust-chamber

def _solve_regen_circuit(
    thrust_chamber,
    circuit_index,
    n_nodes,
    T_coolant_in,
    p_coolant_in,
    mdot_coolant,
):
    """Run one steady heating calculation and return its axial result."""

    boundary_conditions = psf.regen.BoundaryConditions(
        T_coolant_in=T_coolant_in,
        p_coolant_in=p_coolant_in,
        mdot_coolant=mdot_coolant,
    )

    return psf.regen.coupled_steady_heating_analysis(
        thrust_chamber,
        nodes=n_nodes,
        circuit_index=circuit_index,
        boundary_conditions=boundary_conditions,
        solver="newton",
        output=True,
    )


# tutorial:start:regen-only
def run_regen_only(params, thrust_chamber):
    """Solve the sequential fuel passes and independent oxidizer pass."""

    aerothermodynamics = thrust_chamber.combustion_transport
    mdot_fu = aerothermodynamics.mdot_fu
    mdot_ox = aerothermodynamics.mdot_ox

    cooling_data = {}

    fuel_passes = (
        ("fuel_throat", 0, 100),
        ("fuel_nozzle_cocurrent", 1, 40),
        ("fuel_nozzle_countercurrent", 2, 40),
    )

    T_in = params["T_regen_fu_in"]
    p_in = params["p_regen_fu_in"]

    for result_name, circuit_index, n_nodes in fuel_passes:
        result = _solve_regen_circuit(
            thrust_chamber=thrust_chamber,
            circuit_index=circuit_index,
            n_nodes=n_nodes,
            T_coolant_in=T_in,
            p_coolant_in=p_in,
            mdot_coolant=mdot_fu,
        )
        cooling_data[result_name] = result
        T_in = result.T_stagnation[-1]
        p_in = result.p_stagnation[-1]

    cooling_data["oxidizer_copper"] = _solve_regen_circuit(
        thrust_chamber=thrust_chamber,
        circuit_index=3,
        n_nodes=40,
        T_coolant_in=params["T_regen_ox_in"],
        p_coolant_in=params["p_regen_ox_in"],
        mdot_coolant=mdot_ox,
    )

    return cooling_data


# tutorial:end:regen-only

# tutorial:start:initial-stations
def setup_initial_stations(params, thrust_chamber):
    """Create the initial station guesses used by the full-cycle solver."""

    stations = {}
    mdot_fu_est = thrust_chamber.combustion_transport.mdot_fu
    stations["fu_engine_in"]            = psf.common.Station(params["p_tank_fu"], params["T_tank_fu"], mdot_fu_est*0.6,)
    stations["fu_pump_in"]              = psf.common.Station(stations["fu_engine_in"].p, stations["fu_engine_in"].T, mdot_fu_est*0.6,)
    stations["fu_pump_out"]             = psf.common.Station(params["p_c"] * 1.4, params["T_tank_fu"] + 10, mdot_fu_est*0.6,)
    stations["fu_shaft_recirc"]         = psf.common.Station(params["p_c"] * 1.4, params["T_tank_fu"] + 10, 0.1,)
    stations["fu_regen_duct_in"]        = psf.common.Station(params["p_c"] * 1.4, params["T_tank_fu"] + 10, mdot_fu_est*0.6,)
    stations["fu_regen_in"]             = psf.common.Station(params["p_c"] * 1.5, params["T_tank_fu"] + 10, mdot_fu_est*0.6,)
    stations["fu_regen_interstage_1"]   = psf.common.Station(params["p_c"] * 1.4, params["T_tank_fu"] + 100, mdot_fu_est*0.6,)
    stations["fu_regen_interstage_2"]   = psf.common.Station(params["p_c"] * 1.4, params["T_tank_fu"] + 180, mdot_fu_est*0.6,)
    stations["fu_regen_out"]            = psf.common.Station(params["p_c"] * 1.4, params["T_tank_fu"] + 250, mdot_fu_est*0.6,)
    stations["fu_turbine_inlet_split"]  = psf.common.Station(params["p_c"] * 1.3, params["T_tank_fu"] + 250, mdot_fu_est,)
    stations["fu_bypass_valve"]         = psf.common.Station(params["p_c"] * 1.3, params["T_tank_fu"] + 250, mdot_fu_est,)
    stations["fu_turbine_in"]           = psf.common.Station(params["p_c"] * 1.2, params["T_tank_fu"] + 250, mdot_fu_est,)
    stations["fu_turbine_out"]          = psf.common.Station(params["p_c"] * 1.1, params["T_tank_fu"] + 200, mdot_fu_est,)
    stations["fu_turbine_outlet_merge"] = psf.common.Station(params["p_c"] * 1.1, params["T_tank_fu"] + 200, mdot_fu_est,)
    stations["fu_injector_plenum"]      = psf.common.Station(params["p_c"] * 1.1, params["T_tank_fu"] + 200, mdot_fu_est,)
    stations["fu_chamber_in"]           = psf.common.Station(params["p_c"], params["T_tank_fu"] + 200, mdot_fu_est,)

    mdot_ox_est = thrust_chamber.combustion_transport.mdot_ox
    stations["ox_engine_in"]            = psf.common.Station(params["p_tank_ox"], params["T_tank_ox"], mdot_ox_est*0.8,)
    stations["ox_pump_in"]              = psf.common.Station(stations["ox_engine_in"].p, stations["ox_engine_in"].T, mdot_ox_est*0.8,)
    stations["ox_pump_out"]             = psf.common.Station(params["p_c"] * 1.6, params["T_tank_ox"] + 10, mdot_ox_est*0.8,)
    stations["ox_shaft_recirc"]         = psf.common.Station(params["p_c"] * 1.6, params["T_tank_ox"] + 10, 0.1,)
    stations["ox_regen_duct_in"]        = psf.common.Station(params["p_c"] * 1.6, params["T_tank_ox"] + 10, mdot_ox_est,)
    stations["ox_regen_in"]             = psf.common.Station(params["p_c"] * 1.5, params["T_tank_ox"] + 10, mdot=mdot_ox_est,)
    stations["ox_regen_out"]            = psf.common.Station(params["p_c"] * 1.5, params["T_tank_ox"] + 250, mdot_ox_est / 2,)
    stations["ox_turbine_inlet_split"]  = psf.common.Station(params["p_c"] * 1.5, params["T_tank_ox"] + 250, mdot_ox_est,)
    stations["ox_bypass_valve"]         = psf.common.Station(params["p_c"] * 1.3, params["T_tank_ox"] + 250, mdot_ox_est,)
    stations["ox_turbine_in"]           = psf.common.Station(params["p_c"] * 1.3, params["T_tank_ox"] + 250, mdot_ox_est,)
    stations["ox_turbine_out"]          = psf.common.Station(params["p_c"] * 1.2, params["T_tank_ox"] + 200, mdot_ox_est,)
    stations["ox_turbine_outlet_merge"] = psf.common.Station(params["p_c"] * 1.1, params["T_tank_ox"] + 200, mdot_ox_est,)
    stations["ox_injector_plenum"]      = psf.common.Station(params["p_c"] * 1.1, params["T_tank_ox"] + 200, mdot_ox_est,)
    stations["ox_chamber_in"]           = psf.common.Station(params["p_c"], params["T_tank_ox"] + 200, mdot_ox_est,)

    return stations


# tutorial:end:initial-stations

# tutorial:start:initial-signals
def setup_initial_signals(params):
    """Create the scalar initial guesses for the full-cycle solver."""

    return {
        "p_c": params["p_c"],
        "P_fuel_pump": 2.8e5,
        "P_ox_pump": 2.0e5,
        "P_fuel_turbine_required": 2.8e5,
        "P_ox_turbine_required": 2.0e5,
    }


# tutorial:end:initial-signals

def engine_sizer(params, thrust_chamber):
    """Build and converge the full engine-cycle network."""

    # tutorial:start:network-setup
    stations = setup_initial_stations(params, thrust_chamber)
    signals = setup_initial_signals(params)
    blocks = []

    fuel_medium = params["coolprop_fu"].propellants[0]
    ox_medium = params["coolprop_ox"].propellants[0]
    # tutorial:end:network-setup

    # tutorial:start:fuel-side-blocks
    # Fuel-side blocks
    blocks.append(
        psf.common.MassFlowMergerBlock(
            name="Fuel Inlet Merge",
            st_in=["fu_engine_in", "fu_shaft_recirc"],
            st_out="fu_pump_in",
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.PumpBlock(
            name="Fuel Pump",
            st_in="fu_pump_in",
            st_out="fu_pump_out",
            overcome=[
                "Duct Pump-Regen Fuel",
                "Regen Throat Pass",
                "Duct Regen-Turbine Fuel",
                "Fuel Turbine",
                "Duct Turbine-Injector Fuel",
                "Fu Injector",
            ],
            load_fraction=1.0,
            p_base=params["p_c"],
            input_p=params["p_tank_fu"],
            eta=params["eta_pump_fu"],
            n=params["n_fu"],
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.MassFlowSplitterBlock(
            name="Fuel Recirc. Split",
            st_in="fu_pump_out",
            st_out=["fu_regen_duct_in", "fu_shaft_recirc"],
            fractions=[
                1 - params["zeta_fu_recirc"],
                params["zeta_fu_recirc"],
            ],
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.SimpleDuctBlock(
            name="Duct Pump-Regen Fuel",
            st_in="fu_regen_duct_in",
            st_out="fu_regen_in",
            pressure_ratio=params["eta_pump_regen_fu"],
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.RegenBlock(
            name="Regen Throat Pass",
            st_in="fu_regen_in",
            st_out="fu_regen_interstage_1",
            circuit_index=0,
            thrust_chamber=thrust_chamber,
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.RegenBlock(
            name="Regen Cocurrent Nozzle Pass",
            st_in="fu_regen_interstage_1",
            st_out="fu_regen_interstage_2",
            circuit_index=1,
            thrust_chamber=thrust_chamber,
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.RegenBlock(
            name="Regen Countercurrent Nozzle Pass",
            st_in="fu_regen_interstage_2",
            st_out="fu_regen_out",
            circuit_index=2,
            thrust_chamber=thrust_chamber,
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.SimpleDuctBlock(
            name="Duct Regen-Turbine Fuel",
            st_in="fu_regen_out",
            st_out="fu_turbine_inlet_split",
            pressure_ratio=params["eta_regen_turbine_fu"],
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.MassFlowSplitterBlock(
            name="Fuel Split Turbine Bypass",
            st_in="fu_turbine_inlet_split",
            st_out=["fu_turbine_in", "fu_bypass_valve"],
            fractions=[
                1 - params["zeta_fu_turbine_bypass"],
                params["zeta_fu_turbine_bypass"],
            ],
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.TurbineBlock(
            name="Fuel Turbine",
            st_in="fu_turbine_in",
            st_out="fu_turbine_out",
            P_req_key="P_fuel_turbine_required",
            eta=params["eta_turbine_fu"],
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.MassFlowMergerBlock(
            name="Merge Turbine Bypass",
            st_in=["fu_turbine_out", "fu_bypass_valve"],
            st_out="fu_turbine_outlet_merge",
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.SimpleDuctBlock(
            name="Duct Turbine-Injector Fuel",
            st_in="fu_turbine_outlet_merge",
            st_out="fu_injector_plenum",
            pressure_ratio=params["eta_turbine_injector_fu"],
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.SimpleDuctBlock(
            name="Fu Injector",
            st_in="fu_injector_plenum",
            st_out="fu_chamber_in",
            pressure_ratio=params["eta_fu_injector"],
            medium=fuel_medium,
        )
    )
    blocks.append(
        psf.common.TransmissionBlock(
            name="Fuel Shaft",
            sink_keys=["P_Fuel Pump"],
            source_keys=["P_fuel_turbine_required"],
        )
    )

    # tutorial:end:fuel-side-blocks

    # tutorial:start:oxidizer-side-blocks
    # Oxidizer-side blocks
    blocks.append(
        psf.common.MassFlowMergerBlock(
            name="Merge Ox Recirc",
            st_in=["ox_engine_in", "ox_shaft_recirc"],
            st_out="ox_pump_in",
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.PumpBlock(
            name="Ox Pump",
            st_in="ox_pump_in",
            st_out="ox_pump_out",
            overcome=[
                "Duct Pump-Regen Ox",
                "Ox Regen",
                "Duct Regen-Turbine Ox",
                "Ox Turbine",
                "Duct Turbine-Injector Ox",
                "Ox Injector",
            ],
            load_fraction=1.0,
            p_base=params["p_c"],
            input_p=params["p_tank_ox"],
            n=params["n_ox"],
            eta=params["eta_pump_ox"],
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.MassFlowSplitterBlock(
            name="Split Ox Recirc",
            st_in="ox_pump_out",
            st_out=["ox_regen_duct_in", "ox_shaft_recirc"],
            fractions=[
                1 - params["zeta_ox_recirc"],
                params["zeta_ox_recirc"],
            ],
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.SimpleDuctBlock(
            name="Duct Pump-Regen Ox",
            st_in="ox_regen_duct_in",
            st_out="ox_regen_in",
            pressure_ratio=params["eta_pump_regen_ox"],
            medium=ox_medium,
        )
    )

    # The duplicated name is deliberate: these are identical parallel channels,
    # so the pump should see a single branch pressure drop.
    blocks.append(
        psf.common.RegenBlock(
            name="Ox Regen",
            st_in="ox_regen_in",
            st_out="ox_regen_out",
            circuit_index=3,
            thrust_chamber=thrust_chamber,
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.SimpleDuctBlock(
            name="Duct Regen-Turbine Ox",
            st_in="ox_regen_out",
            st_out="ox_turbine_inlet_split",
            pressure_ratio=params["eta_regen_turbine_ox"],
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.MassFlowSplitterBlock(
            name="Ox Split Turbine Bypass",
            st_in="ox_turbine_inlet_split",
            st_out=["ox_turbine_in", "ox_bypass_valve"],
            fractions=[
                1 - params["zeta_ox_turbine_bypass"],
                params["zeta_ox_turbine_bypass"],
            ],
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.TurbineBlock(
            name="Ox Turbine",
            st_in="ox_turbine_in",
            st_out="ox_turbine_out",
            P_req_key="P_ox_turbine_required",
            eta=params["eta_turbine_ox"],
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.MassFlowMergerBlock(
            name="Merge Ox Bypass",
            st_in=["ox_turbine_out", "ox_bypass_valve"],
            st_out="ox_turbine_outlet_merge",
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.SimpleDuctBlock(
            name="Duct Turbine-Injector Ox",
            st_in="ox_turbine_outlet_merge",
            st_out="ox_injector_plenum",
            pressure_ratio=params["eta_turbine_injector_ox"],
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.SimpleDuctBlock(
            name="Ox Injector",
            st_in="ox_injector_plenum",
            st_out="ox_chamber_in",
            pressure_ratio=params["eta_ox_injector"],
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.TransmissionBlock(
            name="Ox Shaft",
            sink_keys=["P_Ox Pump"],
            source_keys=["P_ox_turbine_required"],
        )
    )

    # tutorial:end:oxidizer-side-blocks

    # tutorial:start:pressure-drop-signals
    # Initialise pressure-drop signals from the station guesses. # TODO internalise this process
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

    # tutorial:end:pressure-drop-signals

    # tutorial:start:network-solve
    net = psf.common.EngineNetwork(stations, signals, blocks)
    net.run_fixed_point(tol=1e-2, max_iter=100)

    return {
        "net": net,
        "stations": net.stations,
        "signals": net.signals,
        "residuals": net.residuals,
        "block_results": net.block_results,
    }

    # tutorial:end:network-solve

# tutorial:start:full-cycle-cooling-data
def cooling_data_from_full_cycle(block_results):
    """Normalize full-cycle block results to the common cooling-data contract."""

    return {
        "fuel_throat": block_results["Regen Throat Pass"],
        "fuel_nozzle_cocurrent": block_results["Regen Cocurrent Nozzle Pass"],
        "fuel_nozzle_countercurrent": block_results["Regen Countercurrent Nozzle Pass"],
        "oxidizer_copper": block_results["Ox Regen"],
    }

# tutorial:end:full-cycle-cooling-data

def main(output_dir: Path | None = None):
    # tutorial:start:run-and-save
    if RUN_MODE not in {"regen_only", "full_cycle"}:
        raise ValueError(
            f"RUN_MODE must be 'regen_only' or 'full_cycle', not {RUN_MODE!r}."
        )

    params = make_params()
    thrust_chamber = setup_thrust_chamber(params)

    start_time = time.time()

    if RUN_MODE == "regen_only":
        cooling_data = run_regen_only(params, thrust_chamber)
        cycle_results = None
    else:
        cycle_results = engine_sizer(params, thrust_chamber)
        cooling_data = cooling_data_from_full_cycle(
            cycle_results["block_results"],
        )

    duration = time.time() - start_time
    print(f"{RUN_MODE} simulation completed in {duration:.2f} s")

    results = psf.common.Results()
    results.add(name="mode", obj=RUN_MODE)
    results.add(name="params", obj=params)
    results.add(name="thrust_chamber", obj=thrust_chamber)
    results.add(name="cooling_data", obj=cooling_data)

    if cycle_results is not None:
        for name, value in cycle_results.items():
            results.add(name=name, obj=value)

    if output_dir is None:
        output_dir = Path(__file__).resolve().parent
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / RESULTS_FILENAME
    results.save(output_path)
    print(f"Results saved to {output_path}")
    # tutorial:end:run-and-save


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for generated simulation results.",
    )

    args = parser.parse_args()
    main(output_dir=args.output_dir)
