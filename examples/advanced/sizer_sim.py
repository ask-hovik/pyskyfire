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
import numpy as np
import pyskyfire as psf
from scipy.integrate import cumulative_trapezoid


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
        AR_c=2.0,
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
        fuel_channel_count=140,
        nozzle_channel_count=120,
        channel_blockage_ratio=0.1,
        fuel_hot_gas_surface_area_multiplier=1.0,
        ox_outbound_hot_gas_surface_area_multiplier=1.2,
        ox_return_hot_gas_surface_area_multiplier=1.2,
        fuel_throat_channel_height=2.0e-3,
        fuel_channel_height_slope=6.0e-3,
        ox_channel_height=6e-3,
        ox_single_pass_channel_height=2.5e-3,

        # All cooling circuits use the same uncoated copper wall.
        copper_thickness=0.6e-3,

        # Materials
        copper=psf.common.solids.GRCop42,

        # Pump parameters
        eta_pump_fu=0.72,
        n_fu=50000,  # rpm
        eta_pump_ox_stage1=0.76,
        eta_pump_ox_stage2=0.72,
        n_ox=25000,  # rpm

        # Turbine efficiencies
        eta_turbine_fu=0.75,
        eta_turbine_ox=0.72,

        # Duct pressure ratios: fuel side
        eta_pump_regen_fu=0.95,
        eta_regen_turbine_fu=0.98,
        eta_turbine_injector_fu=0.94,
        eta_fu_injector=0.88,

        # Duct pressure ratios: oxidizer side
        eta_pump_regen_ox=0.98,
        eta_regen_turbine_ox=0.98,
        eta_turbine_injector_ox=0.88,
        eta_ox_injector=0.88,

        # Mass-flow leakage fractions
        zeta_fu_recirc=0.005,
        zeta_ox_recirc=0.005,
        ox_regen_flow_fraction=0.4,
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
        minimum_cea_temperature=305.0,
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

    copper_wall = psf.regen.Wall(
        material=params["copper"],
        thickness=params["copper_thickness"],
    )

    def fuel_channel_height(x):
        return (
            params["fuel_throat_channel_height"]
            + params["fuel_channel_height_slope"] * abs(x)
        )

    ox_outbound_span = (0.11, 1.0)
    ox_return_span = (1.0, 0.21)

    def span_fraction_to_x(span_fraction):
        if span_fraction >= 0.0:
            return span_fraction * float(contour.xs[-1])
        return span_fraction * -float(contour.xs[0])

    ox_outbound_bounds = tuple(map(span_fraction_to_x, ox_outbound_span))
    ox_return_bounds = tuple(map(span_fraction_to_x, ox_return_span))
    ox_overlap_start = max(min(ox_outbound_bounds), min(ox_return_bounds))

    def ox_channel_height(x):
        if x < ox_overlap_start:
            return params["ox_single_pass_channel_height"]
        return params["ox_channel_height"]

    cross_section = psf.regen.CrossSectionSquared(
        blockage_ratio=params["channel_blockage_ratio"],
    )
    fuel_transport = psf.skycea.CoolantTransport(params["coolprop_fu"])
    ox_transport = psf.skycea.CoolantTransport(params["coolprop_ox"])

    fuel_channel_placement = psf.regen.SurfacePlacement(n_channel_positions=params["fuel_channel_count"],)
    nozzle_channel_placement = psf.regen.SurfacePlacement(n_channel_positions=params["nozzle_channel_count"],)

    fuel_chamber_pass = psf.regen.CoolingCircuit(
        name="Fuel Chamber Pass",
        contour=contour,
        coolant_transport=fuel_transport,
        cross_section=cross_section,
        span=[0.1, -1.0],
        placement=fuel_channel_placement,
        walls=[copper_wall],
        roughness=params["copper_roughness_height"],
        channel_height=fuel_channel_height,
        hot_gas_surface_area_multiplier=params[
            "fuel_hot_gas_surface_area_multiplier"
        ],
    )

    ox_nozzle_pass_outbound = psf.regen.CoolingCircuit(
        name="Oxidizer Nozzle Pass Outbound",
        contour=contour,
        coolant_transport=ox_transport,
        cross_section=cross_section,
        span=ox_outbound_span,
        placement=nozzle_channel_placement,
        walls=[copper_wall],
        roughness=params["copper_roughness_height"],
        channel_height=ox_channel_height,
        hot_gas_surface_area_multiplier=params[
            "ox_outbound_hot_gas_surface_area_multiplier"
        ],
    )

    ox_nozzle_pass_return = psf.regen.CoolingCircuit(
        name="Oxidizer Nozzle Pass Return",
        contour=contour,
        coolant_transport=ox_transport,
        cross_section=cross_section,
        span=ox_return_span,
        placement=nozzle_channel_placement,
        walls=[copper_wall],
        roughness=params["copper_roughness_height"],
        channel_height=ox_channel_height,
        hot_gas_surface_area_multiplier=params[
            "ox_return_hot_gas_surface_area_multiplier"
        ],
    )

    return psf.regen.ThrustChamber(
        contour=contour,
        combustion_transport=aerothermodynamics,
        cooling_circuits=[fuel_chamber_pass, ox_nozzle_pass_outbound, ox_nozzle_pass_return],
        h_gas_corr=1.0,
        h_cold_corr=1.0,
        n_nodes=100,
    )


# tutorial:end:thrust-chamber

def gaussian_density_nodes(
    start,
    end,
    count,
    center=0.0,
    amplitude=2.0,
    standard_deviation=0.04,
):
    """Distribute nodes using a Gaussian spatial-density weighting."""
    if count < 3:
        raise ValueError("count must be at least 3")
    if not start < center < end:
        raise ValueError("center must lie strictly between start and end")
    if amplitude < 0.0 or standard_deviation <= 0.0:
        raise ValueError(
            "amplitude must be nonnegative and standard_deviation must be positive"
        )

    samples = np.linspace(start, end, 2001)
    density = 1.0 + amplitude * np.exp(
        -0.5 * ((samples - center) / standard_deviation) ** 2
    )
    cumulative_density = cumulative_trapezoid(density, samples, initial=0.0)
    center_cumulative_density = np.interp(center, samples, cumulative_density)

    interval_count = count - 1
    left_fraction = center_cumulative_density / cumulative_density[-1]
    left_intervals = int(
        np.clip(round(left_fraction * interval_count), 1, interval_count - 1)
    )
    right_intervals = interval_count - left_intervals

    left_targets = np.linspace(
        0.0,
        center_cumulative_density,
        left_intervals + 1,
    )
    right_targets = np.linspace(
        center_cumulative_density,
        cumulative_density[-1],
        right_intervals + 1,
    )
    left_nodes = np.interp(left_targets, cumulative_density, samples)
    right_nodes = np.interp(right_targets, cumulative_density, samples)
    return np.concatenate((left_nodes[:-1], [center], right_nodes[1:]))


def make_regen_node_spacings(thrust_chamber):
    """Return the deliberately coarse grids used by this example.

    The fuel chamber pass gets 15 throat-biased stations, including the throat
    itself. The two oxidizer nozzle passes use eight uniform stations each.
    """
    throat_circuit = thrust_chamber.cooling_circuits[0]
    x_start = float(np.min(throat_circuit.x_domain))
    x_end = float(np.max(throat_circuit.x_domain))
    if not x_start < 0.0 < x_end:
        raise ValueError("The throat circuit must span x=0")

    throat_nodes = gaussian_density_nodes(
        start=x_start,
        end=x_end,
        count=15,
        center=0.0,
        amplitude=2.0,
        standard_deviation=0.04,
    )
    throat_grid = {
        "wall": throat_nodes.copy(),
        "heat_flux": throat_nodes.copy(),
        "coolant": throat_nodes.copy(),
    }

    return {0: throat_grid, 1: 8, 2: 8}


def _solve_regen_circuit(
    thrust_chamber,
    circuit_index,
    nodes,
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
        nodes=nodes,
        circuit_index=circuit_index,
        boundary_conditions=boundary_conditions,
        solver="newton",
        output=True,
        heat_curvature_correction=False,
    )


# tutorial:start:regen-only
def run_regen_only(params, thrust_chamber):
    """Solve the fuel pass and the two serial oxidizer passes."""

    aerothermodynamics = thrust_chamber.combustion_transport
    mdot_fu = aerothermodynamics.mdot_fu
    ox_regen_flow_fraction = float(params["ox_regen_flow_fraction"])
    if not 0.0 < ox_regen_flow_fraction < 1.0:
        raise ValueError("ox_regen_flow_fraction must be strictly between 0 and 1")
    mdot_ox = aerothermodynamics.mdot_ox * ox_regen_flow_fraction

    cooling_data = {}
    node_spacings = make_regen_node_spacings(thrust_chamber)

    cooling_data["fuel_chamber"] = _solve_regen_circuit(
        thrust_chamber=thrust_chamber,
        circuit_index=0,
        nodes=node_spacings[0],
        T_coolant_in=params["T_regen_fu_in"],
        p_coolant_in=params["p_regen_fu_in"],
        mdot_coolant=mdot_fu,
    )

    T_in = params["T_regen_ox_in"]
    p_in = params["p_regen_ox_in"]
    oxidizer_passes = (
        ("oxidizer_nozzle_outbound", 1),
        ("oxidizer_nozzle_return", 2),
    )
    for result_name, circuit_index in oxidizer_passes:
        result = _solve_regen_circuit(
            thrust_chamber=thrust_chamber,
            circuit_index=circuit_index,
            nodes=node_spacings[circuit_index],
            T_coolant_in=T_in,
            p_coolant_in=p_in,
            mdot_coolant=mdot_ox,
        )
        cooling_data[result_name] = result
        T_in = result.T_stagnation[-1]
        p_in = result.p_stagnation[-1]

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
    stations["fu_regen_out"]            = psf.common.Station(params["p_c"] * 1.4, params["T_tank_fu"] + 250, mdot_fu_est*0.6,)
    stations["fu_turbine_inlet_split"]  = psf.common.Station(params["p_c"] * 1.3, params["T_tank_fu"] + 250, mdot_fu_est,)
    stations["fu_bypass_valve"]         = psf.common.Station(params["p_c"] * 1.3, params["T_tank_fu"] + 250, mdot_fu_est,)
    stations["fu_turbine_in"]           = psf.common.Station(params["p_c"] * 1.2, params["T_tank_fu"] + 250, mdot_fu_est,)
    stations["fu_turbine_out"]          = psf.common.Station(params["p_c"] * 1.1, params["T_tank_fu"] + 200, mdot_fu_est,)
    stations["fu_turbine_outlet_merge"] = psf.common.Station(params["p_c"] * 1.1, params["T_tank_fu"] + 200, mdot_fu_est,)
    stations["fu_injector_plenum"]      = psf.common.Station(params["p_c"] * 1.1, params["T_tank_fu"] + 200, mdot_fu_est,)
    stations["fu_chamber_in"]           = psf.common.Station(params["p_c"], params["T_tank_fu"] + 200, mdot_fu_est,)

    ox_regen_flow_fraction = float(params["ox_regen_flow_fraction"])
    if not 0.0 < ox_regen_flow_fraction < 1.0:
        raise ValueError("ox_regen_flow_fraction must be strictly between 0 and 1")

    mdot_ox_est = thrust_chamber.combustion_transport.mdot_ox
    mdot_ox_main_est = mdot_ox_est * 0.8
    mdot_ox_regen_est = mdot_ox_main_est * ox_regen_flow_fraction
    mdot_ox_direct_est = mdot_ox_main_est - mdot_ox_regen_est
    mdot_ox_bypass_est = mdot_ox_regen_est * params["zeta_ox_turbine_bypass"]
    mdot_ox_turbine_est = mdot_ox_regen_est - mdot_ox_bypass_est
    stations["ox_engine_in"]            = psf.common.Station(params["p_tank_ox"], params["T_tank_ox"], mdot_ox_est*0.8,)
    stations["ox_pump_in"]              = psf.common.Station(stations["ox_engine_in"].p, stations["ox_engine_in"].T, mdot_ox_est*0.8,)
    stations["ox_stage1_pump_out"]      = psf.common.Station(params["p_c"] * 1.3, params["T_tank_ox"] + 10, mdot_ox_est*0.8,)
    stations["ox_shaft_recirc"]         = psf.common.Station(params["p_c"] * 1.3, params["T_tank_ox"] + 10, 0.1,)
    stations["ox_main_flow"]            = psf.common.Station(params["p_c"] * 1.3, params["T_tank_ox"] + 10, mdot_ox_main_est,)
    stations["ox_direct_branch"]        = psf.common.Station(params["p_c"] * 1.3, params["T_tank_ox"] + 10, mdot_ox_direct_est,)
    stations["ox_stage2_pump_in"]       = psf.common.Station(params["p_c"] * 1.3, params["T_tank_ox"] + 10, mdot_ox_regen_est,)
    stations["ox_regen_duct_in"]        = psf.common.Station(params["p_c"] * 1.8, params["T_tank_ox"] + 15, mdot_ox_regen_est,)
    stations["ox_regen_in"]             = psf.common.Station(params["p_c"] * 1.7, params["T_tank_ox"] + 15, mdot=mdot_ox_regen_est,)
    stations["ox_regen_interstage"]     = psf.common.Station(params["p_c"] * 1.7, params["T_tank_ox"] + 130, mdot_ox_regen_est,)
    stations["ox_regen_out"]            = psf.common.Station(params["p_c"] * 1.7, params["T_tank_ox"] + 250, mdot_ox_regen_est,)
    stations["ox_turbine_inlet_split"]  = psf.common.Station(params["p_c"] * 1.6, params["T_tank_ox"] + 250, mdot_ox_regen_est,)
    stations["ox_bypass_valve"]         = psf.common.Station(params["p_c"] * 1.6, params["T_tank_ox"] + 250, mdot_ox_bypass_est,)
    stations["ox_turbine_in"]           = psf.common.Station(params["p_c"] * 1.6, params["T_tank_ox"] + 250, mdot_ox_turbine_est,)
    stations["ox_turbine_out"]          = psf.common.Station(params["p_c"] * 1.3, params["T_tank_ox"] + 200, mdot_ox_turbine_est,)
    stations["ox_heated_branch"]        = psf.common.Station(params["p_c"] * 1.3, params["T_tank_ox"] + 200, mdot_ox_regen_est,)
    stations["ox_main_flow_merge"]      = psf.common.Station(params["p_c"] * 1.3, params["T_tank_ox"] + 100, mdot_ox_main_est,)
    stations["ox_injector_plenum"]      = psf.common.Station(params["p_c"] * 1.1, params["T_tank_ox"] + 100, mdot_ox_main_est,)
    stations["ox_chamber_in"]           = psf.common.Station(params["p_c"], params["T_tank_ox"] + 100, mdot_ox_main_est,)

    return stations


# tutorial:end:initial-stations

# tutorial:start:initial-signals
def setup_initial_signals(params):
    """Create the scalar initial guesses for the full-cycle solver."""

    return {
        "p_c": params["p_c"],
        "P_fuel_turbine_required": 2.8e5,
        "P_ox_turbine_required": 2.0e5,
    }


# tutorial:end:initial-signals

def engine_sizer(params, thrust_chamber):
    """Build and converge the full engine-cycle network."""

    if not 0.0 < params["ox_regen_flow_fraction"] < 1.0:
        raise ValueError("ox_regen_flow_fraction must be strictly between 0 and 1")

    # tutorial:start:network-setup
    stations = setup_initial_stations(params, thrust_chamber)
    signals = setup_initial_signals(params)
    blocks = []
    regen_nodes = make_regen_node_spacings(thrust_chamber)

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
                "Regen Fuel Chamber Pass",
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
            name="Regen Fuel Chamber Pass",
            st_in="fu_regen_in",
            st_out="fu_regen_out",
            circuit_index=0,
            thrust_chamber=thrust_chamber,
            medium=fuel_medium,
            nodes=regen_nodes[0],
            post_process_nodes=regen_nodes[0],
            heat_curvature_correction=False,
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
            name="Ox Pump Stage 1",
            st_in="ox_pump_in",
            st_out="ox_stage1_pump_out",
            overcome=[
                "Duct Turbine-Injector Ox",
                "Ox Injector",
            ],
            load_fraction=1.0,
            p_base=params["p_c"],
            input_p=params["p_tank_ox"],
            n=params["n_ox"],
            eta=params["eta_pump_ox_stage1"],
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.MassFlowSplitterBlock(
            name="Split Ox Recirc",
            st_in="ox_stage1_pump_out",
            st_out=["ox_main_flow", "ox_shaft_recirc"],
            fractions=[
                1 - params["zeta_ox_recirc"],
                params["zeta_ox_recirc"],
            ],
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.MassFlowSplitterBlock(
            name="Split Ox Main Flow",
            st_in="ox_main_flow",
            st_out=["ox_stage2_pump_in", "ox_direct_branch"],
            fractions=[
                params["ox_regen_flow_fraction"],
                1 - params["ox_regen_flow_fraction"],
            ],
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.PumpBlock(
            name="Ox Pump Stage 2",
            st_in="ox_stage2_pump_in",
            st_out="ox_regen_duct_in",
            overcome=[
                "Duct Pump-Regen Ox",
                "Regen Ox Nozzle Outbound",
                "Regen Ox Nozzle Return",
                "Duct Regen-Turbine Ox",
                "Ox Turbine",
            ],
            load_fraction=1.0,
            # Stage 2 supplies only the additional pressure required by the
            # regenerative-cooling and turbine branch.
            p_base=0.0,
            input_p=0.0,
            n=params["n_ox"],
            eta=params["eta_pump_ox_stage2"],
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

    blocks.append(
        psf.common.RegenBlock(
            name="Regen Ox Nozzle Outbound",
            st_in="ox_regen_in",
            st_out="ox_regen_interstage",
            circuit_index=1,
            thrust_chamber=thrust_chamber,
            medium=ox_medium,
            nodes=regen_nodes[1],
            post_process_nodes=regen_nodes[1],
            heat_curvature_correction=False,
        )
    )
    blocks.append(
        psf.common.RegenBlock(
            name="Regen Ox Nozzle Return",
            st_in="ox_regen_interstage",
            st_out="ox_regen_out",
            circuit_index=2,
            thrust_chamber=thrust_chamber,
            medium=ox_medium,
            nodes=regen_nodes[2],
            post_process_nodes=regen_nodes[2],
            heat_curvature_correction=False,
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
            st_out="ox_heated_branch",
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.MassFlowMergerBlock(
            name="Merge Ox Main Flow",
            st_in=["ox_direct_branch", "ox_heated_branch"],
            st_out="ox_main_flow_merge",
            medium=ox_medium,
        )
    )
    blocks.append(
        psf.common.SimpleDuctBlock(
            name="Duct Turbine-Injector Ox",
            st_in="ox_main_flow_merge",
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
            sink_keys=["P_Ox Pump Stage 1", "P_Ox Pump Stage 2"],
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
    net.run_fixed_point(tol=1e-3, max_iter=200)

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
        "fuel_chamber": block_results["Regen Fuel Chamber Pass"],
        "oxidizer_nozzle_outbound": block_results["Regen Ox Nozzle Outbound"],
        "oxidizer_nozzle_return": block_results["Regen Ox Nozzle Return"],
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
