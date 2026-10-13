from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, List, Any

import CoolProp.CoolProp as CP
import numpy as np

# internal imports
from pyskyfire.network.engine_network import Station
from pyskyfire.ducts.orifice import Orifice, OrificePlate
from pyskyfire.ducts.pipe import Pipe
from pyskyfire.ducts.volute import Volute
from pyskyfire.injector import Injector
from pyskyfire.pump.impeller import G0
from pyskyfire.thermo.state import FluidState, PHASE_UNKNOWN
from pyskyfire.regen.coupled_solver import (
    BoundaryConditions,
    RegenResult,
    coupled_steady_heating_analysis,
)


class FluidBlock(ABC):
    """Abstract base for blocks that transform **fluid** states.

    Each :class:`FluidBlock` consumes and/or produces :class:`Station`
    objects and may also consume or emit scalar *signals*
    (e.g., power, pressure drop).

    Parameters
    ----------
    medium : str
        CoolProp fluid identifier used for property calls
        (e.g., ``"N2O"`` or ``"Water"``).

    Attributes
    ----------
    station_inputs : list[str]
        Names of input station keys expected by the block.
    station_outputs : list[str]
        Names of output station keys produced by the block.
    signal_inputs : list[str]
        Names of scalar signal keys read by the block.
    signal_outputs : list[str]
        Names of scalar signal keys written by the block.
    medium : str
        CoolProp fluid string used by this instance.

    See Also
    --------
    :class:`~pyskyfire.network.engine_network.Station`
        Thermodynamic state container used by the network.
    :class:`~pyskyfire.network.blocks.SignalBlock`
        Abstract base for blocks that operate on **scalar signals** only.

    Notes
    -----
    .. note::
        Subclasses must implement :meth:`compute`. Missing declared inputs
        should be treated as an error; the engine network orchestrator is
        responsible for providing them.
    """

    station_inputs: list[str]
    station_outputs: list[str]
    signal_inputs: list[str]
    signal_outputs: list[str]

    def __init__(self, medium: str):
        super().__init__()
        self.medium = medium

    @abstractmethod
    def compute(
        self,
        stations_in : dict[str, Station],
        signals_in  : dict[str, float]
    ) -> tuple[dict[str, Station], dict[str, float]]:
        """Compute one pass of the block.

        Parameters
        ----------
        stations_in : dict[str, Station]
            Input network stations keyed by name.
        signals_in : dict[str, float]
            Input scalar signals keyed by name.

        Returns
        -------
        stations_out : dict[str, Station]
            Updated or newly created stations produced by the block.
        signals_out : dict[str, float]
            Updated or newly created scalar signals produced by the block.

        Notes
        -----
        .. Note::
            Implementations should treat missing inputs as an error. The engine network orchestrator is expected to provide the declared inputs.
        """
    
    def post_process(
        self,
        stations: dict[str, "Station"],
        signals : dict[str, float],
    ) -> Any:
        """Optional finalization hook run **after convergence**.

        Parameters
        ----------
        stations : dict[str, Station]
            Final converged stations.
        signals : dict[str, float]
            Final converged scalar signals.

        Returns
        -------
        dict[str, Any]
            Arbitrary post-processed results to be collected by
            the network (e.g., axial profiles, derived scalars). Default is empty.
        """
        return {}

class SignalBlock(ABC):
    """Abstract base for blocks that operate on **scalar signals** only.

    Attributes
    ----------
    station_inputs : list[str]
        Always empty for pure signal blocks.
    station_outputs : list[str]
        Always empty for pure signal blocks.
    signal_inputs : list[str]
        Names of scalar signal keys read by the block.
    signal_outputs : list[str]
        Names of scalar signal keys written by the block.

    See Also
    --------
    FluidBlock
        Base class for blocks that read/write fluid *stations*.
    """
        
    # metadata
    station_inputs : list[str];  station_outputs : list[str]
    signal_inputs  : list[str];  signal_outputs  : list[str]

    @abstractmethod
    def compute(
        self,
        stations_in : dict[str, Station],
        signals_in  : dict[str, float]
    ) -> tuple[dict[str, Station], dict[str, float]]:
        """Compute one pass of the block (signals only).

        Parameters
        ----------
        stations_in : dict[str, Station]
            Unused for pure signal blocks (present for interface uniformity).
        signals_in : dict[str, float]
            Input scalar signals keyed by name.

        Returns
        -------
        stations_out : dict[str, Station]
            Always empty for pure signal blocks.
        signals_out : dict[str, float]
            Signals produced by the block.
        """
    
    def post_process(
        self,
        stations: dict[str, "Station"],
        signals : dict[str, float],
    ) -> Any:
        """Optional finalization hook run **after convergence**.

        Parameters
        ----------
        stations : dict[str, Station]
            Final converged stations (unused).
        signals : dict[str, float]
            Final converged scalar signals (unused).

        Returns
        -------
        dict[str, Any]
            Arbitrary post-processed results. Default is empty.
        """
        return {}


class PumpBlock(FluidBlock):
    """Pump stage charged component by component.

    The block raises inlet pressure to meet a target set by the load
    (sum of downstream pressure drops) and computes shaft power and outlet state.

    Each loss is charged once, where it happens. The outlet is the impeller
    outlet, so a volute or diffuser after it is its own block and its loss is
    part of the load. The blades do ``g*H/eta_h`` of work on every kilogram,
    with ``hydraulic_efficiency`` the impeller's own. Disk friction on the
    outside of the wheel is added on top, and both end up as heat in the
    fluid. The mechanical efficiency then covers bearings and shaft seals,
    whose heat leaves through the housing rather than the flow. No leakage is
    modelled yet: every kilogram through the impeller is delivered.

    Parameters
    ----------
    name : str
        Block name used to form signal keys (e.g., ``P_{name}``).
    st_in : str
        Inlet station key.
    st_out : str
        Outlet station key.
    overcome : list[str]
        Names of blocks whose ``dp_*`` signals form the load to overcome.
    p_base : float
        Baseline pressure added to the load [Pa].
    input_p : float
        Upstream pressure offset already present [Pa].
    load_fraction : float
        Fraction of total load to apply in this pump (0..1).
    n : float
        Rotational speed [rpm].
    hydraulic_efficiency : float
        The impeller's hydraulic efficiency, eta_h,La (0..1].
    medium : str
        CoolProp fluid string.
    mechanical_efficiency : float, optional
        Bearings and shaft seals (0..1]. Default 1.
    disk_friction : pyskyfire.pump.DiskFriction, optional
        Friction on the wheel's outer faces, evaluated at the inlet state.

    Attributes
    ----------
    dp_key : str
        Name of emitted power signal (``P_{name}``).
    station_inputs, station_outputs, signal_inputs, signal_outputs : list[str]
        Network I/O metadata.

    See Also
    --------
    TurbineBlock
        Consumes power to satisfy aggregate shaft load.
    """

    def __init__(self, name, st_in, st_out, overcome, p_base, input_p, load_fraction, n,
                 hydraulic_efficiency, medium, *, mechanical_efficiency=1.0, disk_friction=None):
        for label, value in (("hydraulic_efficiency", hydraulic_efficiency),
                             ("mechanical_efficiency", mechanical_efficiency)):
            if not 0 < value <= 1:
                raise ValueError(f"{label} must lie in (0, 1]")

        self.overcome      = overcome
        self.load_fraction = load_fraction
        self.p_base = p_base
        
        self.name = name
        self.n = n
        self.hydraulic_efficiency = hydraulic_efficiency
        self.mechanical_efficiency = mechanical_efficiency
        self.disk_friction = disk_friction

        self.input_p = input_p
        super().__init__(medium)

        self.st_in = st_in # TODO: why double up
        self.st_out = st_out

        # metadata
        self.station_inputs = [st_in]
        self.station_outputs = [st_out]
        self.signal_inputs = []
        self.signal_outputs  = [f"P_{name}"]

    def _evaluate(self, stations, signals):
        st_i = stations[self.st_in]           # Station object from upstream
        p_in = st_i.p
        T_in = st_i.T
        mdot_in = st_i.mdot                      # now taken from the station!

        dp_total = sum(signals[f"dp_{blk_name}"] for blk_name in self.overcome) + self.p_base - self.input_p
        p_target = stations[self.st_in].p + self.load_fraction * dp_total

        # offset T by slight amount to avoid saturation line.
        rho_in   = CP.PropsSI('Dmass', 'T', T_in-1e-3, 'P', p_in,  self.medium)
        h_in     = CP.PropsSI('Hmass','T', T_in-1e-3, 'P', p_in,  self.medium)

        w_ideal  = (p_target - p_in) / rho_in           # J/kg, g*H
        w_blade  = w_ideal / self.hydraulic_efficiency  # J/kg, g*H_th
        P_blade  = mdot_in * w_blade
        P_disk   = 0.0
        if self.disk_friction is not None:
            mu_in = CP.PropsSI('V', 'T', T_in-1e-3, 'P', p_in, self.medium)
            P_disk = self.disk_friction.power(rho=rho_in, mu=mu_in, n=self.n)

        # Blade work and disk friction both stay in the fluid; the bearings'
        # share of the shaft power leaves through the housing.
        h_out    = h_in + (P_blade + P_disk) / mdot_in
        T_out    = CP.PropsSI('T', 'P', p_target, 'Hmass', h_out, self.medium)
        P_shaft  = (P_blade + P_disk) / self.mechanical_efficiency

        return Station(p_target, T_out, mdot_in), dict(
            shaft_power=P_shaft, blade_power=P_blade, disk_friction_power=P_disk,
            mechanical_loss=P_shaft - P_blade - P_disk, head=w_ideal / G0,
            theoretical_head=w_blade / G0, hydraulic_efficiency=self.hydraulic_efficiency,
            mechanical_efficiency=self.mechanical_efficiency, inlet_density=rho_in)

    def compute(self, stations, signals):
        """Compute outlet state and required pump power.

        Parameters
        ----------
        stations : dict[str, Station]
            Network stations (must contain ``st_in``).
        signals : dict[str, float]
            Scalar signals (must contain all ``dp_{blk}`` listed in ``overcome``).

        Returns
        -------
        stations_out : dict[str, Station]
            ``{st_out: Station}`` with updated pressure and temperature.
        signals_out : dict[str, float]
            ``{f"P_{name}": float}`` required shaft power [W].

        Raises
        ------
        KeyError
            If required stations or signals are missing.

        Notes
        -----
        Temperature property calls are offset by ``1e-3 K`` to avoid
        saturation-line singularities in CoolProp.
        """
        outlet, result = self._evaluate(stations, signals)
        return {self.st_out: outlet}, {f"P_{self.name}": result["shaft_power"]}

    def post_process(self, stations, signals):
        """The shaft power split by where it went, at the converged state."""
        return self._evaluate(stations, signals)[1]

    
class RegenBlock(FluidBlock):
    """Regenerative-cooling segment (single circuit).

    Runs a 1-D steady heating analysis for a specified cooling circuit and
    returns the outlet state and the pressure drop across the circuit.

    Parameters
    ----------
    name : str
        Block name; used to form ``dp_{name}``.
    st_in : str
        Inlet station key (coolant entering the circuit).
    st_out : str
        Outlet station key (coolant leaving the circuit).
    circuit_index : int
        Index of the cooling circuit to simulate.
    thrust_chamber
        Geometry/physics object passed to the solver.
    medium : str
        CoolProp fluid string for the coolant.
    nodes : int or explicit grids, optional
        Axial discretization used by the iterative network solve.
    post_process_nodes : int or explicit grids, optional
        Axial discretization used to collect the final profiles.
    heat_curvature_correction : bool, optional
        Whether to apply the gas-side heat-transfer curvature correction.
    pressure_curvature_correction : bool, optional
        Whether to apply the coolant pressure-drop curvature correction.
    friction_model : {'native', 'fluids'}, optional
        Backend supplying the coolant friction factor and its curvature
        multiplier. See :mod:`pyskyfire.common.friction`.

    Attributes
    ----------
    dp_key : str
        Name of the emitted pressure-drop signal.
    station_inputs, station_outputs, signal_inputs, signal_outputs : list[str]
        Network I/O metadata.
    """

    def __init__(self,
                 name: str,
                 st_in: str,
                 st_out: str,
                 circuit_index: int,
                 thrust_chamber, 
                 medium,
                 nodes=100,
                 post_process_nodes=50,
                 heat_curvature_correction=True,
                 pressure_curvature_correction=True,
                 friction_model="native"):

        self.name = name
        self.st_in = st_in
        self.st_out = st_out
        self.circuit_index = circuit_index
        self.thrust_chamber = thrust_chamber
        self.nodes = nodes
        self.post_process_nodes = post_process_nodes
        self.heat_curvature_correction = bool(heat_curvature_correction)
        self.pressure_curvature_correction = bool(pressure_curvature_correction)
        self.friction_model = friction_model
        super().__init__(medium)

        # ------------- metadata for EngineNetwork ---------------------
        self.station_inputs  = [st_in]
        self.station_outputs = [st_out]

        self.signal_inputs   = []                          # no scalars needed
        self.dp_key          = f"dp_{name}"
        self.signal_outputs  = [self.dp_key]

    # -----------------------------------------------------------------
    def compute(self,
                stations: dict[str, "Station"],
                signals : dict[str, float]
               ) -> tuple[dict[str, "Station"], dict[str, float]]:
        
        """Run steady heating analysis and emit outlet state and Δp.

        Parameters
        ----------
        stations : dict[str, Station]
            Network stations (must contain ``st_in``).
        signals : dict[str, float]
            Scalar signals (unused).

        Returns
        -------
        stations_out : dict[str, Station]
            ``{st_out: Station}`` at circuit exit.
        signals_out : dict[str, float]
            ``{dp_key: float}`` pressure loss across the circuit [Pa].
        """

        st_i = stations[self.st_in]           # Station object from upstream
        p_in = st_i.p
        T_in = st_i.T
        mdot = st_i.mdot                      # now taken from the station!

        # ---- call your PSF steady-state solver ----------------------
        #print(f"T_coolant_in: {st_i.T}")
        #print(f"p_coolant_in: {st_i.p}")
        #print(f"mdot_coolant_in: {st_i.mdot}")
        bc = BoundaryConditions(
                 T_coolant_in = T_in,
                 p_coolant_in = p_in,
                 mdot_coolant = mdot
             )


        cooling_data = coupled_steady_heating_analysis(
                           self.thrust_chamber,
                           nodes          = self.nodes,
                           circuit_index  = self.circuit_index,
                           boundary_conditions = bc,
                           solver         = "newton",
                           output         = False,
                           heat_curvature_correction = self.heat_curvature_correction,
                           pressure_curvature_correction = self.pressure_curvature_correction,
                           friction_model = self.friction_model,
                       )

        # downstream thermo state
        T_out = cooling_data.T_stagnation[-1]
        p_out = cooling_data.p_stagnation[-1]

        # pressure loss across the circuit (positive number)
        dp = p_in - p_out

        # ---- build return dicts -------------------------------------
        stations_out = {
            self.st_out: Station(p = p_out,
                                 T = T_out,
                                 mdot = mdot)
        }

        signals_out = {
            self.dp_key: dp          # <- goes to PressureSumBlock list
        }

        return stations_out, signals_out
    
    def post_process(
        self,
        stations: dict[str, "Station"],
        signals : dict[str, float],
    ) -> RegenResult:
        
        """Re-run the solver on a finer grid to collect detailed outputs.

        Parameters
        ----------
        stations : dict[str, Station]
            Final converged stations.
        signals : dict[str, float]
            Final converged scalars.

        Returns
        -------
        dict[str, Any]
            Solver output dictionary (profiles and scalars) suitable
            for reporting/plotting.
        """

        st_i = stations[self.st_in]

        bc = BoundaryConditions(
            T_coolant_in = st_i.T,
            p_coolant_in = st_i.p,
            mdot_coolant = st_i.mdot,
        )


        # Use a finer axial grid for the final report
        cooling_data = coupled_steady_heating_analysis(
            self.thrust_chamber,
            nodes          = self.post_process_nodes,
            circuit_index  = self.circuit_index,
            boundary_conditions = bc,
            solver         = "newton",
            output         = False,
            heat_curvature_correction = self.heat_curvature_correction,
            pressure_curvature_correction = self.pressure_curvature_correction,
            friction_model = self.friction_model,
        )

        # Handy scalar that users often want
        #cooling_data["dp"] = st_i.p - cooling_data["p_stagnation"][-1]

        return cooling_data          # any dict structure is fine

# ---------------------------------------------------------------------------
#  TURBINE
# ---------------------------------------------------------------------------
class TurbineBlock(FluidBlock):
    """Single-stage turbine model.

    Consumes a required shaft-power signal and computes an outlet state
    that provides that work at given efficiency.

    Parameters
    ----------
    name : str
        Block name.
    st_in : str
        Inlet station key.
    st_out : str
        Outlet station key.
    P_req_key : str
        Name of the required power signal [W].
    eta : float
        Turbine isentropic efficiency (0..1).
    medium : str
        CoolProp fluid string.

    Attributes
    ----------
    dp_key : str
        Name of the emitted pressure-drop signal.
    station_inputs, station_outputs, signal_inputs, signal_outputs : list[str]
        Network I/O metadata.

    See Also
    --------
    TransmissionBlock
        Aggregates power demands into the required shaft power.
    """

    def __init__(self,
                 name: str,
                 st_in: str,
                 st_out: str,
                 P_req_key: str,
                 eta: float,
                 medium: str):

        self.name       = name
        self.st_in      = st_in
        self.st_out     = st_out
        self.P_req_key  = P_req_key
        self.eta        = eta
        super().__init__(medium)

        # metadata -----------------------------------------------------
        self.station_inputs   = [st_in]
        self.station_outputs  = [st_out]
        self.signal_inputs    = [P_req_key]
        self.dp_key           = f"dp_{name}"
        self.signal_outputs   = [self.dp_key]

    # ----------------------------------------------------------------
    def compute(self,
                stations: Dict[str, "Station"],
                signals : Dict[str, float]
               ) -> tuple[Dict[str, "Station"], Dict[str, float]]:
        
        """Compute outlet state delivering the requested shaft power.

        Parameters
        ----------
        stations : dict[str, Station]
            Network stations (must contain ``st_in``).
        signals : dict[str, float]
            Scalar signals (must contain ``P_req_key``).

        Returns
        -------
        stations_out : dict[str, Station]
            ``{st_out: Station}`` with updated pressure and temperature.
        signals_out : dict[str, float]
            ``{dp_key: float}`` pressure drop across the turbine [Pa].

        Raises
        ------
        ValueError
            If inlet mass flow is non-positive, or if no physical expansion
            can supply the requested power.

        Notes
        -----
        Total-to-total expansion evaluated on enthalpy and entropy rather
        than on temperature and pressure, so no ideal-gas or constant-``c_p``
        assumption is made. The isentropic enthalpy drop ``w_req / η`` sets
        the exit pressure at constant inlet entropy; the actual drop
        ``w_req`` then sets the exit temperature at that pressure.
        """

        st_i   = stations[self.st_in]
        P_req  = signals[self.P_req_key]           # [W]

        mdot   = st_i.mdot
        if mdot <= 0.0:
            raise ValueError(f"{self.name}: mdot must be positive")

        w_req  = P_req / mdot                      # J kg⁻¹

        # inlet state (stagnation; the network carries total conditions)
        h_in   = CP.PropsSI("Hmass", "T", st_i.T,
                            "P", st_i.p, self.medium)
        s_in   = CP.PropsSI("Smass", "T", st_i.T,
                            "P", st_i.p, self.medium)

        # η relates the actual work to the isentropic work, so the isentropic
        # enthalpy drop is the larger of the two and is what sets the exit
        # pressure. The actual drop equals the extracted work and sets the
        # exit temperature.
        h_out_s = h_in - w_req / self.eta          # isentropic exit enthalpy
        h_out   = h_in - w_req                     # actual exit enthalpy

        try:
            p_out = CP.PropsSI("P", "Hmass", h_out_s,
                               "Smass", s_in, self.medium)
            # allow a hair of numerical overshoot so that w_req -> 0 (which
            # gives p_out == p_in) is a valid pass-through rather than an error
            ok = np.isfinite(p_out) and 0.0 < p_out <= st_i.p * (1.0 + 1e-9)
            p_out = min(p_out, st_i.p)
        except ValueError:
            ok = False

        if not ok:
            raise ValueError(
                f"{self.name}: cannot extract {P_req:.4g} W from {mdot:.4g} kg/s "
                f"of {self.medium} at {st_i.T:.4g} K / {st_i.p:.4g} Pa "
                f"(eta={self.eta:.3g}); no physical expansion satisfies the "
                f"required power"
            )

        T_out  = CP.PropsSI("T", "Hmass", h_out,
                            "P", p_out, self.medium)

        dp     = st_i.p - p_out                    # positive number

        st_o = Station(p = p_out,
                       T = T_out,
                       mdot = mdot)

        return {self.st_out: st_o}, {self.dp_key: dp}



# ---------------------------------------------------------------------------
#  SIMPLE DUCT
# ---------------------------------------------------------------------------
'''class SimpleDuctBlock(FluidBlock):
    """Homogeneous duct loss model with fixed efficiency.

    Applies a fixed pressure efficiency :math:`\\eta` via
    ``p_out = η * p_in`` while keeping temperature unchanged.

    Parameters
    ----------
    name : str
        Block name.
    st_in : str
        Inlet station key.
    st_out : str
        Outlet station key.
    eta : float
        Pressure efficiency (``0 < η ≤ 1``).
    medium : str
        CoolProp fluid string.

    Attributes
    ----------
    dp_key : str
        Name of the emitted pressure-drop signal.
    station_inputs, station_outputs, signal_inputs, signal_outputs : list[str]
        Network I/O metadata.

    Raises
    ------
    ValueError
        If ``eta`` is not in ``(0, 1]``.
    """

    def __init__(self,
                 name: str,
                 st_in: str,
                 st_out: str,
                 eta: float, 
                 medium):
        
        if not (0.0 < eta <= 1.0):
            raise ValueError("eta must be 0 < η ≤ 1")
        self.name  = name
        self.st_in = st_in
        self.st_out = st_out
        self.eta   = eta
        super().__init__(medium)

        # metadata -----------------------------------------------------
        self.station_inputs   = [st_in]
        self.station_outputs  = [st_out]
        self.signal_inputs    = []
        self.dp_key           = f"dp_{name}"
        self.signal_outputs   = [self.dp_key]

    # ----------------------------------------------------------------
    def compute(self,
                stations: Dict[str, "Station"],
                signals : Dict[str, float]
               ) -> tuple[Dict[str, "Station"], Dict[str, float]]:

        """Apply fixed loss and emit outlet station and Δp.

        Parameters
        ----------
        stations : dict[str, Station]
            Network stations (must contain ``st_in``).
        signals : dict[str, float]
            Scalar signals (unused).

        Returns
        -------
        stations_out : dict[str, Station]
            ``{st_out: Station}`` with ``p_out = η * p_in`` and unchanged ``T``.
        signals_out : dict[str, float]
            ``{dp_key: float}`` with ``dp = p_in - p_out`` [Pa].
        """

        st_i = stations[self.st_in]
        p_in, T_in, mdot = st_i.p, st_i.T, st_i.mdot

        p_out = self.eta * p_in
        dp    = p_in - p_out                     # Pa

        st_o = Station(p = p_out,
                       T = T_in,                 # no heat pick-up modelled
                       mdot = mdot)

        return {self.st_out: st_o}, {self.dp_key: dp}

from CoolProp.CoolProp import PropsSI'''

class SimpleDuctBlock(FluidBlock):
    """Adiabatic homogeneous duct loss model with fixed pressure ratio.

    Applies a fixed pressure ratio via

        p_out = pressure_ratio * p_in

    and computes the outlet temperature from constant specific enthalpy:

        h_out = h_in

    This is a better approximation for an adiabatic low-Mach duct or
    concentrated pressure-loss element than holding temperature constant.

    Notes
    -----
    This is not an isentropic expansion. Entropy increases. For real fluids,
    especially near saturation or the critical point, T_out may differ strongly
    from T_in.
    """

    def __init__(self,
                 name: str,
                 st_in: str,
                 st_out: str,
                 pressure_ratio: float,
                 medium: str):

        if not (0.0 < pressure_ratio <= 1.0):
            raise ValueError("pressure_ratio must be 0 < pressure_ratio ≤ 1")

        self.name = name
        self.st_in = st_in
        self.st_out = st_out
        self.pressure_ratio = pressure_ratio

        super().__init__(medium)

        self.station_inputs = [st_in]
        self.station_outputs = [st_out]
        self.signal_inputs = []

        self.dp_key = f"dp_{name}"
        self.signal_outputs = [self.dp_key]

    def compute(self,
                stations: Dict[str, "Station"],
                signals: Dict[str, float]
                ) -> tuple[Dict[str, "Station"], Dict[str, float]]:

        st_i = stations[self.st_in]

        p_in = st_i.p
        T_in = st_i.T
        mdot = st_i.mdot

        p_out = self.pressure_ratio * p_in
        dp = p_in - p_out

        # Inlet specific enthalpy [J/kg]
        h_in = CP.PropsSI("H", "P", p_in, "T", T_in, self.medium)

        # Adiabatic low-speed duct / throttling approximation:
        # h_out = h_in
        T_out = CP.PropsSI("T", "P", p_out, "H", h_in, self.medium)

        st_o = Station(
            p=p_out,
            T=T_out,
            mdot=mdot,
        )

        return {self.st_out: st_o}, {self.dp_key: dp}


class PipeBlock(FluidBlock):
    """Engine-network adapter for one explicitly resolved :class:`Pipe`.

    The contained pipe remains independently usable through :meth:`Pipe.solve`.
    This block only translates network stations into that public solver API and
    emits the resulting outlet station and ``dp_<name>`` loss signal.

    Parameters
    ----------
    name : str
        Block name used for the pressure-loss signal.
    st_in, st_out : str
        Inlet and outlet network station keys.
    pipe : pyskyfire.ducts.Pipe
        Physical pipe model, with or without a core.
    transport : object
        Fluid-property provider accepted by :meth:`Pipe.solve`.
    friction_model : {'native', 'fluids'}, optional
        Distributed-friction backend. Default ``'fluids'``.
    inlet_swirl : float, optional
        Swirl at the inlet mean radius [m/s], as :meth:`Pipe.solve` takes it.
        Network stations carry no swirl, so it is fixed here. Default 0.
    """

    def __init__(
        self,
        name: str,
        st_in: str,
        st_out: str,
        pipe: Pipe,
        transport,
        friction_model: str = "fluids",
        inlet_swirl: float = 0.0,
    ):
        if not isinstance(pipe, Pipe):
            raise TypeError("pipe must be a pyskyfire.ducts.Pipe")

        self.name = name
        self.st_in = st_in
        self.st_out = st_out
        self.pipe = pipe
        self.transport = transport
        self.friction_model = friction_model
        self.inlet_swirl = float(inlet_swirl)
        super().__init__(getattr(transport, "fluid", type(transport).__name__))

        self.station_inputs = [st_in]
        self.station_outputs = [st_out]
        self.signal_inputs = []
        self.dp_key = f"dp_{name}"
        self.signal_outputs = [self.dp_key]

    def _solve(self, station: Station):
        inlet_state = FluidState(
            T=station.T,
            p=station.p,
            h=float("nan"),
            quality=float("nan"),
            rho=float("nan"),
            phase=PHASE_UNKNOWN,
            T_sat=float("nan"),
        )
        return self.pipe.solve(
            inlet_state=inlet_state,
            mass_flow=station.mdot,
            transport=self.transport,
            inlet_swirl=self.inlet_swirl,
            friction_model=self.friction_model,
        )

    def compute(
        self,
        stations: Dict[str, "Station"],
        signals: Dict[str, float],
    ) -> tuple[Dict[str, "Station"], Dict[str, float]]:
        inlet = stations[self.st_in]
        result = self._solve(inlet)
        outlet = Station(
            p=float(result.p_stagnation[-1]),
            T=float(result.T_stagnation[-1]),
            mdot=inlet.mdot,
        )
        return {self.st_out: outlet}, {self.dp_key: result.total_pressure_loss}

    def post_process(self, stations, signals):
        """Return the detailed pipe result at the converged inlet condition."""
        return self._solve(stations[self.st_in])


class OrificeBlock(FluidBlock):
    """Engine-network adapter for one :class:`Orifice` or :class:`OrificePlate`.

    The contained restriction remains independently usable through its own
    ``solve``. This block only translates network stations into that public
    API and emits the outlet station and the ``dp_<name>`` loss signal. The
    signal is the irreversible loss: the whole drop for an orifice, and the
    permanent loss after recovery for a plate.

    Parameters
    ----------
    name : str
        Block name used for the pressure-loss signal.
    st_in, st_out : str
        Inlet and outlet network station keys.
    orifice : pyskyfire.ducts.Orifice or pyskyfire.ducts.OrificePlate
        Physical restriction.
    fluid : str or Fluid or EquationOfState
        The fluid, as :func:`~pyskyfire.thermo.eos.as_equation_of_state` reads it.
    viscosity : float, optional
        Passed to :meth:`OrificePlate.solve`; only a plate uses it.

    Notes
    -----
    The flow is adiabatic, so the outlet carries the inlet enthalpy.
    """

    def __init__(self, name: str, st_in: str, st_out: str, orifice, fluid,
                 viscosity: float | None = None):
        if not isinstance(orifice, (Orifice, OrificePlate)):
            raise TypeError("orifice must be a pyskyfire.ducts.Orifice or OrificePlate")
        if viscosity is not None and not isinstance(orifice, OrificePlate):
            raise ValueError("viscosity is only used by an OrificePlate")

        self.name = name
        self.st_in = st_in
        self.st_out = st_out
        self.orifice = orifice
        self.fluid = fluid
        self.viscosity = viscosity
        super().__init__(getattr(fluid, "fluid", fluid if isinstance(fluid, str) else repr(fluid)))

        self.station_inputs = [st_in]
        self.station_outputs = [st_out]
        self.signal_inputs = []
        self.dp_key = f"dp_{name}"
        self.signal_outputs = [self.dp_key]

    def _solve(self, station: Station):
        inlet_state = FluidState(
            T=station.T,
            p=station.p,
            h=float("nan"),
            quality=float("nan"),
            rho=float("nan"),
            phase=PHASE_UNKNOWN,
            T_sat=float("nan"),
        )
        kwargs = {} if self.viscosity is None else {"viscosity": self.viscosity}
        return self.orifice.solve(inlet_state=inlet_state, mass_flow=station.mdot,
                                  fluid=self.fluid, **kwargs)

    def compute(
        self,
        stations: Dict[str, "Station"],
        signals: Dict[str, float],
    ) -> tuple[Dict[str, "Station"], Dict[str, float]]:
        inlet = stations[self.st_in]
        result = self._solve(inlet)
        outlet = Station(
            p=float(result.outlet_state.p),
            T=float(result.outlet_state.T),
            mdot=inlet.mdot,
        )
        return {self.st_out: outlet}, {self.dp_key: result.pressure_loss}

    def post_process(self, stations, signals):
        """Return the detailed orifice result at the converged inlet condition."""
        return self._solve(stations[self.st_in])


class InjectorBlock(FluidBlock):
    """Engine-network adapter for one feed of an :class:`~pyskyfire.injector.Injector`.

    An injector with several feeds (a main fuel plenum, a film ring and an
    oxidizer plenum, say) takes one block per feed, each with its own plenum
    station. The block hands the plenum station and its mass flow to
    :meth:`Injector.solve_feed` and emits the face station and the
    ``dp_<name>`` pressure drop from plenum to face.

    Parameters
    ----------
    name : str
        Block name used for the pressure-loss signal.
    st_in, st_out : str
        Plenum and face station keys.
    injector : pyskyfire.injector.Injector
        The injector.
    feed : hashable
        Which of the injector's feeds this block carries.
    fluid : str or Fluid or EquationOfState
        The fluid, as :func:`~pyskyfire.thermo.eos.as_equation_of_state` reads it.

    Notes
    -----
    The flow is adiabatic, so the face station carries the plenum enthalpy.
    """

    def __init__(self, name: str, st_in: str, st_out: str, injector, feed, fluid):
        if not isinstance(injector, Injector):
            raise TypeError("injector must be a pyskyfire.injector.Injector")
        if feed not in injector.feeds:
            raise ValueError(f"{feed!r} is not a feed of this injector")

        self.name = name
        self.st_in = st_in
        self.st_out = st_out
        self.injector = injector
        self.feed = feed
        self.fluid = fluid
        super().__init__(getattr(fluid, "fluid", fluid if isinstance(fluid, str) else repr(fluid)))

        self.station_inputs = [st_in]
        self.station_outputs = [st_out]
        self.signal_inputs = []
        self.dp_key = f"dp_{name}"
        self.signal_outputs = [self.dp_key]

    def _solve(self, station: Station):
        inlet_state = FluidState(
            T=station.T,
            p=station.p,
            h=float("nan"),
            quality=float("nan"),
            rho=float("nan"),
            phase=PHASE_UNKNOWN,
            T_sat=float("nan"),
        )
        return self.injector.solve_feed(self.feed, inlet_state=inlet_state,
                                        mass_flow=station.mdot, fluid=self.fluid)

    def compute(
        self,
        stations: Dict[str, "Station"],
        signals: Dict[str, float],
    ) -> tuple[Dict[str, "Station"], Dict[str, float]]:
        inlet = stations[self.st_in]
        result = self._solve(inlet)
        face = result.orifice.outlet_state
        outlet = Station(p=float(face.p), T=float(face.T), mdot=inlet.mdot)
        return {self.st_out: outlet}, {self.dp_key: result.pressure_drop}

    def post_process(self, stations, signals):
        """Return the feed result at the converged plenum condition."""
        return self._solve(stations[self.st_in])

class VoluteBlock(FluidBlock):
    """Engine-network adapter for one :class:`Volute` or collector.

    The contained volute remains independently usable through
    :meth:`Volute.solve`. This block only translates network stations into that
    public solver API and emits the resulting outlet station and ``dp_<name>``
    loss signal. Its geometry is resolved once, when the volute is built, so
    the network may iterate the operating point freely.

    Parameters
    ----------
    name : str
        Block name used for the pressure-loss signal.
    st_in, st_out : str
        Inlet and outlet network station keys.
    volute : pyskyfire.ducts.Volute
        Physical volute model.
    transport : object
        Fluid-property provider accepted by :meth:`Volute.solve`.
    exit_diameter : float, optional
        Diameter of the round pipe the exit cone discharges to [m]. Defaults to
        the volute's attached :class:`~pyskyfire.ducts.ExitCone` final
        diameter, and is required when the volute has none.
    roughness : float, optional
        Absolute wall roughness [m].
    friction_model : {'native', 'fluids'}, optional
        Skin-friction backend. Default ``'fluids'``.

    Notes
    -----
    The volute is adiabatic, so the outlet station carries the inlet
    temperature. The inlet velocity components are a design input of the volute
    itself rather than something the network sets, so a mass flow that
    disagrees with them raises
    :class:`~pyskyfire.ducts.volute.VoluteConsistencyWarning` on every solve.
    """

    def __init__(
        self,
        name: str,
        st_in: str,
        st_out: str,
        volute: Volute,
        transport,
        exit_diameter: float | None = None,
        roughness: float = 0.0,
        friction_model: str = "fluids",
    ):
        if not isinstance(volute, Volute):
            raise TypeError("volute must be a pyskyfire.ducts.Volute")
        if exit_diameter is None and volute.exit_geometry is None:
            raise ValueError(
                "exit_diameter is required when the volute has no exit cone")

        self.name = name
        self.st_in = st_in
        self.st_out = st_out
        self.volute = volute
        self.transport = transport
        self.exit_diameter = None if exit_diameter is None else float(exit_diameter)
        self.roughness = float(roughness)
        self.friction_model = friction_model
        super().__init__(getattr(transport, "fluid", type(transport).__name__))

        self.station_inputs = [st_in]
        self.station_outputs = [st_out]
        self.signal_inputs = []
        self.dp_key = f"dp_{name}"
        self.signal_outputs = [self.dp_key]

    def _solve(self, station: Station):
        inlet_state = FluidState(
            T=station.T,
            p=station.p,
            h=float("nan"),
            quality=float("nan"),
            rho=float("nan"),
            phase=PHASE_UNKNOWN,
            T_sat=float("nan"),
        )
        return self.volute.solve(
            inlet_state=inlet_state,
            mass_flow=station.mdot,
            transport=self.transport,
            exit_diameter=self.exit_diameter,
            roughness=self.roughness,
            friction_model=self.friction_model,
        )

    def compute(
        self,
        stations: Dict[str, "Station"],
        signals: Dict[str, float],
    ) -> tuple[Dict[str, "Station"], Dict[str, float]]:
        inlet = stations[self.st_in]
        result = self._solve(inlet)
        outlet = Station(
            p=float(result.states[-1].p),
            T=inlet.T,
            mdot=inlet.mdot,
        )
        return {self.st_out: outlet}, {self.dp_key: result.total_pressure_loss}

    def post_process(self, stations, signals):
        """Return the detailed volute result at the converged inlet condition."""
        return self._solve(stations[self.st_in])


class DiffuserBlock(FluidBlock):
    """Engine-network adapter for one :class:`~pyskyfire.pump.VanedDiffuser`.

    Sits between a pump block, whose outlet is the impeller outlet, and the
    volute that collects the diffuser's discharge. It hands the network's
    inlet station to :meth:`VanedDiffuser.solve` and emits the outlet station
    and a ``dp_<name>`` loss signal. The geometry is resolved when the
    diffuser is designed, so the network may iterate the operating point
    freely. The diffuser is adiabatic; the outlet carries the inlet
    temperature, as a volute's does.

    Parameters
    ----------
    name : str
        Block name used for the pressure-loss signal.
    st_in, st_out : str
        Inlet and outlet network station keys.
    diffuser : pyskyfire.pump.VanedDiffuser
        The designed diffuser.
    transport : object
        Fluid-property provider with ``get_rho`` and ``get_mu``.
    """

    def __init__(self, name, st_in, st_out, diffuser, transport):
        from pyskyfire.pump.diffuser import VanedDiffuser
        if not isinstance(diffuser, VanedDiffuser):
            raise TypeError("diffuser must be a pyskyfire.pump.VanedDiffuser")
        self.name = name
        self.st_in = st_in
        self.st_out = st_out
        self.diffuser = diffuser
        self.transport = transport
        super().__init__(getattr(transport, "fluid", type(transport).__name__))
        self.station_inputs = [st_in]
        self.station_outputs = [st_out]
        self.signal_inputs = []
        self.dp_key = f"dp_{name}"
        self.signal_outputs = [self.dp_key]

    def _solve(self, station: Station):
        inlet_state = FluidState(
            T=station.T, p=station.p, h=float("nan"), quality=float("nan"),
            rho=float("nan"), phase=PHASE_UNKNOWN, T_sat=float("nan"),
        )
        return self.diffuser.solve(inlet_state=inlet_state, mass_flow=station.mdot,
                                   transport=self.transport)

    def compute(self, stations, signals):
        inlet = stations[self.st_in]
        result = self._solve(inlet)
        outlet = Station(p=float(result.states[-1].p), T=inlet.T, mdot=inlet.mdot)
        return {self.st_out: outlet}, {self.dp_key: result.total_pressure_loss}

    def post_process(self, stations, signals):
        """Return the detailed diffuser result at the converged inlet condition."""
        return self._solve(stations[self.st_in])


class MassFlowSplitterBlock(FluidBlock):
    """Split a single inlet stream into multiple outlet branches.

    Fractions can be fixed (sum to 1.0) or read each pass from scalar
    signals to enable dynamic splits.

    Parameters
    ----------
    name : str
        Block name.
    st_in : str
        Inlet station key.
    st_out : list[str]
        Outlet station keys.
    medium : str
        CoolProp fluid string.
    fractions : list[float] | None, optional
        Fixed fractions (sum to 1.0) if provided.
    frac_keys : list[str] | None, optional
        Signal keys providing fractions at runtime.

    Attributes
    ----------
    station_inputs, station_outputs, signal_inputs, signal_outputs : list[str]
        Network I/O metadata.

    Raises
    ------
    ValueError
        If both or neither of ``fractions`` and ``frac_keys`` are given,
        or if fixed ``fractions`` do not sum to 1.0.
    """

    def __init__(self,
                 name      : str,
                 st_in     : str,
                 st_out   : list[str],
                 medium,
                 fractions : list[float] | None = None,
                 frac_keys : list[str] | None = None, 
                 ):
        

        if (fractions is None) == (frac_keys is None):
            raise ValueError("Specify *either* fixed fractions *or* "
                             "signal keys, not both.")

        if fractions is not None:
            if not np.isclose(sum(fractions), 1.0, atol=1e-8):
                raise ValueError("fractions must sum to 1.0")
            self.fractions = fractions
            self.frac_keys = None
        else:
            self.fractions = None
            self.frac_keys = frac_keys

        self.name           = name
        self.st_in          = st_in
        self.st_out        = st_out
        super().__init__(medium)

        # ─── metadata ──────────────────────────────────────────────
        self.station_inputs  = [st_in]
        self.station_outputs = list(st_out)
        self.signal_inputs   = [] if frac_keys is None else list(frac_keys)
        self.signal_outputs  = []               # no new scalars emitted

    # --------------------------------------------------------------
    def compute(self, stations, signals):

        """Create one outlet station per branch based on split fractions.

        Parameters
        ----------
        stations : dict[str, Station]
            Network stations (must contain ``st_in``).
        signals : dict[str, float]
            Scalar signals supplying fractions when ``frac_keys`` is used.

        Returns
        -------
        stations_out : dict[str, Station]
            Output stations for each branch.
        signals_out : dict[str, float]
            Empty dict (no signals produced).

        Raises
        ------
        ValueError
            If dynamic fractions do not sum to 1.0.
        """

        st_i = stations[self.st_in]
        mdot = st_i.mdot

        # choose fixed or dynamic fractions
        if self.fractions is not None:
            fracs = self.fractions
        else:
            fracs = [signals[k] for k in self.frac_keys]
            if not np.isclose(sum(fracs), 1.0, atol=1e-6):
                raise ValueError(f"{self.name}: supplied fractions "
                                 "do not sum to 1.0")

        # build one Station per branch
        stations_out = {}
        for f, st_key in zip(fracs, self.st_out):
            stations_out[st_key] = Station(
                p    = st_i.p,         # no Δp inside the node itself
                T    = st_i.T,
                mdot = f * mdot
            )

        return stations_out, {}

class MassFlowMergerBlock(FluidBlock):
    """Merge multiple inlet streams into a single outlet.

    Assumes identical fluid species and mixes by enthalpy at a reference
    pressure equal to the minimum inlet pressure.

    Parameters
    ----------
    name : str
        Block name.
    st_in : list[str]
        Inlet station keys.
    st_out : str
        Outlet station key.
    medium : str
        CoolProp fluid string.

    Attributes
    ----------
    station_inputs, station_outputs, signal_inputs, signal_outputs : list[str]
        Network I/O metadata.
    """

    def __init__(self, name: str, st_in: list[str], st_out: str, medium):

        self.name   = name
        self.st_in = st_in
        self.st_out = st_out

        super().__init__(medium)

        self.station_inputs  = list(st_in)
        self.station_outputs = [st_out]
        self.signal_inputs   = []
        self.signal_outputs  = []

    # --------------------------------------------------------------
    def compute(self, stations, signals):

        """Mix inlet streams by mass and enthalpy to produce the outlet.

        Parameters
        ----------
        stations : dict[str, Station]
            Network stations (must contain all in ``st_in``).
        signals : dict[str, float]
            Scalar signals (unused).

        Returns
        -------
        stations_out : dict[str, Station]
            Outlet station dict.
        signals_out : dict[str, float]
            Empty dict.

        Raises
        ------
        ValueError
            If the merged total mass flow is non-positive.

        Notes
        -----
        .. Note::
            Property calls use a small temperature offset (``-1e-3 K``) to avoid saturation issues in CoolProp.
        """

        mdot_tot = 0.0
        h_sum    = 0.0
        p_ref    = min(stations[s].p for s in self.st_in)  # safe side
        

        for key in self.st_in:
            st = stations[key]
            mdot_tot += st.mdot

            h_i = CP.PropsSI("Hmass", "T", st.T-1e-3, "P", st.p, self.medium)  # ∗
            h_sum += st.mdot * h_i


        if mdot_tot <= 0:
            raise ValueError(f"{self.name}: merged mdot must be positive")

        h_mix = h_sum / mdot_tot
        T_mix = CP.PropsSI("T", "Hmass", h_mix, "P", p_ref, self.medium)  # ∗

        return {
            self.st_out: Station(p = p_ref, T = T_mix, mdot = mdot_tot)
        }, {}


# Signal blocks
class TransmissionBlock(SignalBlock):
    """Sum multiple shaft-power signals into one required-power signal.

    Useful for aggregating pump/compressor power demands prior to a turbine
    block that must supply the total shaft power.

    Parameters
    ----------
    name : str
        Block name.
    sink_keys : list[str]
        Input power signal keys to sum.
    out_key : str, optional
        Output key for the total required power (default ``"P_required"``).

    Attributes
    ----------
    station_inputs, station_outputs, signal_inputs, signal_outputs : list[str]
        Network I/O metadata.
    """

    def __init__(self, name, sink_keys, source_keys): # TODO: don't hardcode P_required
                
        self.station_inputs   = []
        self.station_outputs  = []
        self.signal_inputs    = list(sink_keys)
        self.signal_outputs   = [source_keys]
        self.sink_keys = sink_keys
        self.source_keys = source_keys
        self.name=name

    def compute(self, st, sg):
        """Sum input power signals.

        Parameters
        ----------
        st : dict[str, Station]
            Stations dict (unused).
        sg : dict[str, float]
            Signals dict containing all ``sink_keys`` [W].

        Returns
        -------
        stations_out : dict[str, Station]
            Empty dict.
        signals_out : dict[str, float]
            ``{out_key: total_power}`` in watts.
        """
                
        P = sum(sg[k] for k in self.sink_keys)
        return {}, {self.source_keys[0]: P}

