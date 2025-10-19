from __future__ import annotations
import CoolProp.CoolProp as CP
import numpy as np
from copy import deepcopy
from abc import ABC, abstractmethod
from typing import Dict, List

# internal imports
from pyskyfire.common.engine_network import Station
from pyskyfire.regen.solver import BoundaryConditions, steady_heating_analysis

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
    :class:`~pyskyfire.common.engine_network.Station`
        Thermodynamic state container used by the network.
    :class:`~pyskyfire.common.blocks.SignalBlock`
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
    ) -> dict[str, any]:
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
    ) -> dict[str, any]:
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
    """Centrifugal/positive-displacement pump model (lumped).

    The block raises inlet pressure to meet a target set by the load
    (sum of downstream pressure drops) and computes shaft power and outlet state.

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
        Rotational speed (metadata).
    eta : float
        Pump efficiency (0..1).
    medium : str
        CoolProp fluid string.

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

    def __init__(self, name, st_in, st_out, overcome, p_base, input_p, load_fraction, n, eta, medium):

        self.overcome      = overcome
        self.load_fraction = load_fraction
        self.p_base = p_base
        
        self.name = name
        self.n = n
        self.eta = eta

        self.input_p = input_p
        super().__init__(medium)

        self.st_in = st_in # TODO: why double up
        self.st_out = st_out

        # metadata
        self.station_inputs = [st_in]
        self.station_outputs = [st_out]
        self.signal_inputs = []
        self.signal_outputs  = [f"P_{name}"]


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
            ``{f"P_{name}": float}`` required pump power [W].

        Raises
        ------
        KeyError
            If required stations or signals are missing.

        Notes
        -----
        Temperature property calls are offset by ``1e-3 K`` to avoid
        saturation-line singularities in CoolProp.
        """

        st_i = stations[self.st_in]           # Station object from upstream
        p_in = st_i.p
        T_in = st_i.T
        mdot_in = st_i.mdot                      # now taken from the station!

        dp_total = sum(signals[f"dp_{blk_name}"] for blk_name in self.overcome) + self.p_base - self.input_p
        p_target = stations[self.st_in].p + self.load_fraction * dp_total
        #print(f"p_target: {p_target}")



        rho_in   = CP.PropsSI('Dmass', 'T', T_in-1e-3, 'P', p_in,  self.medium) # offset T by slight amount to avoid saturation line. 
        h_in     = CP.PropsSI('Hmass','T', T_in-1e-3, 'P', p_in,  self.medium)

        w_ideal  = (p_target - p_in) / rho_in          # J/kg
        dh_act   =  w_ideal / self.eta                 # J/kg      (η is 0–1)

        h_out    = h_in + dh_act
        T_out    = CP.PropsSI('T', 'P', p_target, 'Hmass', h_out, self.medium)

        P_pump   = mdot_in * dh_act                    # W

        stations_out = {self.st_out: Station(p_target, T_out, mdot_in)}
        signals_out  = {f"P_{self.name}": P_pump}

        return stations_out, signals_out

    
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
                 medium):

        self.name = name
        self.st_in = st_in
        self.st_out = st_out
        self.circuit_index = circuit_index
        self.thrust_chamber = thrust_chamber
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


        cooling_data = steady_heating_analysis(
                           self.thrust_chamber,
                           n_nodes        = 100,
                           circuit_index  = self.circuit_index,
                           boundary_conditions = bc,
                           solver         = "newton",
                           output         = False
                       )

        # downstream thermo state
        T_out = cooling_data["T_stagnation"][-1]
        p_out = cooling_data["p_stagnation"][-1]

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
    ) -> dict[str, any]:
        
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
        cooling_data = steady_heating_analysis(
            self.thrust_chamber,
            n_nodes        = 50,
            circuit_index  = self.circuit_index,
            boundary_conditions = bc,
            solver         = "newton",
            output         = False,
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
            If inlet mass flow is non-positive.
        """

        st_i   = stations[self.st_in]
        P_req  = signals[self.P_req_key]           # [W]

        mdot   = st_i.mdot
        if mdot <= 0.0:
            raise ValueError(f"{self.name}: mdot must be positive")

        w_req  = P_req / mdot                      # J kg⁻¹

        # thermodynamic properties at inlet
        c_p    = CP.PropsSI("Cpmass", "T", st_i.T,
                            "P", st_i.p, self.medium)
        c_v    = CP.PropsSI("Cvmass", "T", st_i.T,
                            "P", st_i.p, self.medium)
        gamma  = c_p / c_v

        # isentropic outlet temperature drop
        T_out  = st_i.T - w_req / (self.eta * c_p)
        # pressure ratio from ideal-gas isentropic relation
        TPR    = (T_out / st_i.T) ** (gamma / (gamma - 1.0))
        p_out  = st_i.p * TPR

        dp     = st_i.p - p_out                    # positive number

        st_o = Station(p = p_out,
                       T = T_out,
                       mdot = mdot)

        return {self.st_out: st_o}, {self.dp_key: dp}



# ---------------------------------------------------------------------------
#  SIMPLE DUCT
# ---------------------------------------------------------------------------
class SimpleDuctBlock(FluidBlock):
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

    def __init__(self, name, sink_keys, out_key="P_required"): # TODO: don't hardcode P_required
                
        self.station_inputs   = []
        self.station_outputs  = []
        self.signal_inputs    = list(sink_keys)
        self.signal_outputs   = [out_key]
        self.sink_keys = sink_keys
        self.out_key   = out_key
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
        return {}, {self.out_key: P}

