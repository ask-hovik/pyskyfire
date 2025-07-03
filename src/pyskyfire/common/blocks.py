from __future__ import annotations
import pyskyfire as psf
import CoolProp.CoolProp as CP
import numpy as np
from copy import deepcopy

g0 = 9.81

from abc import ABC, abstractmethod
from typing import Dict, List
from dataclasses import dataclass

class FluidBlock(ABC):
    # metadata
    station_inputs : list[str];  station_outputs : list[str]
    signal_inputs  : list[str];  signal_outputs  : list[str]

    def __init__(self, medium: str):
        """
        Base constructor for fluid-processing blocks.
        """
        super().__init__()
        self.medium = medium

    @abstractmethod
    def compute(
        self,
        stations_in : dict[str, psf.common.Station],
        signals_in  : dict[str, float]
    ) -> tuple[dict[str, psf.common.Station], dict[str, float]]:
        """
        Return two dicts:
            • updated/created Station objects
            • updated/created scalar signals
        """
    
    def post_process(
        self,
        stations: dict[str, "Station"],
        signals : dict[str, float],
    ) -> dict[str, any]:
        """Called once after convergence; override when you have extras."""
        return {}

class SignalBlock(ABC):
    # metadata
    station_inputs : list[str];  station_outputs : list[str]
    signal_inputs  : list[str];  signal_outputs  : list[str]

    @abstractmethod
    def compute(
        self,
        stations_in : dict[str, psf.common.Station],
        signals_in  : dict[str, float]
    ) -> tuple[dict[str, psf.common.Station], dict[str, float]]:
        """
        Return two dicts:
            • updated/created Station objects
            • updated/created scalar signals
        """
    
    def post_process(
        self,
        stations: dict[str, "Station"],
        signals : dict[str, float],
    ) -> dict[str, any]:
        """Called once after convergence; override when you have extras."""
        return {}


class PumpBlock(FluidBlock):
    """
    Encapsulates pump sizing and performance logic.
    """

    def __init__(self, name, st_in, st_out, overcome, p_base, input_p, load_fraction, n, eta, medium):
        """
        params is a dictionary containing keys like:
          - n_fu (pump rotational speed)
          - eta_pump_fu (pump efficiency)
          - rho_fu_tank (liquid density)
          ...
        """
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
        """
        Given pump inlet/outlet conditions, compute the required power, 
        outlet temperature, mass, etc.

        Returns a dict with keys:
          'P_pump' (pump power, W)
          'T_pump_out' (pump outlet temperature, K)

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

        stations_out = {self.st_out: psf.common.Station(p_target, T_out, mdot_in)}
        signals_out  = {f"P_{self.name}": P_pump}

        return stations_out, signals_out

    
class RegenBlock(FluidBlock):
    """
    Regenerative-cooling segment (single circuit).

    Parameters
    ----------
    name : str
        Unique tag; used to form the scalar Δp key.
    st_in : str
        Station key for coolant entering the circuit.
    st_out : str
        Station key written back to the network (exit of circuit).
    circuit_index : int
        Which cooling circuit of the engine model this block represents.
    engine : psf.Engine
        Your pre-built rocket engine object (geometry + prop data).
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

        st_i = stations[self.st_in]           # Station object from upstream
        p_in = st_i.p
        T_in = st_i.T
        mdot = st_i.mdot                      # now taken from the station!

        # ---- call your PSF steady-state solver ----------------------
        #print(f"T_coolant_in: {st_i.T}")
        #print(f"p_coolant_in: {st_i.p}")
        #print(f"mdot_coolant_in: {st_i.mdot}")
        bc = psf.regen.BoundaryConditions(
                 T_coolant_in = T_in,
                 p_coolant_in = p_in,
                 mdot_coolant = mdot
             )


        cooling_data = psf.regen.steady_heating_analysis(
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
            self.st_out: psf.common.Station(p = p_out,
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

        st_i = stations[self.st_in]

        bc = psf.regen.BoundaryConditions(
            T_coolant_in = st_i.T,
            p_coolant_in = st_i.p,
            mdot_coolant = st_i.mdot,
        )


        # Use a finer axial grid for the final report
        cooling_data = psf.regen.steady_heating_analysis(
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
    """
    Single-stage impulse/expander turbine.
    Produces an outlet Station and its pressure drop; consumes the
    shaft-power demand that the TransmissionBlock summed.
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

        st_o = psf.common.Station(p = p_out,
                       T = T_out,
                       mdot = mdot)

        return {self.st_out: st_o}, {self.dp_key: dp}



# ---------------------------------------------------------------------------
#  SIMPLE DUCT
# ---------------------------------------------------------------------------
class SimpleDuctBlock(FluidBlock):
    """
    Applies a fixed efficiency η to simulate a homogeneous duct loss:
        p_out = η · p_in
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

        st_i = stations[self.st_in]
        p_in, T_in, mdot = st_i.p, st_i.T, st_i.mdot

        p_out = self.eta * p_in
        dp    = p_in - p_out                     # Pa

        st_o = psf.common.Station(p = p_out,
                       T = T_in,                 # no heat pick-up modelled
                       mdot = mdot)

        return {self.st_out: st_o}, {self.dp_key: dp}


class FractionalMassFlowLossBlock(FluidBlock):
    """
    Peel off a fixed fraction of mdot from st_in and inject it into st_out
    *without* showing the underlying split / merge to normal users.

    Visible metadata:
        station_inputs  = [st_in, st_out]
        station_outputs = [st_in, st_out]
    Internally:
        ┌─ MassFlowSplitterBlock ─┐
        │  st_in  ──► main + bleed│
        └─────────────────────────┘
                         │
                         ▼
        ┌─ MassFlowMergerBlock  ─┐
        │  st_out + bleed ──► out│
        └─────────────────────────┘
    """

    def __init__(self, *, name: str,
                 st_in: str, st_out: str,
                 fraction: float, medium: str):

        if not 0.0 < fraction < 1.0:
            raise ValueError("fraction must be between 0 and 1 (exclusive)")

        self.name      = name
        self.st_in     = st_in
        self.st_out    = st_out
        self.fraction  = fraction
        super().__init__(medium)

        # hidden station keys
        self._st_main   = f"{st_in}__remain__{name}"
        self._st_bleed  = f"{st_out}__bleed__{name}"
        self._st_merged = f"{st_out}__merged__{name}"

        # hidden blocks
        self._splitter = MassFlowSplitterBlock(
            name      = f"{name}__split",
            st_in     = st_in,
            st_out    = [self._st_main, self._st_bleed],
            fractions = [1.0 - fraction, fraction],
            medium = medium
        )

        self._merger   = MassFlowMergerBlock(
            name   = f"{name}__merge",
            st_in  = [st_out, self._st_bleed],
            st_out = self._st_merged,
            medium = medium,
        )

        # what the *network* sees
        self.station_inputs   = [st_in, st_out]
        self.station_outputs  = [st_in, st_out]   # same keys as the user
        self.signal_inputs    = []
        self.signal_outputs   = []

    # --------------------------------------------------------------
    def compute(self, stations, signals):

        # 1. split the inlet stream
        d_st, _  = self._splitter.compute(stations, signals)

        # 2. overwrite original inlet key with the remaining flow
        d_st[self.st_in] = d_st.pop(self._st_main)

        # 3. merge bleed into st_out
        d_st2, _ = self._merger.compute({**stations, **d_st}, signals)

        # 4. overwrite original outlet key with merged result
        d_st[self.st_out] = d_st2[self._st_merged]

        return d_st, {}



# not updated bellow

# ──────────────────────────────────────────────────────────────────
# ---------------------------------------------------------------------------
#  ORIFICE PLATE
# ---------------------------------------------------------------------------
class OrificePlateBlock(FluidBlock):
    r"""
    Thin-plate orifice with a *permanent* pressure loss.

    The usual incompressible relation is applied

        ṁ = C_d · A · √(2 ρ Δp)   ⟹   Δp = (ṁ / (C_d A))² / (2 ρ)

    Parameters
    ----------
    name        : str
        Unique tag – shows up in the network’s scalar list as ``dp_<name>``.
    st_in, st_out : str
        Keys of the inlet and outlet *Station* objects.
    Cd          : float
        Discharge coefficient (0 < Cₙ ≤ 1).  Typical sharp-edged orifice ≈ 0.6.
    A           : float
        Flow area of the hole [m²].
    medium      : str
        CoolProp fluid identifier – must match the other stations.
    """

    def __init__(self, *, name: str,
                 st_in: str, st_out: str,
                 Cd: float, A: float,
                 medium: str):

        if Cd <= 0.0:
            raise ValueError("Cd must be > 0")
        if A  <= 0.0:
            raise ValueError("A must be > 0")

        self.name   = name
        self.st_in  = st_in
        self.st_out = st_out
        self.Cd     = Cd
        self.A      = A

        super().__init__(medium)

        # ─── metadata for EngineNetwork ────────────────────────────
        self.station_inputs   = [st_in]
        self.station_outputs  = [st_out]
        self.signal_inputs    = []
        self.dp_key           = f"dp_{name}"
        self.signal_outputs   = [self.dp_key]

    # ------------------------------------------------------------------
    def compute(self, stations: dict[str, "psf.common.Station"],
                signals : dict[str, float]
               ) -> tuple[dict[str, "psf.common.Station"], dict[str, float]]:

        st_i = stations[self.st_in]
        p_in, T_in, mdot = st_i.p, st_i.T, st_i.mdot

        if mdot <= 0.0:
            raise ValueError(f"{self.name}: mdot must be positive")

        # density at the *upstream* stagnation state
        rho = CP.PropsSI("Dmass", "T", T_in - 1.0e-3, "P", p_in, self.medium)  # small T offset keeps CoolProp away from the sat. line

        # incompressible orifice formula
        dp = (mdot / (self.Cd * self.A)) ** 2 / (2.0 * rho)         # Pa
        p_out = p_in - dp
        if p_out <= 0.0:
            raise RuntimeError(f"{self.name}: computed p_out ≤ 0 Pa – geometry/Cd/mdot inconsistent?")

        st_o = psf.common.Station(p=p_out, T=T_in, mdot=mdot)

        return {self.st_out: st_o}, {self.dp_key: dp}


class JunctionBlock(FluidBlock):
    """A generalised mass‑flow junction.

    The block can *merge* an arbitrary number of inlet streams and *split*
    the combined flow into an arbitrary number of outlet streams **without**
    hard‑coded mass‑flow ratios.  Instead, a small one‑dimensional root‑finder
    allocates the branch mass‑flows such that user‑supplied pressure
    constraints (or pressure equalisation between the branches) are met.

    Parameters
    ----------
    name
        A human‑readable unique tag – used only for bookkeeping.
    st_in
        List of input *Station* keys that feed the junction.
    st_out
        List of output *Station* keys produced by the junction.  Each element
        must also appear as a key in *equilibrate*.
    medium
        CoolProp fluid name – *all* inlets are assumed to be of identical
        composition.
    equilibrate
        Defines how the total mass‑flow is distributed among the outlets.  It
        is a list of one‑element dictionaries so that *order* is preserved::

            [
                {"fu_leg": ["duct_fu", "turbine"]},
                {"by_leg": "diff"}
            ]

        • If the *value* is ``"diff"`` the branch simply receives whatever is
          left after the other branches have been solved.
        • Otherwise, the value must be a *list of block names* that constitute
          the hydraulic path of that branch **in the EngineNetwork’s block
          order**.  During *compute* the junction varies the branch mass‑flow
          until the exit‑pressure of the final block matches the *current*
          network pressure of its outlet station.

        If *two* branches meet again downstream (→ equal outlet station), pass
        two dictionaries, both with a *list* of blocks.  The solver will then
        iterate to make the *exit pressures equal* between the two branches.
    block_lookup
        Mapping *block‑name → block‑instance* so that the junction can call
        ``block.compute`` internally while it searches for a pressure match.
        The simplest way is to build ``{blk.name: blk for blk in blocks}`` once
        when assembling the network and pass it to every *JunctionBlock*.
    tol
        Relative tolerance used by the in‑house bisection routine.
    """

    # ------------------------------------------------------------------
    def __init__(self, *, name: str, st_in: Sequence[str], st_out: Sequence[str],
                 medium: str, equilibrate: List[Dict[str, Any]] | None = None,
                 block_lookup: Dict[str, "FluidBlock"] | None = None,
                 tol: float = 1e-6):

        if equilibrate is None:
            # default: split *equally* between the outlets
            equilibrate = [{k: "diff"} for k in st_out]
        else:
            # sanity: every outlet must appear once and only once
            eq_keys = [list(d.keys())[0] for d in equilibrate]
            if sorted(eq_keys) != sorted(st_out):
                raise ValueError("every st_out must appear once in *equilibrate*")

        self.name = name
        self._st_in = list(st_in)
        self._st_out = list(st_out)
        self._medium = medium
        self._equilibrate = equilibrate
        self._blkmap = block_lookup or {}
        self._tol = tol

        super().__init__(medium)

        # ───── metadata expected by *EngineNetwork* ──────────────────
        self.station_inputs = list(st_in)
        self.station_outputs = list(st_out)
        self.signal_inputs: List[str] = []
        self.signal_outputs: List[str] = []

    # ------------------------------------------------------------------
    #  PRIVATE HELPERS
    # ------------------------------------------------------------------
    def _enthalpy(self, T: float, p: float) -> float:
        """Convenience wrapper around CoolProp."""
        return CP.PropsSI("Hmass", "T", T - 1.0e-3, "P", p, self._medium)

    # .................................................................
    def _simulate_path(self, blk_names: List[str], st_key: str, p_mix: float,
                        T_mix: float, mdot: float, stations: Dict[str, "psf.common.Station"],
                        signals: Dict[str, float]) -> float:
        """Return the *exit pressure* of one branch for a trial mdot.

        A deep‑copy of the dictionaries isolates the internal probe from the
        live network state.  Only the *p* field of **all** stations that the
        blocks touch is needed, so we keep the implementation lightweight.
        """

        # check that all block names are known
        try:
            path_blocks = [self._blkmap[n] for n in blk_names]
        except KeyError as missing:
            raise KeyError(f"JunctionBlock '{self.name}': unknown block '{missing.args[0]}'")

        loc_st = deepcopy(stations)
        loc_sg = deepcopy(signals)

        loc_st[st_key] = psf.common.Station(p=p_mix, T=T_mix, mdot=mdot)

        for blk in path_blocks:
            ds, sg = blk.compute(loc_st, loc_sg)
            loc_st.update(ds)
            loc_sg.update(sg)

        # the final block’s first advertised outlet is taken as the branch exit
        end_key = path_blocks[-1].station_outputs[0]
        return loc_st[end_key].p

    # .................................................................
    def _bisect(self, func, lo: float, hi: float, *, max_iter: int = 60) -> float:
        f_lo = func(lo)
        f_hi = func(hi)
        if f_lo == 0.0:
            return lo
        if f_hi == 0.0:
            return hi
        if f_lo * f_hi > 0.0:
            raise RuntimeError("Bisection failing to bracket the root.")

        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            f_mid = func(mid)
            if abs(f_mid) < self._tol:
                return mid
            if f_lo * f_mid < 0.0:
                hi, f_hi = mid, f_mid
            else:
                lo, f_lo = mid, f_mid
        return 0.5 * (lo + hi)  # best effort

    # ------------------------------------------------------------------
    #  MAIN ENTRY POINT
    # ------------------------------------------------------------------
    def compute(self, stations: Dict[str, "psf.common.Station"],
                signals: Dict[str, float]):

        # ─── 1. *merge* all inlet streams ────────────────────────────
        mdot_tot = 0.0
        h_sum = 0.0
        p_mix = min(stations[k].p for k in self._st_in)

        for k in self._st_in:
            st = stations[k]
            mdot_tot += st.mdot
            h_sum += st.mdot * self._enthalpy(st.T, st.p)

        if mdot_tot <= 0.0:
            raise ValueError(f"{self.name}: total mdot must be positive")

        h_mix = h_sum / mdot_tot
        T_mix = CP.PropsSI("T", "Hmass", h_mix, "P", p_mix, self._medium)

        # ─── 2. decide *split* of mdot_tot across the outlets ────────
        # classify branches
        cfg_by_out = {list(d.keys())[0]: list(d.values())[0] for d in self._equilibrate}

        branch_mdots: Dict[str, float] = {}

        # CASE A – exactly two *hydraulic* branches that re‑merge further
        lists = [v for v in cfg_by_out.values() if isinstance(v, list)]
        if len(lists) == 2 and "diff" not in cfg_by_out.values():
            # equalise exit pressure between the two lists
            out1, out2 = [k for k, v in cfg_by_out.items() if isinstance(v, list)]
            blks1, blks2 = cfg_by_out[out1], cfg_by_out[out2]

            def g(m1):
                p1 = self._simulate_path(blks1, out1, p_mix, T_mix, m1, stations, signals)
                p2 = self._simulate_path(blks2, out2, p_mix, T_mix, mdot_tot - m1, stations, signals)
                return p1 - p2

            # bracket by forcing m1→0 and m1→mdot_tot
            m1 = self._bisect(g, 1.0e-12, mdot_tot - 1.0e-12)
            branch_mdots[out1] = m1
            branch_mdots[out2] = mdot_tot - m1

        # CASE B – one hydraulic list + remainder branch(es) marked "diff"
        else:
            mdot_left = mdot_tot
            diff_outlets = []
            for out_key, method in cfg_by_out.items():
                if isinstance(method, list):
                    blks = method
                    p_target = stations[out_key].p  # current pressure of downstream node

                    def f(m):
                        p_exit = self._simulate_path(blks, out_key, p_mix, T_mix, m, stations, signals)
                        return p_exit - p_target

                    try:
                        m_solution = self._bisect(f, 1.0e-12, mdot_tot)
                    except RuntimeError:
                        m_solution = 0.5 * mdot_tot  # fall back to half‑split
                    branch_mdots[out_key] = m_solution
                    mdot_left -= m_solution
                else:  # "diff"
                    diff_outlets.append(out_key)

            # spread remaining flow over the diff branches (usually just one)
            if diff_outlets:
                share = mdot_left / len(diff_outlets)
                for o in diff_outlets:
                    branch_mdots[o] = share
            if abs(sum(branch_mdots.values()) - mdot_tot) / mdot_tot > 1e-6:
                raise RuntimeError("mass‑flow bookkeeping error inside JunctionBlock")

        # ─── 3. build output Station objects ─────────────────────────
        st_out = {}
        for out_key, m in branch_mdots.items():
            st_out[out_key] = psf.common.Station(p=p_mix, T=T_mix, mdot=m)

        # this block does not emit scalar signals
        return st_out, {}

    # ------------------------------------------------------------------
    def post_process(self, stations, signals):
        return {}  # nothing special here
 

class MassFlowSplitterBlock(FluidBlock):
    """
    Split one incoming stream into N outgoing branches.

    Parameters
    ----------
    name        : str
                 Unique tag.  Used only for bookkeeping.
    st_in       : str
                 Key of the inlet Station.
    st_outs     : list[str]
                 Keys of the outlet Stations (length = N).
    fractions   : list[float] | None
                 Fixed mass-flow fractions that sum to 1.0, OR
                 None if you want to supply them at run-time
                 through scalar signals (see *dynamic split* below).
    frac_keys   : list[str] | None
                 If you choose a dynamic split, give one signal key
                 per outlet.  The network will read those each pass.
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
            stations_out[st_key] = psf.common.Station(
                p    = st_i.p,         # no Δp inside the node itself
                T    = st_i.T,
                mdot = f * mdot
            )

        return stations_out, {}

class MassFlowMergerBlock(FluidBlock):
    """
    Combine multiple inlet streams into one outlet.

    NOTE:  Assumes identical fluid species.
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
            self.st_out: psf.common.Station(p = p_ref, T = T_mix, mdot = mdot_tot)
        }, {}


# Signal blocks
class TransmissionBlock(SignalBlock):
    """Sums shaft-power keys → one P_required key for the turbine."""
    def __init__(self, name, sink_keys, out_key="P_required"): # TODO: don't hardcode P_required
        self.station_inputs   = []
        self.station_outputs  = []
        self.signal_inputs    = list(sink_keys)
        self.signal_outputs   = [out_key]
        self.sink_keys = sink_keys
        self.out_key   = out_key
        self.name=name

    def compute(self, st, sg):
        P = sum(sg[k] for k in self.sink_keys)
        return {}, {self.out_key: P}

