"""Coupled film + regenerative cooling solver.

The hot-side boundary condition can come from either the hot combustion gas,
a film cooling solution, or a mix of the two along the chamber.

Three usage modes fall out of the same machinery:

* **Regenerative only** -- call :func:`coupled_steady_heating_analysis` with
  ``film=False`` (or on a chamber without ``film_cooling``).
* **Film only** -- call :func:`solve_film_cooling`, which runs the Grisson
  film model on its own and returns the liquid/gaseous film results.
* **Both** -- call :func:`coupled_steady_heating_analysis` on a chamber that
  carries a :class:`~pyskyfire.regen.thrust_chamber.FilmCooling` object. The
  film model is solved first, and the resulting wall heat flux replaces the
  hot-gas heat flux everywhere the film covers the wall.

How the coupling works
----------------------

Every hot-side model in here reduces to a wall heat flux ``qpp_hot``, so the
wall/coolant stack downstream never needs to know which model produced it.
:class:`FilmHeatFluxModel` decides, per station, which one does:

==================  =====================================================
Region              Wall heat flux
==================  =====================================================
upstream of the     enthalpy-driven Bartz, unchanged:
injector            :math:`q'' = h_g (H_\\text{aw} - H_\\text{hw})`
liquid film         conduction into the saturated film:
                    :math:`q'' = h_\\text{liq} (T_\\text{sat} - T_\\text{hw})`
gaseous film        Bartz coefficient, film recovery temperature:
                    :math:`q'' = h_\\text{hot} (T_\\text{aw,film} - T_\\text{hw}) + q''_\\text{rad}`
==================  =====================================================

In the gaseous-film region only the *driving* state is replaced, by the film's
adiabatic wall temperature from the Grisson march. The gas-side coefficient
itself is left at its bare-wall Bartz value, which is conservative: a real film
both thickens the boundary layer and blows through it, and both effects reduce
it.

Note in particular that Grisson's ``h_conv`` for the gaseous film is *not*
used as the wall coefficient. It is :math:`C_p \\, dM_\\text{bl,e}/dx`, the
rate at which free-stream gas is entrained into the boundary layer, and it
runs several times larger than the wall-side Bartz value; feeding it to the
wall would make the film appear to *heat* the chamber.

The coupling is **one-way**: the Grisson model assumes every joule reaching the
liquid film goes into evaporating it, so the film solution does not respond to
heat lost through the wall into the coolant. Dryout location is therefore a
slight under-estimate of the film's reach. Radiation transmitted *through* the
liquid film to the wall (the complement of ``liquid_absorptivity``) is not
modelled either, so the liquid-film region is optimistic by of order
0.1 MW m⁻².

The one empirical knob: ``h_liquid_wall``
-----------------------------------------

Behind a liquid film the wall is bathed in saturated liquid, so the driving
temperature collapses from the gas recovery temperature (thousands of kelvin)
to :math:`T_\\text{sat}` (hundreds). What is *not* pinned down by the Grisson
model is the film-to-wall conductance, because the model never resolves the
film thickness. ``h_liquid_wall`` supplies it, defaulting to
:data:`DEFAULT_H_LIQUID_WALL` -- an order-of-magnitude estimate for a thin
saturated film (:math:`k_l/\\delta` with :math:`k_l \\sim 0.15`
W m⁻¹ K⁻¹ and :math:`\\delta \\sim 30` µm). Set it very large (``1e6``) to
recover the textbook idealisation :math:`T_\\text{hw} = T_\\text{sat}`, or pass
a callable ``h(x)`` for a spatially varying estimate. It is the number to
calibrate first if you have wall-temperature data.

Note also that the film mass flow is *not* subtracted from the regenerative
circuit's mass flow here -- set ``BoundaryConditions.mdot_coolant`` to whatever
actually reaches the channels.

Coolant properties are evaluated at the coolant pressure
--------------------------------------------------------

:meth:`CoupledHeatExchangerPhysics.cold_side_coefficients` evaluates the
coolant film properties at the coolant pressure. The coolant pressure is what
the march already tracks and what ``coolant_friction_rate`` and
``coolant_enthalpy_rate`` already use, so it is threaded through here too.

The coolant march: enthalpy where possible, temperature otherwise
-----------------------------------------------------------------

The coolant state is advanced on **specific enthalpy**, so the march stays
valid through the boiling dome and reports vapour quality. Temperature is not a
usable state variable there -- inside the dome it is pinned to
:math:`T_\\text{sat}(p)` while the fluid keeps absorbing heat -- and a
temperature march stepping past saturation silently lands on the superheated
vapour branch, with a ~19x density error. See
:mod:`pyskyfire.regen.coolant_state` for the full argument and for what the
two-phase treatment does and does not cover.

Two coolants keep the old temperature march, via
:func:`~pyskyfire.regen.coolant_state.probe_coolant`:

* **mixtures**, because a mixture ``(h, p)`` flash costs ~374 ms against ~57 us
  for a pure fluid and would put a 100-node run into the hours;
* **backends without saturation data**, such as ``INCOMP::``.

On those the solver behaves exactly as it did before -- no quality, no dome --
and says so through a ``RuntimeWarning`` plus a printed note when ``output`` is
enabled.

What is still missing is **two-phase heat transfer**. Inside the dome
``h_cold`` remains the single-phase Colburn value on saturated-liquid
properties, so it captures neither the flow-boiling enhancement nor dryout, and
nucleate boiling and CHF are not modelled at all. Read the quality profile as
telling you *where* the coolant boils, not how hot the wall gets there. Nodes
where the energy balance fails to close are flagged by a ``RuntimeWarning``, so
check for those before trusting a profile.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from functools import wraps
from numbers import Integral
from typing import Callable, Optional, Sequence, Union

import numpy as np
from scipy.optimize import least_squares

from . import physics
from .coolant_state import (
    _SATURATION_NUDGE,
    PHASE_SUPERCRITICAL,
    PHASE_TWO_PHASE,
    CoolantState,
    _off_saturation,
)
from .film_solver import (
    GaseousFilmResults,
    GrissonFilmCoolingModel,
    LiquidFilmResults,
)

ResidualResult = tuple[np.ndarray, np.ndarray] | tuple[None, None]


def _clear_gas_cache_after(func):
    """Bound combustion-property cache lifetime to one simulation call."""
    @wraps(func)
    def wrapped(thrust_chamber, *args, **kwargs):
        try:
            return func(thrust_chamber, *args, **kwargs)
        finally:
            combustion_transport = getattr(
                thrust_chamber,
                "combustion_transport",
                None,
            )
            clear_cache = getattr(combustion_transport, "clear_cache", None)
            if callable(clear_cache):
                clear_cache()

    return wrapped


@dataclass
class BoundaryConditions:
    """Boundary conditions for the coolant-side inlet.

    Parameters
    ----------
    T_coolant_in : float
        Coolant static temperature at the channel inlet [K].
    p_coolant_in : float
        Coolant static pressure at the channel inlet [Pa].
    mdot_coolant : float
        Coolant mass-flow rate through the cooling channel [kg s⁻¹].
    """
    T_coolant_in: float
    p_coolant_in: float
    mdot_coolant: float


@dataclass
class RegenResult:
    """Regenerative-cooling result with optional film-cooling data.

    The coolant-state fields (``coolant_enthalpy`` onwards) describe the bulk
    coolant along the circuit. ``coolant_enthalpy`` and ``coolant_quality`` are
    all-``nan`` when the coolant fell back to the temperature march, and
    ``coolant_quality`` is ``nan`` at any station outside the boiling dome --
    including supercritical stations, where no dome exists. Use
    ``coolant_phase`` to tell those cases apart.
    """

    circuit_name: str
    circuit_index: int
    x: np.ndarray
    T: np.ndarray
    T_static: np.ndarray
    T_stagnation: np.ndarray
    p_static: np.ndarray
    p_stagnation: np.ndarray
    dQ_dA: np.ndarray
    velocity: np.ndarray
    h_hot: np.ndarray
    h_hot_enthalpy: np.ndarray
    h_cold: np.ndarray
    T_aw_hot: np.ndarray
    residuals: ResidualResult
    qpp_hot: Optional[np.ndarray] = None
    T_drive: Optional[np.ndarray] = None
    film_regime: Optional[np.ndarray] = None
    liquid_film: Optional[LiquidFilmResults] = None
    gaseous_film: Optional[GaseousFilmResults] = None
    #: Bulk coolant specific enthalpy [J kg⁻¹]; ``nan`` under the T-march.
    coolant_enthalpy: Optional[np.ndarray] = None
    #: Bulk vapour mass fraction; ``nan`` outside the dome.
    coolant_quality: Optional[np.ndarray] = None
    #: Local coolant saturation temperature [K]; ``nan`` where undefined.
    coolant_T_sat: Optional[np.ndarray] = None
    #: Per-station ``PHASE_*`` label from :mod:`pyskyfire.regen.coolant_state`.
    coolant_phase: Optional[np.ndarray] = None
    #: True when the coolant advanced on enthalpy rather than temperature.
    enthalpy_march: bool = False
    #: Axial stations used for the nonlinear wall-temperature balances [m].
    x_wall: Optional[np.ndarray] = None
    #: Axial stations where heat flux and hot-side diagnostics are reported [m].
    x_heat_flux: Optional[np.ndarray] = None
    #: Axial stations used for the lumped coolant-state march [m].
    x_coolant: Optional[np.ndarray] = None
    #: Final scaled energy-balance residual at every wall node [-].
    wall_residual_scaled: Optional[np.ndarray] = None
    #: True where ``wall_residual_scaled <= RESIDUAL_TOL``.
    wall_converged: Optional[np.ndarray] = None

__all__ = [
    "BoundaryConditions",
    "RegenResult",
    "DEFAULT_H_LIQUID_WALL",
    "RESIDUAL_TOL",
    "SATURATION_NUDGE",
    "FilmHeatFluxModel",
    "CoupledHeatExchangerPhysics",
    "solve_film_cooling",
    "solve_coupled_heat_exchanger",
    "coupled_steady_heating_analysis",
]


#: Default film-to-wall conductance behind a liquid film [W m⁻² K⁻¹].
DEFAULT_H_LIQUID_WALL = 5.0e3

#: Temperature step used to leave CoolProp's saturation-line exclusion band [K].
#: Kept as an alias so the two modules cannot drift.
SATURATION_NUDGE = _SATURATION_NUDGE

#: Scaled energy-balance residual above which a node is reported as suspect.
#: The residual is normalised by the node's reference heat load, so this is a
#: relative energy imbalance of 0.01%.
RESIDUAL_TOL = 1.0e-4

#: Maximum fixed-point sweeps used to settle one coolant pressure segment.
#: Friction, acceleration and the equation of state are mutually coupled across
#: a segment: the downstream density sets the velocity change, which sets the
#: pressure, which sets the density.
PRESSURE_SWEEPS = 8

#: Pressure change below which a segment's fixed-point sweep is converged [Pa].
PRESSURE_SWEEP_TOL = 1.0

# Regime labels reported per station.
REGIME_GAS = "gas"
REGIME_LIQUID_FILM = "liquid_film"
REGIME_GASEOUS_FILM = "gaseous_film"


class FilmHeatFluxModel:
    """Turn a Grisson film solution into a hot-side boundary condition.

    Parameters
    ----------
    liquid_results : LiquidFilmResults
        Output of the liquid-film march.
    gaseous_results : GaseousFilmResults or None
        Output of the gaseous-film march, empty/``None`` when the liquid film
        never dries out inside the domain.
    x_injection : float
        Axial station where the film is injected [m]. No film upstream of it.
    h_liquid_wall : float or callable, optional
        Film-to-wall conductance behind the liquid film [W m⁻² K⁻¹], either a
        constant or ``h(x)``. See the module docstring.

    Notes
    -----
    ``x`` outside the solved film range is clamped to the nearest solved
    station, so a regenerative grid finer than the film grid is fine.
    """

    def __init__(
        self,
        liquid_results: LiquidFilmResults,
        gaseous_results: Optional[GaseousFilmResults],
        x_injection: float,
        h_liquid_wall: Union[float, Callable[[float], float]] = DEFAULT_H_LIQUID_WALL,
    ):
        self.liquid_results = liquid_results
        self.gaseous_results = gaseous_results
        self.x_injection = float(x_injection)
        self.h_liquid_wall = h_liquid_wall

        # --- liquid branch ------------------------------------------------
        self._x_liq = np.asarray(liquid_results.x, dtype=float)
        self._T_sat = np.asarray(liquid_results.T_film, dtype=float)
        self._has_liquid = self._x_liq.size > 0

        if self._has_liquid:
            # Everything up to dryout is wetted; without dryout the liquid film
            # simply survives to the end of the solved range.
            self._x_liquid_end = (
                float(liquid_results.x_dryout)
                if liquid_results.x_dryout is not None
                else float(self._x_liq[-1])
            )
        else:
            self._x_liquid_end = self.x_injection

        # --- gaseous branch -----------------------------------------------
        x_gas = np.asarray(
            gaseous_results.x if gaseous_results is not None else [],
            dtype=float,
        )
        if gaseous_results is not None and x_gas.size > 0:
            self._has_gaseous = True
            self._x_gas = x_gas
            # The film-cooled adiabatic wall temperature is the whole point of
            # the gaseous march; its h_conv is an entrainment conductance and
            # is deliberately not used here (see the module docstring).
            self._T_aw_gas = np.asarray(gaseous_results.T_aw, dtype=float)
            self._q_rad_gas = np.asarray(gaseous_results.Q_rad, dtype=float)
            self._x_gas_end = float(self._x_gas[-1])
        else:
            self._has_gaseous = False
            self._x_gas_end = self._x_liquid_end

    # ------------------------------------------------------------------
    # construction helpers
    # ------------------------------------------------------------------
    @classmethod
    def from_thrust_chamber(
        cls,
        thrust_chamber,
        boundary_conditions,
        x_array: Optional[Sequence[float]] = None,
        h_liquid_wall: Union[float, Callable[[float], float]] = DEFAULT_H_LIQUID_WALL,
    ) -> "FilmHeatFluxModel":
        """Solve the Grisson film model for ``thrust_chamber`` and wrap it.

        Parameters
        ----------
        thrust_chamber : Any
            Chamber carrying a ``film_cooling`` object.
        boundary_conditions : BoundaryConditions
            Film injection temperature/pressure are taken from here.
        x_array : sequence of float, optional
            Grid for the film march. Defaults to the contour's own ``xs``.
        h_liquid_wall : float or callable, optional
            See the class docstring.
        """
        liquid, gaseous = solve_film_cooling(
            thrust_chamber, boundary_conditions, x_array=x_array
        )
        return cls(
            liquid_results=liquid,
            gaseous_results=gaseous,
            x_injection=float(thrust_chamber.film_cooling.x),
            h_liquid_wall=h_liquid_wall,
        )

    # ------------------------------------------------------------------
    # queries
    # ------------------------------------------------------------------
    def regime(self, x: float) -> str:
        """Return the hot-side regime at ``x``.

        Returns
        -------
        str
            ``'gas'`` where the wall sees the combustion gas directly,
            ``'liquid_film'`` behind the liquid film, ``'gaseous_film'``
            behind its vapour downstream of dryout.
        """
        x = float(x)
        if x < self.x_injection:
            return REGIME_GAS
        if self._has_liquid and x <= self._x_liquid_end:
            return REGIME_LIQUID_FILM
        if self._has_gaseous and x <= self._x_gas_end:
            return REGIME_GASEOUS_FILM
        # Past the end of the solved film domain the wall is bare again.
        return REGIME_GAS

    def covers(self, x: float) -> bool:
        """True when the film -- liquid or gaseous -- sets the wall flux at ``x``."""
        return self.regime(x) != REGIME_GAS

    def _h_liquid(self, x: float) -> float:
        h = self.h_liquid_wall
        return float(h(x)) if callable(h) else float(h)

    def state(self, x: float) -> Optional[dict]:
        """Film boundary condition at ``x``, or ``None`` where the wall is bare.

        Parameters
        ----------
        x : float
            Axial coordinate [m].

        Returns
        -------
        dict or None
            ``regime``, ``T_drive`` [K], ``q_rad`` [W m⁻²] and ``h_wall``
            [W m⁻² K⁻¹]. ``h_wall`` is ``None`` in the gaseous-film region,
            where the caller supplies the gas-side coefficient instead.
            ``None`` when the station is not film-cooled at all.
        """
        regime = self.regime(x)
        if regime == REGIME_GAS:
            return None

        if regime == REGIME_LIQUID_FILM:
            return {
                "regime": regime,
                "T_drive": float(np.interp(x, self._x_liq, self._T_sat)),
                # Radiation reaching the wall through the film is not modelled.
                "q_rad": 0.0,
                "h_wall": self._h_liquid(x),
            }

        return {
            "regime": regime,
            "T_drive": float(np.interp(x, self._x_gas, self._T_aw_gas)),
            "q_rad": float(np.interp(x, self._x_gas, self._q_rad_gas)),
            "h_wall": None,
        }


class CoupledHeatExchangerPhysics:
    """Heat-exchanger physics whose hot side may be gas- or film-driven.

    The bare-wall hot side uses the enthalpy-driven Bartz model. Where a film
    covers the wall, :meth:`hot_side_coefficients` substitutes the appropriate
    liquid- or gaseous-film boundary condition. The shared wall, coolant, and
    pressure models consume the resulting ``qpp_hot`` regardless of its source.

    Parameters
    ----------
    thrust_chamber, boundary_conditions, circuit_index
        As for the base class.
    film_model : FilmHeatFluxModel or None, optional
        Film boundary condition. ``None`` reproduces the base class exactly.
    heat_curvature_correction : bool, optional
        Apply the RL10 curvature multiplier to coolant-side heat transfer.
    pressure_curvature_correction : bool, optional
        Apply the RTE/Ito curvature multiplier to Darcy pressure loss.
    """

    def __init__(
        self,
        thrust_chamber,
        boundary_conditions,
        circuit_index,
        film_model=None,
        coolant_state: Optional[CoolantState] = None,
        heat_curvature_correction: bool = True,
        pressure_curvature_correction: bool = True,
    ):
        self.thrust_chamber = thrust_chamber
        self.boundary_conditions = boundary_conditions
        self.circuit_index = circuit_index
        self.film_model = film_model
        self.heat_curvature_correction = bool(heat_curvature_correction)
        self.pressure_curvature_correction = bool(
            pressure_curvature_correction
        )

        circuit = thrust_chamber.cooling_circuits[circuit_index]
        self.coolant_state = coolant_state or CoolantState(
            circuit.coolant_transport, boundary_conditions.p_coolant_in
        )

    def _gas_side_coefficients(self, x, T_hw, T_gr_mode="reference"):
        """Return the bare-gas, enthalpy-driven Bartz heat-transfer state."""
        ct = self.thrust_chamber.combustion_transport

        D_hyd = 2 * self.thrust_chamber.contour.r(x)
        A_chmb = self.thrust_chamber.contour.A(x)

        gas = ct.get_state(x)
        T_g = gas.T
        H_g = gas.h
        M_g = gas.M
        a_g = gas.a

        try:
            H_hw = ct.get_h(x, T=T_hw)
        except Exception:  # noqa: BLE001
            H_hw = H_g

        dynamic_enthalpy = 0.5 * M_g**2 * a_g**2
        H_gr = 0.5 * (H_hw + H_g) + 0.18 * dynamic_enthalpy

        ref = ct.get_state(x, h=H_gr)
        if T_gr_mode == "reference":
            T_gr = ref.T
        elif T_gr_mode == "film":
            T_gr = 0.5 * (T_hw + T_g)
        else:
            raise ValueError("T_gr_mode must be 'reference' or 'film'")

        h_gr = physics.h_gas_bartz_enthalpy_driven(
            ref.k,
            D_hyd,
            ref.cp,
            ref.mu,
            ct.mdot,
            A_chmb,
            T_g,
            T_gr,
        ) * self.thrust_chamber.h_gas_corr
        h_g = h_gr / ref.cp

        recovery_factor = ref.Pr ** (1.0 / 3.0)
        H_aw = H_g + recovery_factor * dynamic_enthalpy

        # Recover T_aw by inverting the equilibrium equation of state at H_aw
        # rather than from the frozen-gamma isentropic relation. The two agree
        # in the chamber but diverge hard in the nozzle: with equilibrium
        # chemistry, gamma is a local property and the isentropic formula is no
        # longer a statement about conserved total enthalpy, so it reports a
        # recovery temperature above the chamber stagnation value. That inflates
        # the T_aw - T_hw span and therefore deflates the reported h_hot.
        try:
            T_aw = ct.get_state(x, h=H_aw).T
        except Exception:  # noqa: BLE001
            T_aw = physics.T_aw(
                gamma=gas.gamma,
                M_inf=M_g,
                T_inf=T_g,
                Pr=ref.Pr,
            )
        qpp_hot = h_g * (H_aw - H_hw)

        delta_T = T_aw - T_hw
        h_hot = qpp_hot / delta_T if abs(delta_T) > 1.0e-9 else float("nan")

        return {
            "h_hot": h_hot,
            "h_g": h_g,
            "h_gr": h_gr,
            "T_aw": T_aw,
            "qpp_hot": qpp_hot,
            "H_gr": H_gr,
            "H_aw": H_aw,
            "H_hw": H_hw,
            "T_gr": T_gr,
        }

    # ------------------------------------------------------------------
    # coolant side: evaluated at the coolant pressure, not the gas pressure
    # ------------------------------------------------------------------
    def _coolant_pressure(self, p_cool):
        """Resolve the pressure to evaluate coolant properties at [Pa]."""
        if p_cool is not None:
            return float(p_cool)
        return float(self.boundary_conditions.p_coolant_in)

    _off_saturation = staticmethod(_off_saturation)

    def bulk_density(self, T_cool, p_cool, quality=None):
        """Bulk coolant density [kg m⁻³], homogeneous inside the boiling dome.

        Parameters
        ----------
        T_cool : float
            Bulk coolant temperature [K].
        p_cool : float
            Coolant pressure [Pa].
        quality : float, optional
            Vapour mass fraction. ``None`` or ``nan`` means single-phase, in
            which case density comes from an ordinary ``(T, p)`` lookup.

        Notes
        -----
        The two-phase value is the homogeneous, no-slip mixture density. It
        captures the large acceleration effect through the dome but not slip,
        so it is an under-estimate of the true momentum flux.
        """
        p = self._coolant_pressure(p_cool)
        if quality is not None and np.isfinite(quality):
            sat = self.coolant_state.saturation(p)
            if sat is not None:
                return sat.homogeneous_density(quality)

        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]
        return _off_saturation(circuit.coolant_transport.get_rho, T_cool, p)

    def cold_side_coefficients(self, x, T_cw, T_cool, p_cool=None, quality=None):
        """Coolant-side heat transfer coefficient at ``x``.

        Coolant properties are evaluated at the local coolant pressure.

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        T_cw : float
            Coolant-side wall temperature [K].
        T_cool : float
            Bulk coolant temperature [K].
        p_cool : float, optional
            Local coolant pressure [Pa]. Defaults to the inlet pressure, which
            is the best available guess outside a march.
        quality : float, optional
            Bulk vapour mass fraction, used for the two-phase bulk density.
            ``None``/``nan`` outside the dome.

        Returns
        -------
        dict
            ``h_cold``, ``phi_curv``, ``Re_c``.

        Notes
        -----
        This is a **single-phase** correlation everywhere, including inside the
        boiling dome, where it is evaluated on saturated-liquid properties. It
        therefore models neither flow-boiling enhancement nor dryout. See
        :mod:`pyskyfire.regen.coolant_state`.
        """
        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]

        p = self._coolant_pressure(p_cool)
        mdot_c_single_channel = self.mdot_per_channel()
        T_coolant_film = 0.5 * (T_cool + T_cw)

        ct = circuit.coolant_transport
        # Clamped to saturated liquid once the film temperature crosses T_sat,
        # rather than falling through to CoolProp's vapour branch.
        k_cf, Cp_cr, mu_cf = self.coolant_state.film_properties(
            T_coolant_film, p, (ct.get_k, ct.get_cp, ct.get_mu)
        )

        D_c = circuit.Dh_coolant(x)
        A_channel = circuit.A_coolant(x)

        rho_bulk = self.bulk_density(T_cool, p, quality)
        if quality is not None and np.isfinite(quality):
            sat = self.coolant_state.saturation(p)
            mu_bulk = sat.mu_f if sat is not None else _off_saturation(ct.get_mu, T_cool, p)
        else:
            mu_bulk = _off_saturation(ct.get_mu, T_cool, p)
        u_c = physics.u_coolant(rho_bulk, mdot_c_single_channel, A_channel)
        Re_c = physics.reynolds(rho_bulk, u_c, D_c, mu_bulk)

        if self.heat_curvature_correction:
            R_curv = circuit.radius_of_curvature(x)
            phi_curv = physics.phi_curv(Re_c, D_c, R_curv)
        else:
            phi_curv = 1.0

        h_c = physics.h_coolant_colburn(
            k_cf,
            D_c,
            Cp_cr,
            mu_cf,
            mdot_c_single_channel,
            A_channel,
            phi_curv=phi_curv,
        ) * self.thrust_chamber.h_cold_corr

        return {
            "h_cold": h_c,
            "phi_curv": phi_curv,
            "Re_c": Re_c,
        }

    def dQ_cold_dx(self, x, T_cw, T_cool, p_cool=None, quality=None):
        """Heat removed by the coolant per unit length [W m⁻¹].

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        T_cw : float
            Coolant-side wall temperature [K].
        T_cool : float
            Bulk coolant temperature [K].
        p_cool : float, optional
            Local coolant pressure [Pa], forwarded to
            :meth:`cold_side_coefficients`.
        quality : float, optional
            Bulk vapour mass fraction, forwarded to
            :meth:`cold_side_coefficients`.
        """
        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]
        h_c = self.cold_side_coefficients(
            x, T_cw, T_cool, p_cool=p_cool, quality=quality
        )["h_cold"]

        T_rep = 0.5 * (T_cw + T_cool)
        R_cool_per_len = circuit.R_coolant_per_len(x, h_c=h_c, T_wall_rep=T_rep)

        return (T_cw - T_cool) / R_cool_per_len

    def hot_side_coefficients(self, x, T_hw, T_gr_mode="reference"):
        """Hot-side coefficients, from the film where there is one.

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        T_hw : float
            Hot-side wall temperature [K].
        T_gr_mode : {'reference', 'film'}, optional
            Passed through to the gas-side model; ignored behind a liquid film.

        Returns
        -------
        dict
            Same keys as the base implementation, plus ``regime`` and
            ``q_rad``. Behind a liquid film the wall never sees the gas, so the
            Bartz diagnostics (``h_g``, ``h_gr``, ``H_gr``, ``H_aw``, ``H_hw``,
            ``T_gr``) are ``nan`` and no CEA solve is done.
        """
        film = self.film_model.state(x) if self.film_model is not None else None

        # --- behind a liquid film: the wall is bathed in saturated liquid ---
        if film is not None and film["regime"] == REGIME_LIQUID_FILM:
            T_drive = film["T_drive"]
            h_wall = film["h_wall"]
            nan = float("nan")
            return {
                "h_hot": h_wall,
                "h_g": nan,
                "h_gr": nan,
                "T_aw": T_drive,
                "qpp_hot": h_wall * (T_drive - float(T_hw)) + film["q_rad"],
                "H_gr": nan,
                "H_aw": nan,
                "H_hw": nan,
                "T_gr": nan,
                "q_rad": film["q_rad"],
                "regime": REGIME_LIQUID_FILM,
            }

        coeffs = self._gas_side_coefficients(x, T_hw, T_gr_mode=T_gr_mode)
        coeffs["regime"] = REGIME_GAS
        coeffs["q_rad"] = 0.0

        # --- behind the vapour: same coefficient, film recovery state -------
        if film is not None:
            T_drive = film["T_drive"]
            q_rad = film["q_rad"]

            # Temperature-driven rather than enthalpy-driven, deliberately.
            # The enthalpy form would need H at the film recovery temperature,
            # and equilibrium CEA enthalpy is badly behaved down there -- cp
            # can swing by an order of magnitude across a few tens of kelvin as
            # species freeze out, which kinks the residual and stalls the wall
            # solve. The base coefficient h_hot already carries a *mean* cp
            # over the whole wall-to-recovery span, which is smooth.
            qpp_hot = coeffs["h_hot"] * (T_drive - float(T_hw)) + q_rad

            coeffs.update(
                {
                    "T_aw": T_drive,
                    "qpp_hot": qpp_hot,
                    "H_aw": float("nan"),
                    "q_rad": q_rad,
                    "regime": REGIME_GASEOUS_FILM,
                }
            )

        return coeffs

    def dQ_hot_dx(self, x, T_hw):
        """Return hot-side heat flow per unit axial length [W m⁻¹]."""
        qpp_hot = self.hot_side_coefficients(x, T_hw)["qpp_hot"]
        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]
        return qpp_hot * circuit.dA_dx_thermal_exhaust(x)

    def dQ_cond_dx(self, x, T_hw, T_cw):
        """Return conduction through the wall stack per unit length [W m⁻¹]."""
        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]
        dA_dx_hot = circuit.dA_dx_thermal_exhaust(x)
        T_wall = 0.5 * (T_hw + T_cw)
        resistance = sum(
            wall.thickness(x) / (wall.material.get_k(T_wall) * dA_dx_hot)
            for wall in circuit.walls
        )
        return (T_hw - T_cw) / resistance

    def mdot_per_channel(self):
        """Coolant mass flow through a single channel [kg s⁻¹]."""
        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]
        return self.boundary_conditions.mdot_coolant / circuit.n_channels

    def coolant_enthalpy_rate(self, Q_cold):
        """Return the coolant specific-enthalpy change from ``Q_cold`` [J kg⁻¹].

        The state variable of the coolant march. Unlike
        :meth:`coolant_temperature_rate` this stays valid through the boiling
        dome, where added heat raises enthalpy at constant temperature.
        """
        return Q_cold / self.mdot_per_channel()

    def coolant_temperature_rate(self, T_cool, p_cool, Q_cold):
        """Return the coolant temperature change caused by ``Q_cold`` [K].

        Single-phase only -- it assumes all the heat is sensible. Used for
        coolants that cannot support the enthalpy march (mixtures, and backends
        without saturation data); see :func:`.coolant_state.probe_coolant`.
        """
        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]
        cp = circuit.coolant_transport.get_cp(T_cool, p_cool)
        return Q_cold / (self.mdot_per_channel() * cp)

    def bulk_velocity(self, x, T_cool, p_cool, quality=None):
        """Return the bulk coolant density and single-channel velocity at ``x``.

        Returns
        -------
        tuple of float
            ``(rho, u)`` in [kg m⁻³] and [m s⁻¹].
        """
        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]
        rho = self.bulk_density(T_cool, p_cool, quality)
        velocity = physics.u_coolant(
            rho, self.mdot_per_channel(), circuit.A_coolant(x)
        )
        return rho, velocity

    def coolant_friction_rate(self, x, T_cool, p_cool, quality=None):
        """Return the friction-driven coolant-pressure gradient [Pa m⁻¹].

        This is the irreversible term of the 1-D momentum equation only. The
        reversible acceleration term :math:`-\\rho u\\,du/dx` is applied by the
        march station to station, because it depends on the state at both ends
        of a segment rather than on a local gradient.

        Parameters
        ----------
        quality : float, optional
            Bulk vapour mass fraction. Inside the dome the homogeneous mixture
            density is used, which captures the acceleration through the dome
            but not the two-phase friction multiplier, so :math:`\\Delta p` is
            under-predicted at appreciable quality.
        """
        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]
        mdot_channel = self.mdot_per_channel()
        coolant = circuit.coolant_transport

        rho_cool = self.bulk_density(T_cool, p_cool, quality)
        two_phase = quality is not None and np.isfinite(quality)
        sat = self.coolant_state.saturation(p_cool) if two_phase else None
        mu_cool = sat.mu_f if sat is not None else coolant.get_mu(T_cool, p_cool)

        area = circuit.A_coolant(x)
        velocity = physics.u_coolant(rho_cool, mdot_channel, area)
        hydraulic_diameter = circuit.Dh_coolant(x)
        reynolds = physics.reynolds(
            rho_cool,
            velocity,
            hydraulic_diameter,
            mu_cool,
        )
        friction = physics.f_darcy(
            reynolds,
            hydraulic_diameter,
            x,
            circuit.roughness,
        )
        if self.pressure_curvature_correction:
            friction *= physics.phi_curv_friction(
                reynolds,
                hydraulic_diameter,
                circuit.radius_of_curvature(x),
            )

        return (
            -friction
            / hydraulic_diameter
            * rho_cool
            * velocity**2
            / 2
            * circuit.ds_dx(x)
        )

    def interface_temperatures(self, x, T_hw, T_cw):
        """Return wall-interface temperatures from the hot to the cold face."""
        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]
        dA_dx_hot = circuit.dA_dx_thermal_exhaust(x)
        qdx = self.dQ_cond_dx(x, T_hw, T_cw)
        T_wall = 0.5 * (T_hw + T_cw)
        temperatures = [T_hw]

        for wall in circuit.walls:
            resistance = wall.thickness(x) / (
                wall.material.get_k(T_wall) * dA_dx_hot
            )
            temperatures.append(temperatures[-1] - qdx * resistance)

        return temperatures


@_clear_gas_cache_after
def solve_film_cooling(
    thrust_chamber,
    boundary_conditions,
    x_array: Optional[Sequence[float]] = None,
) -> tuple[LiquidFilmResults, GaseousFilmResults]:
    """Run the film cooling model on its own.

    Parameters
    ----------
    thrust_chamber : Any
        Chamber carrying a ``film_cooling`` object.
    boundary_conditions : BoundaryConditions
        Film injection temperature/pressure come from here.
    x_array : sequence of float, optional
        Axial grid for the film march. Defaults to ``thrust_chamber.contour.xs``.

    Returns
    -------
    tuple
        ``(liquid_results, gaseous_results)``. The gaseous results are empty
        when the liquid film survives to the end of the grid.

    Raises
    ------
    ValueError
        If the chamber has no ``film_cooling`` attached.
    """
    if getattr(thrust_chamber, "film_cooling", None) is None:
        raise ValueError(
            "thrust_chamber.film_cooling is None; attach a FilmCooling object "
            "to run a film cooling analysis."
        )

    if x_array is None:
        x_array = np.asarray(thrust_chamber.contour.xs, dtype=float)

    model = GrissonFilmCoolingModel(thrust_chamber, boundary_conditions)
    return model.solve(x_array)


_NODE_GRID_KEYS = ("wall", "heat_flux", "coolant")


def _prepare_node_grids(nodes, circuit) -> list[np.ndarray]:
    """Return wall, heat-flux, and coolant grids in coolant-march order.

    An integer expands to three identical uniform grids. Explicit grids may be
    supplied as a ragged three-item sequence or as a mapping whose keys are
    ``wall``, ``heat_flux``, and ``coolant``.
    """
    if isinstance(nodes, Integral) and not isinstance(nodes, (bool, np.bool_)):
        count = int(nodes)
        if count < 2:
            raise ValueError("nodes must be at least 2")
        start, end = (
            (circuit.x_domain[0], circuit.x_domain[-1])
            if circuit.direction == 1
            else (circuit.x_domain[-1], circuit.x_domain[0])
        )
        uniform = np.linspace(float(start), float(end), count)
        return [uniform.copy() for _ in _NODE_GRID_KEYS]

    if isinstance(nodes, Mapping):
        missing = [key for key in _NODE_GRID_KEYS if key not in nodes]
        extra = [key for key in nodes if key not in _NODE_GRID_KEYS]
        if missing or extra:
            raise ValueError(
                "nodes mapping must contain exactly 'wall', 'heat_flux', and "
                f"'coolant'; missing={missing}, extra={extra}"
            )
        raw_grids = [nodes[key] for key in _NODE_GRID_KEYS]
    else:
        try:
            raw_grids = list(nodes)
        except TypeError as exc:
            raise TypeError(
                "nodes must be an integer, a three-item ragged sequence, or "
                "a mapping with wall/heat_flux/coolant keys"
            ) from exc
        if len(raw_grids) != 3:
            raise ValueError(
                "explicit nodes must contain three grids in the order "
                "[wall, heat_flux, coolant]"
            )

    grids = []
    for key, values in zip(_NODE_GRID_KEYS, raw_grids):
        grid = np.asarray(values, dtype=float)
        if grid.ndim != 1:
            raise ValueError(f"nodes[{key!r}] must be one-dimensional")
        if grid.size < 2:
            raise ValueError(f"nodes[{key!r}] must contain at least two stations")
        if not np.all(np.isfinite(grid)):
            raise ValueError(f"nodes[{key!r}] contains a non-finite station")
        if np.unique(grid).size != grid.size:
            raise ValueError(f"nodes[{key!r}] contains duplicate stations")
        grid = np.sort(grid)
        if circuit.direction == -1:
            grid = grid[::-1]
        grids.append(grid)

    wall, heat_flux, coolant = grids
    direction = float(circuit.direction)
    s_wall = direction * wall
    s_heat = direction * heat_flux
    s_cool = direction * coolant
    span = float(s_cool[-1] - s_cool[0])
    if span <= 0.0:
        raise ValueError("coolant nodes do not advance in the circuit direction")

    # Digitised reference curves commonly differ by a fraction of one plot
    # pixel at nominally shared end stations. Accommodate that without allowing
    # a thermal grid that genuinely extends beyond the coolant model.
    endpoint_tol = max(1.0e-9, 1.0e-3 * span)
    for key, s_grid in (("wall", s_wall), ("heat_flux", s_heat)):
        if s_grid[0] < s_cool[0] - endpoint_tol or s_grid[-1] > s_cool[-1] + endpoint_tol:
            raise ValueError(
                f"nodes[{key!r}] must lie within the coolant-node span "
                f"({coolant[0]:.6g} to {coolant[-1]:.6g} m)"
            )

    # Each lumped coolant segment must own at least one metal node. This is the
    # hybrid arrangement used by the RL10 system model and avoids silently
    # creating extra nonlinear wall solves at hidden stations.
    owners = np.searchsorted(s_cool, s_wall, side="right") - 1
    owners = np.clip(owners, 0, coolant.size - 2)
    missing_segments = [
        i for i in range(coolant.size - 1) if not np.any(owners == i)
    ]
    if missing_segments:
        raise ValueError(
            "every coolant segment must contain at least one wall node; "
            f"segments without wall nodes: {missing_segments}"
        )

    return grids


def _interp_axial(x_query, x, values):
    """One-dimensional interpolation independent of march direction."""
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float)
    order = np.argsort(x)
    return np.interp(x_query, x[order], values[order])


def _acceleration_pressure_drop(rho_up, velocity_up, rho_down, velocity_down):
    """Static-pressure drop across one segment from :math:`\\int\\rho u\\,du` [Pa].

    Trapezoidal in the mass flux :math:`G=\\rho u`, which makes the discrete
    form exact in both limits the term has to span in a heated, tapering
    channel: Bernoulli :math:`\\tfrac{1}{2}\\rho(u_2^2-u_1^2)` at constant
    density, and :math:`G^2(1/\\rho_2-1/\\rho_1)` at constant area.

    Note that at constant area this is twice the change in dynamic head, so
    reconstructing static pressure as :math:`p_0-\\tfrac{1}{2}\\rho u^2` from a
    friction-only stagnation march counts the heating-driven acceleration at
    half strength.
    """
    return (
        0.5
        * (rho_up * velocity_up + rho_down * velocity_down)
        * (velocity_down - velocity_up)
    )


@_clear_gas_cache_after
def solve_coupled_heat_exchanger(
    thrust_chamber,
    boundary_conditions,
    nodes,
    circuit_index,
    output,
    film_model: Optional[FilmHeatFluxModel] = None,
    log_residuals: bool = True,
    heat_curvature_correction: bool = True,
    pressure_curvature_correction: bool = True,
) -> RegenResult:
    """Solve wall heating on detailed nodes coupled to lumped coolant nodes.

    ``nodes`` is normalized immediately to three axial grids in the order
    ``wall``, ``heat_flux``, and ``coolant``. An integer creates three equal
    uniform grids. A ragged three-item sequence permits different counts, and
    a mapping with those three names provides the same capability without
    relying on positional order.

    Each coolant interval owns one or more wall nodes. The wall balances use
    the interval's lumped coolant state; their coolant-side heat rates are
    integrated over metal-node control volumes to advance the next coolant
    state. This mirrors the Binder et al. RL10 hybrid model, where multiple
    metal nodes were connected to each of five coolant segments.
    """
    circuit = thrust_chamber.cooling_circuits[circuit_index]
    x_wall, x_heat, x_cool = _prepare_node_grids(nodes, circuit)
    n_wall, n_heat, n_cool = len(x_wall), len(x_heat), len(x_cool)
    direction = float(circuit.direction)
    s_wall = direction * x_wall
    s_cool = direction * x_cool
    owners = np.searchsorted(s_cool, s_wall, side="right") - 1
    owners = np.clip(owners, 0, n_cool - 2)

    residual_log = [] if log_residuals else None
    iter_counter = np.zeros(n_wall, dtype=int)

    T_hw_arr = np.zeros(n_wall)
    T_cw_arr = np.zeros(n_wall)
    cold_qdx_arr = np.zeros(n_wall)
    wall_residual_scaled = np.full(n_wall, np.nan)

    T_cool_arr = np.zeros(n_cool)
    p_stagnation_arr = np.zeros(n_cool)
    p_static_arr = np.zeros(n_cool)
    velocity_arr = np.zeros(n_cool)
    coolant_enthalpy_arr = np.full(n_cool, np.nan)
    quality_arr = np.full(n_cool, np.nan)
    T_sat_arr = np.full(n_cool, np.nan)
    phase_arr = np.empty(n_cool, dtype=object)

    T_cool_in = boundary_conditions.T_coolant_in
    p_cool_in = boundary_conditions.p_coolant_in
    p_stagnation_arr[0] = p_cool_in

    coolant_state = CoolantState(circuit.coolant_transport, p_cool_in)
    capability = coolant_state.capability
    use_enthalpy = capability.enthalpy_march
    if not use_enthalpy:
        warnings.warn(
            f"coolant march fell back to temperature for circuit "
            f"{circuit.name!r}: {capability.reason} Coolant temperatures past "
            f"the saturation point are not physical.",
            RuntimeWarning,
            stacklevel=2,
        )

    if use_enthalpy:
        coolant_enthalpy_arr[0] = coolant_state.enthalpy(T_cool_in, p_cool_in)
        bulk = coolant_state.from_hp(coolant_enthalpy_arr[0], p_cool_in)
    else:
        bulk = coolant_state.from_tp(T_cool_in, p_cool_in)
    T_cool_arr[0] = bulk.T
    quality_arr[0] = bulk.quality
    T_sat_arr[0] = bulk.T_sat
    phase_arr[0] = bulk.phase

    physics_helper = CoupledHeatExchangerPhysics(
        thrust_chamber,
        boundary_conditions,
        circuit_index,
        film_model=film_model,
        coolant_state=coolant_state,
        heat_curvature_correction=heat_curvature_correction,
        pressure_curvature_correction=pressure_curvature_correction,
    )

    # ``p_coolant_in`` is the inlet *stagnation* pressure; the static value the
    # march carries is that minus the inlet dynamic head.
    rho_in, velocity_arr[0] = physics_helper.bulk_velocity(
        x_cool[0], T_cool_arr[0], p_cool_in, quality_arr[0]
    )
    p_static_arr[0] = p_cool_in - 0.5 * rho_in * velocity_arr[0] ** 2

    T_gas_first = thrust_chamber.combustion_transport.get_T(x_wall[0])
    T_hw_guess = 0.5 * (T_gas_first + T_cool_in)
    T_cw_guess = T_hw_guess
    solved_wall_count = 0

    if output is True:
        if heat_curvature_correction and getattr(circuit, "is_helical", False):
            warnings.warn(
                "The RL10 curvature correction applied to the coolant-side "
                "Colburn equation is being used on a helical channel. Its "
                "local-radius approximation may be inaccurate for large "
                "helix angles or curvature-heavy configurations.",
                RuntimeWarning,
                stacklevel=2,
            )
        covered = (
            sum(1 for x_i in x_heat if film_model.covers(x_i))
            if film_model else 0
        )
        print(
            "Started coupled heat exchanger simulation with "
            f"{n_wall} wall, {n_heat} heat-flux, and {n_cool} coolant nodes "
            f"({covered} film-cooled heat-flux nodes)"
        )
        if use_enthalpy:
            note = (
                "no boiling dome at inlet pressure (supercritical)"
                if bulk.phase == PHASE_SUPERCRITICAL
                else f"T_sat = {bulk.T_sat:.1f} K at inlet"
            )
            print(f"Coolant march: enthalpy, {note}")
        else:
            print(f"Coolant march: temperature -- {capability.reason}")

    def solve_wall_node(wall_index, coolant_index):
        nonlocal T_hw_guess, T_cw_guess, solved_wall_count
        x_i = x_wall[wall_index]
        T_cool_i = T_cool_arr[coolant_index]
        p_cool_i = p_stagnation_arr[coolant_index]
        quality_i = quality_arr[coolant_index]

        probe = physics_helper.hot_side_coefficients(x_i, T_hw_guess)
        T_drive_i = probe["T_aw"]
        if not np.isfinite(T_drive_i):
            T_drive_i = thrust_chamber.combustion_transport.get_T(x_i)
        q_rad_i = probe.get("q_rad", 0.0)
        h_probe = probe["h_hot"]
        if q_rad_i and np.isfinite(h_probe) and h_probe != 0.0:
            T_drive_i += q_rad_i / h_probe

        dA_dx_i = circuit.dA_dx_thermal_exhaust(x_i)
        if np.isfinite(probe["h_hot"]):
            q_ref = abs(probe["h_hot"]) * max(abs(T_drive_i - T_cool_i), 1.0)
        else:
            q_ref = abs(probe["qpp_hot"])
        Q_ref_per_len = max(q_ref * dA_dx_i, 1.0)

        def residuals_scaled(vars_):
            T_cw_trial, dT_trial = vars_
            T_hw_trial = T_cw_trial + dT_trial
            Q_hot_per_len = physics_helper.dQ_hot_dx(x_i, T_hw_trial)
            Q_cond_per_len = physics_helper.dQ_cond_dx(
                x_i, T_hw_trial, T_cw_trial
            )
            Q_cold_per_len = physics_helper.dQ_cold_dx(
                x_i,
                T_cw_trial,
                T_cool_i,
                p_cool=p_cool_i,
                quality=quality_i,
            )
            R1 = Q_hot_per_len - Q_cond_per_len
            R2 = Q_cond_per_len - Q_cold_per_len
            if residual_log is not None:
                iteration = iter_counter[wall_index]
                residual_log.append((wall_index, iteration, R1, R2))
                iter_counter[wall_index] += 1
            return np.array([R1 / Q_ref_per_len, R2 / Q_ref_per_len])

        T_lo = min(T_cool_i, T_drive_i) - 1.0
        T_hi = max(T_cool_i, T_drive_i) + 1.0
        span = max(1.0, T_hi - T_lo)
        dT_lo = 0.0 if T_drive_i >= T_cool_i else -span
        dT_hi = span
        x0 = np.array(
            [
                np.clip(T_cw_guess, T_lo, T_hi),
                np.clip(T_hw_guess - T_cw_guess, dT_lo, dT_hi),
            ],
            dtype=float,
        )
        solution = least_squares(
            residuals_scaled,
            x0=x0,
            bounds=([T_lo, dT_lo], [T_hi, dT_hi]),
            method="trf",
            loss="soft_l1",
            f_scale=1.0,
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
            max_nfev=200,
            diff_step=1e-4,
        )
        if not solution.success:
            raise RuntimeError(
                f"least_squares failed at wall node {wall_index} "
                f"(x={x_i:.6g}): {solution.message}; x={solution.x}"
            )
        final_residual = float(np.max(np.abs(solution.fun)))
        wall_residual_scaled[wall_index] = final_residual
        if final_residual > RESIDUAL_TOL:
            warnings.warn(
                f"wall node {wall_index} (x={x_i:.6g}) converged to a scaled "
                f"energy-balance residual of {final_residual:.3e}; the local "
                "wall temperatures may not be trustworthy.",
                RuntimeWarning,
                stacklevel=2,
            )

        T_cw_sol, dT_sol = solution.x
        T_hw_sol = T_cw_sol + dT_sol
        T_hw_arr[wall_index] = T_hw_sol
        T_cw_arr[wall_index] = T_cw_sol
        cold_qdx_arr[wall_index] = physics_helper.dQ_cold_dx(
            x_i,
            T_cw_sol,
            T_cool_i,
            p_cool=p_cool_i,
            quality=quality_i,
        )
        T_hw_guess, T_cw_guess = T_hw_sol, T_cw_sol
        solved_wall_count += 1
        if output is True:
            print(
                f"\rSimulating wall balances: "
                f"{math.ceil(solved_wall_count / n_wall * 100)}%",
                end="",
                flush=True,
            )

    # March the lumped coolant volumes. Every wall node is solved exactly once
    # against the state of the coolant segment to which it is attached.
    for coolant_index in range(n_cool - 1):
        wall_indices = np.flatnonzero(owners == coolant_index)
        for wall_index in wall_indices:
            solve_wall_node(int(wall_index), coolant_index)

        step = float(s_cool[coolant_index + 1] - s_cool[coolant_index])
        local_s = np.clip(
            s_wall[wall_indices],
            s_cool[coolant_index],
            s_cool[coolant_index + 1],
        )
        local_order = np.argsort(local_s)
        local_s = local_s[local_order]
        local_qdx = cold_qdx_arr[wall_indices][local_order]
        if local_s.size == 1:
            weights = np.array([step])
        else:
            edges = np.concatenate(
                (
                    [s_cool[coolant_index]],
                    0.5 * (local_s[:-1] + local_s[1:]),
                    [s_cool[coolant_index + 1]],
                )
            )
            weights = np.diff(edges)
        Q_cold = float(np.dot(local_qdx, weights))

        x_pressure = 0.5 * (x_cool[coolant_index] + x_cool[coolant_index + 1])
        dp_friction = step * physics_helper.coolant_friction_rate(
            x_pressure,
            T_cool_arr[coolant_index],
            p_stagnation_arr[coolant_index],
            quality=quality_arr[coolant_index],
        )

        # The downstream enthalpy is fixed by the heat load, and in the
        # temperature fallback the rise is evaluated on the upstream state, so
        # neither depends on the pressure being solved for below.
        if use_enthalpy:
            coolant_enthalpy_arr[coolant_index + 1] = (
                coolant_enthalpy_arr[coolant_index]
                + physics_helper.coolant_enthalpy_rate(Q_cold)
            )
        else:
            T_next = T_cool_arr[coolant_index] + (
                physics_helper.coolant_temperature_rate(
                    T_cool_arr[coolant_index],
                    p_stagnation_arr[coolant_index],
                    Q_cold,
                )
            )

        # Static pressure obeys dp/dx = -friction - rho u du/dx. The downstream
        # density, velocity and pressure are mutually dependent, so sweep them
        # to a fixed point starting from the friction-only guess.
        rho_here = physics_helper.bulk_density(
            T_cool_arr[coolant_index],
            p_stagnation_arr[coolant_index],
            quality_arr[coolant_index],
        )
        velocity_here = velocity_arr[coolant_index]
        p_next = p_stagnation_arr[coolant_index] + dp_friction

        for _ in range(PRESSURE_SWEEPS):
            if p_next <= 0.0:
                raise RuntimeError(
                    f"coolant pressure fell to {p_next / 1e5:.3g} bar at "
                    f"coolant node {coolant_index + 1} "
                    f"(x={x_cool[coolant_index + 1]:.4g} m) on circuit "
                    f"{circuit.name!r}"
                )
            if use_enthalpy:
                bulk = coolant_state.from_hp(
                    coolant_enthalpy_arr[coolant_index + 1], p_next
                )
            else:
                bulk = coolant_state.from_tp(T_next, p_next)
            rho_there, velocity_there = physics_helper.bulk_velocity(
                x_cool[coolant_index + 1], bulk.T, p_next, bulk.quality
            )
            dp_momentum = _acceleration_pressure_drop(
                rho_here, velocity_here, rho_there, velocity_there
            )
            p_static_next = (
                p_static_arr[coolant_index] + dp_friction - dp_momentum
            )
            p_previous, p_next = (
                p_next,
                p_static_next + 0.5 * rho_there * velocity_there**2,
            )
            if abs(p_next - p_previous) < PRESSURE_SWEEP_TOL:
                break
        else:
            warnings.warn(
                f"coolant pressure did not settle within {PRESSURE_SWEEPS} "
                f"sweeps at coolant node {coolant_index + 1} "
                f"(x={x_cool[coolant_index + 1]:.4g} m) on circuit "
                f"{circuit.name!r}; last change was "
                f"{abs(p_next - p_previous):.3g} Pa.",
                RuntimeWarning,
                stacklevel=2,
            )

        if p_static_next <= 0.0:
            raise RuntimeError(
                f"coolant static pressure fell to "
                f"{p_static_next / 1e5:.3g} bar at "
                f"coolant node {coolant_index + 1} "
                f"(x={x_cool[coolant_index + 1]:.4g} m) on circuit "
                f"{circuit.name!r}"
            )

        p_static_arr[coolant_index + 1] = p_static_next
        p_stagnation_arr[coolant_index + 1] = p_next
        velocity_arr[coolant_index + 1] = velocity_there
        T_cool_arr[coolant_index + 1] = bulk.T
        quality_arr[coolant_index + 1] = bulk.quality
        T_sat_arr[coolant_index + 1] = bulk.T_sat
        phase_arr[coolant_index + 1] = bulk.phase

    if output is True:
        print("\rSimulating wall balances: 100%\n", end="", flush=True)

    # Hot-side quantities are sampled only on their requested output grid.
    T_hw_heat = _interp_axial(x_heat, x_wall, T_hw_arr)
    h_hot_arr = np.zeros(n_heat)
    h_hot_enthalpy_arr = np.zeros(n_heat)
    T_aw_hot_arr = np.zeros(n_heat)
    dQ_dA_arr = np.zeros(n_heat)
    regime_arr = np.empty(n_heat, dtype=object)
    for i, x_i in enumerate(x_heat):
        hot = physics_helper.hot_side_coefficients(x_i, T_hw_heat[i])
        h_hot_arr[i] = hot["h_hot"]
        h_hot_enthalpy_arr[i] = hot["h_g"]
        T_aw_hot_arr[i] = hot["T_aw"]
        dQ_dA_arr[i] = hot["qpp_hot"]
        regime_arr[i] = hot["regime"]

    # Coolant-side diagnostics stay on the small coolant grid. Pressures and
    # velocities already come from the march and are not recomputed here.
    T_cw_cool = _interp_axial(x_cool, x_wall, T_cw_arr)
    h_cold_arr = np.zeros(n_cool)
    for i, x_i in enumerate(x_cool):
        cold = physics_helper.cold_side_coefficients(
            x_i,
            T_cw_cool[i],
            T_cool_arr[i],
            p_cool=p_stagnation_arr[i],
            quality=quality_arr[i],
        )
        h_cold_arr[i] = cold["h_cold"]

    T_stagnation_arr = np.zeros(n_cool)
    for i in range(n_cool):
        if phase_arr[i] == PHASE_TWO_PHASE:
            T_stagnation_arr[i] = T_cool_arr[i]
        else:
            cp = circuit.coolant_transport.get_cp(
                T_cool_arr[i], p_static_arr[i]
            )
            T_stagnation_arr[i] = (
                T_cool_arr[i] + velocity_arr[i] ** 2 / (2.0 * cp)
            )

    # Wall interfaces live on wall nodes; bulk coolant is interpolated solely
    # for the first plotting column and is not used to advance the solution.
    n_walls = len(circuit.walls)
    T_full = np.zeros((n_wall, 1 + n_walls + 1))
    T_full[:, 0] = _interp_axial(x_wall, x_cool, T_cool_arr)
    for i, x_i in enumerate(x_wall):
        interfaces = physics_helper.interface_temperatures(
            x_i, T_hw_arr[i], T_cw_arr[i]
        )
        T_full[i, 1:] = interfaces[::-1]

    boiling = np.flatnonzero(phase_arr == PHASE_TWO_PHASE)
    if boiling.size:
        first = int(boiling[0])
        warnings.warn(
            f"coolant enters the boiling dome at x={x_cool[first]:.4g} m "
            f"(coolant node {first}) and reaches a quality of "
            f"{np.nanmax(quality_arr):.3f} on circuit {circuit.name!r}. "
            "Two-phase heat transfer is not modelled.",
            RuntimeWarning,
            stacklevel=2,
        )

    global_R, final_R = analyse_residuals(residual_log, n_wall)
    minimum_cea_temperature = getattr(
        thrust_chamber.combustion_transport,
        "minimum_cea_temperature",
        None,
    )
    if minimum_cea_temperature is not None:
        extrapolated = np.flatnonzero(T_hw_arr < minimum_cea_temperature)
        if extrapolated.size:
            locations = ", ".join(
                f"{x_wall[index]:.6g}"
                for index in extrapolated[:5]
            )
            if extrapolated.size > 5:
                locations += ", ..."
            warnings.warn(
                f"{extrapolated.size} final hot-wall temperature(s) on "
                f"circuit {circuit.name!r} are below the configured CEA "
                f"temperature boundary of {minimum_cea_temperature:g} K "
                f"(x={locations} m). Gas properties at those wall states "
                "were C1 tangent-extrapolated below the boundary.",
                RuntimeWarning,
                stacklevel=2,
            )
    return RegenResult(
        circuit_name=circuit.name,
        circuit_index=circuit_index,
        x=x_wall,
        T=T_full,
        T_static=T_cool_arr,
        T_stagnation=T_stagnation_arr,
        p_static=p_static_arr,
        p_stagnation=p_stagnation_arr,
        dQ_dA=dQ_dA_arr,
        velocity=velocity_arr,
        h_hot=h_hot_arr,
        h_hot_enthalpy=h_hot_enthalpy_arr,
        h_cold=h_cold_arr,
        T_aw_hot=T_aw_hot_arr,
        residuals=(global_R, final_R),
        qpp_hot=dQ_dA_arr,
        T_drive=T_aw_hot_arr,
        film_regime=regime_arr.astype(str),
        liquid_film=film_model.liquid_results if film_model else None,
        gaseous_film=film_model.gaseous_results if film_model else None,
        coolant_enthalpy=coolant_enthalpy_arr,
        coolant_quality=quality_arr,
        coolant_T_sat=T_sat_arr,
        coolant_phase=phase_arr.astype(str),
        enthalpy_march=use_enthalpy,
        x_wall=x_wall,
        x_heat_flux=x_heat,
        x_coolant=x_cool,
        wall_residual_scaled=wall_residual_scaled,
        wall_converged=wall_residual_scaled <= RESIDUAL_TOL,
    )


def coupled_steady_heating_analysis(
    thrust_chamber,
    boundary_conditions,
    nodes=100,
    circuit_index: int = 0,
    film: Union[bool, str] = "auto",
    film_model: Optional[FilmHeatFluxModel] = None,
    h_liquid_wall: Union[float, Callable[[float], float]] = DEFAULT_H_LIQUID_WALL,
    film_x_array: Optional[Sequence[float]] = None,
    solver: str = "newton",
    output: bool = True,
    heat_curvature_correction: bool = True,
    pressure_curvature_correction: bool = True,
) -> RegenResult:
    """Run a steady heating analysis with film cooling, regen cooling, or both.

    Parameters
    ----------
    thrust_chamber : Any
        Chamber model exposing the required geometry & property APIs.
    boundary_conditions : BoundaryConditions
        Coolant inlet boundary conditions. Also supplies the film injection
        temperature and pressure.
    nodes : int, sequence of sequences, or mapping, optional
        Axial discretization. An integer creates equal uniform wall,
        heat-flux, and coolant grids. An explicit ragged sequence uses the
        order ``[wall, heat_flux, coolant]``; a mapping may instead use those
        names as keys. Different grid lengths implement lumped coolant
        segments connected to multiple wall nodes. Default 100.
    circuit_index : int, optional
        Cooling-circuit index. Default 0.
    film : bool or 'auto', optional
        Whether to include film cooling. ``'auto'`` (default) enables it when
        ``thrust_chamber.film_cooling`` is present. ``True`` requires it.
    film_model : FilmHeatFluxModel, optional
        Pre-solved film boundary condition. Supplying it skips the film solve
        and overrides ``film``.
    h_liquid_wall : float or callable, optional
        Film-to-wall conductance behind the liquid film [W m⁻² K⁻¹]. See the
        module docstring.
    film_x_array : sequence of float, optional
        Grid for the film march. Defaults to the contour's own ``xs``.
    solver : {'newton'}, optional
        Solver selector. Only ``'newton'`` is implemented.
    output : bool, optional
        If True, print progress. Default True.
    heat_curvature_correction : bool, optional
        Apply the RL10 coolant heat-transfer curvature correction. Default
        True.
    pressure_curvature_correction : bool, optional
        Apply the RTE/Ito Darcy-friction curvature correction. Default True.

    Returns
    -------
    RegenResult

    Raises
    ------
    ValueError
        If ``film=True`` but the chamber carries no ``film_cooling``, or the
        solver name is unknown.
    """
    if solver.lower() != "newton":
        raise ValueError(f"solver name not recognized: {solver!r}")

    has_film = getattr(thrust_chamber, "film_cooling", None) is not None

    if film_model is None:
        if film is True and not has_film:
            raise ValueError(
                "film=True but thrust_chamber.film_cooling is None; attach a "
                "FilmCooling object or pass film=False."
            )
        use_film = has_film if film == "auto" else bool(film)
        if use_film:
            if output is True:
                print("Solving film cooling model")
            film_model = FilmHeatFluxModel.from_thrust_chamber(
                thrust_chamber,
                boundary_conditions,
                x_array=film_x_array,
                h_liquid_wall=h_liquid_wall,
            )

    return solve_coupled_heat_exchanger(
        thrust_chamber,
        boundary_conditions,
        nodes,
        circuit_index,
        output,
        film_model=film_model,
        heat_curvature_correction=heat_curvature_correction,
        pressure_curvature_correction=pressure_curvature_correction,
    )


def analyse_residuals(residual_log, n_cells, p=2) -> ResidualResult:
    """Aggregate local solver residuals into global history and final per-cell vector.

    Parameters
    ----------
    residual_log : list or None
        List of tuples ``(cell, iter, R1, R2)`` recorded during solves.
        If ``None`` or empty, returns ``(None, None)``.
    n_cells : int
        Number of axial cells.
    p : int or float, optional
        Norm order for global residual history: ``2`` for RMS,
        ``np.inf`` for :math:`L_\\infty`, etc. Default is 2.

    Returns
    -------
    history : ndarray or None
        Global residual norm for iterations ``0..k_max``, or ``None``.
    final_per_cell : ndarray or None
        Final residual magnitude per cell at its last local iteration, or ``None``.
    """
    if not residual_log:                          # catches [] and None
        return None, None

    log = np.asarray(residual_log)                # shape (m,4)
    cells = log[:, 0].astype(int)
    iters = log[:, 1].astype(int)
    rmag  = np.hypot(log[:, 2], log[:, 3])        # L2 of (R1,R2)

    # ---- global norm history  -----------------------------
    k_max = iters.max()
    history = np.empty(k_max + 1)

    for k in range(k_max + 1):
        mask = iters == k
        if p == np.inf:
            history[k] = rmag[mask].max()
        else:
            history[k] = (rmag[mask] ** p).mean() ** (1.0 / p)

    # ---- per‑cell residual at final local iteration -------
    final_per_cell = np.full(n_cells, np.nan)
    for c in range(n_cells):
        mask = cells == c
        if mask.any():
            final_per_cell[c] = rmag[mask][-1]    # last entry for cell c

    return history, final_per_cell
