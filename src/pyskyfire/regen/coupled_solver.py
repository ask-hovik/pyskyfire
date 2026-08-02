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
the march already tracks and what ``coolant_pressure_rate`` and
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
from dataclasses import dataclass
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
    """

    def __init__(
        self,
        thrust_chamber,
        boundary_conditions,
        circuit_index,
        film_model=None,
        coolant_state: Optional[CoolantState] = None,
    ):
        self.thrust_chamber = thrust_chamber
        self.boundary_conditions = boundary_conditions
        self.circuit_index = circuit_index
        self.film_model = film_model

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
        n_chan = (
            circuit.placement.n_channel_positions
            * circuit.placement.n_channels_per_leaf
        )

        p = self._coolant_pressure(p_cool)
        mdot_c_single_channel = self.boundary_conditions.mdot_coolant / n_chan
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

        R_curv = circuit.radius_of_curvature(x)
        phi_curv = physics.phi_curv(Re_c, D_c, R_curv)

        h_c = physics.h_coolant_colburn(
            k_cf, D_c, Cp_cr, mu_cf, mdot_c_single_channel, A_channel, phi_curv=1
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
        n_channels = (
            circuit.placement.n_channel_positions
            * circuit.placement.n_channels_per_leaf
        )
        return self.boundary_conditions.mdot_coolant / n_channels

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

    def coolant_pressure_rate(self, x, T_cool, p_cool, quality=None):
        """Return static and stagnation coolant-pressure rates [Pa m⁻¹].

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

        dp_stagnation_dx = (
            -friction
            / hydraulic_diameter
            * rho_cool
            * velocity**2
            / 2
            * circuit.ds_dx(x)
        )
        dp_static_dx = (
            dp_stagnation_dx
            - rho_cool * velocity**2 / area * circuit.dA_dx_coolant(x)
        )
        return dp_static_dx, dp_stagnation_dx

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


def solve_coupled_heat_exchanger(
    thrust_chamber,
    boundary_conditions,
    n_nodes,
    circuit_index,
    output,
    film_model: Optional[FilmHeatFluxModel] = None,
    log_residuals: bool = True,
) -> RegenResult:
    """March the 1-D wall energy balance with an optional film boundary condition.

    Parameters
    ----------
    thrust_chamber : Any
        Chamber model exposing geometry and property methods.
    boundary_conditions : BoundaryConditions
        Coolant inlet conditions.
    n_nodes : int
        Number of axial nodes.
    circuit_index : int
        Which cooling circuit to simulate.
    output : bool
        If True, print progress to stdout.
    film_model : FilmHeatFluxModel or None, optional
        Film-derived hot-side boundary condition. ``None`` gives a pure
        regeneratively cooled solve.
    log_residuals : bool, optional
        If True, record local residuals for every Newton iteration per cell.

    Returns
    -------
    RegenResult

    Notes
    -----
    The wall-temperature bracket is derived from the local driving temperature
    rather than the free-stream gas temperature, and reversed heat flow
    (``T_hw < T_cw``) is admitted wherever the drive falls below the bulk
    coolant temperature -- which does happen behind a liquid film once the
    coolant has heated past :math:`T_\\text{sat}`.
    """

    # 1) Build the axial grid in [x_min, x_max].
    residual_log = [] if log_residuals else None
    iter_counter = np.zeros(n_nodes, dtype=int)

    circuit = thrust_chamber.cooling_circuits[circuit_index]
    orig_x_domain = circuit.x_domain
    # Re-interpolate the x_domain to have exactly n_nodes points.
    if circuit.direction == 1:
        x_domain = np.linspace(orig_x_domain[0], orig_x_domain[-1], n_nodes)
    else:
        x_domain = np.linspace(orig_x_domain[-1], orig_x_domain[0], n_nodes)
    dx = abs((x_domain[-1] - x_domain[0]) / (n_nodes - 1))

    # 2) Prepare arrays to hold results
    T_hw_arr = np.zeros(n_nodes)
    T_cw_arr = np.zeros(n_nodes)
    T_cool_arr = np.zeros(n_nodes)
    p_static_arr = np.zeros(n_nodes)
    p_stagnation_arr = np.zeros(n_nodes)
    dQ_dA_arr = np.zeros(n_nodes)
    regime_arr = np.empty(n_nodes, dtype=object)
    h_cool_arr = np.full(n_nodes, np.nan)
    quality_arr = np.full(n_nodes, np.nan)
    T_sat_arr = np.full(n_nodes, np.nan)
    phase_arr = np.empty(n_nodes, dtype=object)

    T_cool_in = boundary_conditions.T_coolant_in
    p_cool_in = boundary_conditions.p_coolant_in
    p_static_arr[0] = p_cool_in
    p_stagnation_arr[0] = p_cool_in  # assuming no velocity at inlet

    # 3) Decide how the coolant state is advanced. Enthalpy where CoolProp can
    #    support it, temperature otherwise -- see coolant_state.probe_coolant.
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
        h_cool_arr[0] = coolant_state.enthalpy(T_cool_in, p_cool_in)
        bulk = coolant_state.from_hp(h_cool_arr[0], p_cool_in)
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
    )

    # Initial guesses for the wall temperatures
    T_hw_guess = 0.5 * (thrust_chamber.combustion_transport.get_T(x_domain[0]) + T_cool_in)
    T_cw_guess = 0.5 * (thrust_chamber.combustion_transport.get_T(x_domain[0]) + T_cool_in)

    if output is True:
        covered = (
            sum(1 for x_i in x_domain if film_model.covers(x_i)) if film_model else 0
        )
        print(
            f"Started coupled heat exchanger simulation with {n_nodes} nodes "
            f"({covered} film-cooled)"
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

    # ==============
    # --- SOLVER ---
    # ==============
    for i in range(n_nodes):
        if output is True:
            print(f'\rSimulating: {math.ceil(i/n_nodes*100)}%', end='', flush=True)

        x_i = x_domain[i]

        # Driving temperature at this station, evaluated at the running wall
        # guess. It sets both the regime label and the solver bracket.
        probe = physics_helper.hot_side_coefficients(x_i, T_hw_guess)
        regime_arr[i] = probe["regime"]
        T_drive_i = probe["T_aw"]
        if not np.isfinite(T_drive_i):
            T_drive_i = thrust_chamber.combustion_transport.get_T(x_i)
        T_cool_i = T_cool_arr[i]

        # Radiation enters as an additive flux, not through the drive
        # temperature, so behind a gaseous film the wall can settle *above*
        # T_aw and still be heated. Bracket on the wall temperature at which
        # the net hot-side flux vanishes, which is the real ceiling. Without
        # this the bracket excludes the root just downstream of dryout, where
        # the film's radiative load switches on.
        q_rad_i = probe.get("q_rad", 0.0)
        h_probe = probe["h_hot"]
        if q_rad_i and np.isfinite(h_probe) and h_probe != 0.0:
            T_drive_i = T_drive_i + q_rad_i / h_probe

        # Fixed residual scale for this node: roughly the heat the wall would
        # take if it sat at the bulk coolant temperature. Holding it constant
        # through the local solve matters -- scaling by max(|Q_hot|, |Q_cond|,
        # |Q_cold|) instead puts a kink in the residual exactly where the
        # balance closes, which stalls the trust-region solve whenever the
        # hot-side conductance is large (e.g. a stiff liquid-film contact).
        dA_dx_i = circuit.dA_dx_thermal_exhaust(x_i)
        if np.isfinite(probe["h_hot"]):
            q_ref = abs(probe["h_hot"]) * max(abs(T_drive_i - T_cool_i), 1.0)
        else:
            q_ref = abs(probe["qpp_hot"])
        Q_ref = max(q_ref * dA_dx_i * dx, 1.0)

        def residuals_scaled(vars_, cell=i):
            T_cw_trial, dT_trial = vars_
            T_hw_trial = T_cw_trial + dT_trial

            # Use per-cell energy flows
            Q_hot_val = physics_helper.dQ_hot_dx(x_i, T_hw_trial) * dx
            Q_cond_val = physics_helper.dQ_cond_dx(x_i, T_hw_trial, T_cw_trial) * dx
            Q_cold_val = physics_helper.dQ_cold_dx(
                x_i,
                T_cw_trial,
                T_cool_arr[i],
                p_cool=p_stagnation_arr[i],
                quality=quality_arr[i],
            ) * dx

            R1 = Q_hot_val - Q_cond_val
            R2 = Q_cond_val - Q_cold_val

            if residual_log is not None:
                k = iter_counter[cell]
                residual_log.append((cell, k, R1, R2))
                iter_counter[cell] += 1

            return np.array([R1 / Q_ref, R2 / Q_ref], dtype=float)

        # Bracket both wall temperatures between the coolant and the drive.
        T_lo = min(T_cool_i, T_drive_i) - 1.0
        T_hi = max(T_cool_i, T_drive_i) + 1.0
        span = max(1.0, T_hi - T_lo)

        # Heat normally flows inward (T_hw >= T_cw). Behind a film whose drive
        # has fallen below the bulk coolant it flows the other way, so only
        # then is a negative wall gradient allowed.
        dT_lo = 0.0 if T_drive_i >= T_cool_i else -span
        dT_hi = span

        # Initial guess in (T_cw, dT) space, clipped into the bracket.
        T_cw_0 = float(np.clip(T_cw_guess, T_lo, T_hi))
        dT_0 = float(np.clip(T_hw_guess - T_cw_guess, dT_lo, dT_hi))
        x0 = np.array([T_cw_0, dT_0], dtype=float)

        res = least_squares(
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
        )

        if not res.success:
            raise RuntimeError(
                f"least_squares failed at node {i} (x={x_i:.6g}): {res.message}; "
                f"x={res.x}"
            )

        # `success` only reports that a termination criterion fired, not that
        # the energy balance actually closed.
        R_final = float(np.max(np.abs(res.fun)))
        if R_final > RESIDUAL_TOL:
            warnings.warn(
                f"node {i} (x={x_i:.6g}) converged to a scaled energy-balance "
                f"residual of {R_final:.3e}; the local wall temperatures may "
                f"not be trustworthy. A coolant film temperature sitting on "
                f"the saturation line is the usual cause -- see the module "
                f"docstring.",
                RuntimeWarning,
                stacklevel=2,
            )

        T_cw_sol, dT_sol = res.x
        T_hw_sol = T_cw_sol + dT_sol

        T_hw_arr[i] = T_hw_sol
        T_cw_arr[i] = T_cw_sol

        # Update guesses for the next node
        T_hw_guess, T_cw_guess = T_hw_sol, T_cw_sol

        # Update the coolant state and pressure for the next node
        if i < n_nodes - 1:
            Q_cold_val = physics_helper.dQ_cold_dx(
                x_i,
                T_cw_sol,
                T_cool_arr[i],
                p_cool=p_stagnation_arr[i],
                quality=quality_arr[i],
            ) * dx
            dp_static, dp_stagnation = physics_helper.coolant_pressure_rate(
                x_i, T_cool_arr[i], p_stagnation_arr[i], quality=quality_arr[i]
            )
            p_static_arr[i + 1] = p_static_arr[i] + dp_static * dx
            p_stagnation_arr[i + 1] = p_stagnation_arr[i] + dp_stagnation * dx

            # A non-positive pressure is always a modelling error, not a
            # solution: the circuit cannot pass this mass flow. Catch it here,
            # because the property flash downstream would otherwise fail with a
            # CoolProp message that says nothing about the cause.
            if p_stagnation_arr[i + 1] <= 0.0:
                raise RuntimeError(
                    f"coolant stagnation pressure fell to "
                    f"{p_stagnation_arr[i + 1]/1e5:.3g} bar at node {i + 1} "
                    f"(x = {x_domain[i + 1]:.4g} m) on circuit {circuit.name!r}, "
                    f"starting from {p_cool_in/1e5:.3g} bar. The channel cannot "
                    f"pass {boundary_conditions.mdot_coolant:.4g} kg/s: raise the "
                    f"inlet pressure, enlarge the channel, or lower the flow."
                )

            # Advance the state variable, then resolve (T, quality) from it at
            # the *new* pressure -- the dome moves with pressure, so resolving
            # against the upstream pressure would misplace the crossing.
            if use_enthalpy:
                h_cool_arr[i + 1] = h_cool_arr[i] + physics_helper.coolant_enthalpy_rate(
                    Q_cold_val
                )
                bulk = coolant_state.from_hp(h_cool_arr[i + 1], p_stagnation_arr[i + 1])
            else:
                dT = physics_helper.coolant_temperature_rate(
                    T_cool_arr[i], p_stagnation_arr[i], Q_cold_val
                )
                bulk = coolant_state.from_tp(
                    T_cool_arr[i] + dT, p_stagnation_arr[i + 1]
                )

            T_cool_arr[i + 1] = bulk.T
            quality_arr[i + 1] = bulk.quality
            T_sat_arr[i + 1] = bulk.T_sat
            phase_arr[i + 1] = bulk.phase

    if output is True:
        print(f'\rSimulating: {100}%\n', end='', flush=True)

    # ===============
    # --- RESULTS ---
    # ===============

    # 1. Hot-side coefficients and heat flux at the converged wall temperatures.
    h_hot_arr = np.zeros(n_nodes)           # effective hot-side h [W/m²/K]
    h_hot_enthalpy_arr = np.zeros(n_nodes)  # enthalpy-based hot-side coeff [kg/m²/s]
    h_cold_arr = np.zeros(n_nodes)          # coolant-side h [W/m²/K]
    T_aw_hot_arr = np.zeros(n_nodes)        # driving temperature [K]
    for i, x_i in enumerate(x_domain):
        hot = physics_helper.hot_side_coefficients(x_i, T_hw_arr[i])
        cold = physics_helper.cold_side_coefficients(
            x_i,
            T_cw_arr[i],
            T_cool_arr[i],
            p_cool=p_stagnation_arr[i],
            quality=quality_arr[i],
        )

        h_hot_arr[i] = hot["h_hot"]
        h_hot_enthalpy_arr[i] = hot["h_g"]
        h_cold_arr[i] = cold["h_cold"]
        T_aw_hot_arr[i] = hot["T_aw"]
        regime_arr[i] = hot["regime"]
        dQ_dA_arr[i] = hot["qpp_hot"]

    # 2. Coolant velocity at each node
    velocity_arr = np.zeros(n_nodes)
    n_chan = (
        thrust_chamber.cooling_circuits[circuit_index].placement.n_channel_positions
        * thrust_chamber.cooling_circuits[circuit_index].placement.n_channels_per_leaf
    )
    for i, x_i in enumerate(x_domain):
        A_channel = thrust_chamber.cooling_circuits[circuit_index].A_coolant(x_i)
        mdot_c_single = boundary_conditions.mdot_coolant / n_chan
        rho_cool = physics_helper.bulk_density(
            T_cool_arr[i], p_static_arr[i], quality_arr[i]
        )
        velocity_arr[i] = physics.u_coolant(rho_cool, mdot_c_single, A_channel)

    T_stagnation_arr = np.zeros_like(T_cool_arr)
    for i, x_i in enumerate(x_domain):
        # Inside the dome there is no stagnation *temperature* rise to speak
        # of: bringing a saturated mixture to rest condenses vapour at fixed
        # T_sat rather than heating it, and cp is not even finite there.
        if phase_arr[i] == PHASE_TWO_PHASE:
            T_stagnation_arr[i] = T_cool_arr[i]
            continue
        cp_cool = thrust_chamber.cooling_circuits[circuit_index].coolant_transport.get_cp(
            T_cool_arr[i], p_static_arr[i]
        )
        T_stagnation_arr[i] = T_cool_arr[i] + (velocity_arr[i] ** 2) / (2.0 * cp_cool)

    n_walls = len(thrust_chamber.cooling_circuits[circuit_index].walls)
    T_full = np.zeros((n_nodes, 1 + n_walls + 1))  # coolant + (n_walls+1) interfaces

    for i, x_i in enumerate(x_domain):
        # Ts = [T_hot, ..., T_cold], length = n_walls+1
        Ts = physics_helper.interface_temperatures(x_i, T_hw_arr[i], T_cw_arr[i])
        Ts_rev = Ts[::-1]  # reverse so Ts[0]=T_cold, Ts[-1]=T_hot
        T_full[i, 0] = T_cool_arr[i]
        T_full[i, 1:] = Ts_rev

    p_static_corrected = np.zeros_like(p_stagnation_arr)
    for i in range(n_nodes):
        rho_i = physics_helper.bulk_density(
            T_cool_arr[i], p_stagnation_arr[i], quality_arr[i]
        )
        q_dyn = 0.5 * rho_i * velocity_arr[i] ** 2
        p_static_corrected[i] = p_stagnation_arr[i] - q_dyn

    p_static_arr = p_static_corrected

    # The quality profile says where the coolant boils; it does not say how hot
    # the wall gets there, because h_cold is still a single-phase correlation.
    boiling = np.flatnonzero(phase_arr == PHASE_TWO_PHASE)
    if boiling.size:
        i0 = int(boiling[0])
        warnings.warn(
            f"coolant enters the boiling dome at x = {x_domain[i0]:.4g} m "
            f"(node {i0}) and reaches a quality of "
            f"{np.nanmax(quality_arr):.3f} on circuit {circuit.name!r}. "
            f"Two-phase heat transfer is not modelled -- h_cold remains the "
            f"single-phase correlation on saturated-liquid properties, so wall "
            f"temperatures downstream of that station are not trustworthy.",
            RuntimeWarning,
            stacklevel=2,
        )
        if output is True:
            print(
                f"Coolant boils from x = {x_domain[i0]:.4g} m; exit quality "
                f"{quality_arr[-1]:.3f} (two-phase h_cold not modelled)"
            )

    global_R, final_R = analyse_residuals(residual_log, n_nodes)

    return RegenResult(
        circuit_name=thrust_chamber.cooling_circuits[circuit_index].name,
        circuit_index=circuit_index,
        x=x_domain,
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
        coolant_enthalpy=h_cool_arr,
        coolant_quality=quality_arr,
        coolant_T_sat=T_sat_arr,
        coolant_phase=phase_arr.astype(str),
        enthalpy_march=use_enthalpy,
    )


def coupled_steady_heating_analysis(
    thrust_chamber,
    boundary_conditions,
    n_nodes: int = 100,
    circuit_index: int = 0,
    film: Union[bool, str] = "auto",
    film_model: Optional[FilmHeatFluxModel] = None,
    h_liquid_wall: Union[float, Callable[[float], float]] = DEFAULT_H_LIQUID_WALL,
    film_x_array: Optional[Sequence[float]] = None,
    solver: str = "newton",
    output: bool = True,
) -> RegenResult:
    """Run a steady heating analysis with film cooling, regen cooling, or both.

    Parameters
    ----------
    thrust_chamber : Any
        Chamber model exposing the required geometry & property APIs.
    boundary_conditions : BoundaryConditions
        Coolant inlet boundary conditions. Also supplies the film injection
        temperature and pressure.
    n_nodes : int, optional
        Number of axial nodes for the regenerative march. Default 100.
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
        n_nodes,
        circuit_index,
        output,
        film_model=film_model,
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
