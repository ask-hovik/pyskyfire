"""Coolant thermodynamic state for the regenerative march, including the boiling dome.

Why the march is on enthalpy, not temperature
---------------------------------------------

The obvious way to advance the coolant along a channel is
:math:`\\Delta T = Q / (\\dot m\\, c_p)`, and that is what
:meth:`~pyskyfire.regen.coupled_solver.CoupledHeatExchangerPhysics.coolant_temperature_rate`
does. It is only valid while the coolant is single-phase. Once the bulk
reaches saturation, added heat goes into latent heat at constant temperature,
and temperature stops being a state variable at all: inside the dome ``T`` and
``p`` are not independent, so ``(T, p)`` no longer identifies a state.

That failure is silent rather than loud. CoolProp answers a ``(T, p)`` query
just inside the dome by picking a single-phase branch, so a temperature march
that steps past :math:`T_\\text{sat}` lands on *superheated vapour*. For
ethanol at 20 bar, 0.2 K either side of saturation:

===================  ===========================
``T_sat - 0.1 K``    :math:`\\rho = 596.9` kg m⁻³
``T_sat + 0.1 K``    :math:`\\rho = 31.4` kg m⁻³
===================  ===========================

A 19x density jump, and with it a 19x jump in ``rho u^2`` and therefore in the
pressure gradient. No exception is raised; the march simply continues on the
wrong branch and returns a plausible-looking profile.

Specific enthalpy has none of these problems. It is monotonic in the heat
added, it stays a valid state variable through the dome, and ``(h, p)``
identifies the state everywhere -- subcooled, saturated, superheated, or
supercritical. So the march integrates

.. math:: h_{i+1} = h_i + Q_\\text{cold} / \\dot m_\\text{channel}

and :meth:`CoolantState.from_hp` resolves ``(T, quality, rho)`` afterwards.

What this module does *not* model
---------------------------------

**Two-phase heat transfer.** In the dome the coolant-side coefficient is still
the single-phase Colburn correlation, evaluated with *saturated liquid*
properties. That is a placeholder in both directions: it misses the flow-boiling
enhancement (which would raise it, by 2-10x) and it misses dryout (which would
collapse it). Treat the quality profile as telling you *where* the coolant
boils, not how hot the wall gets there. Nucleate boiling and CHF are not
modelled at all.

**Subcooled boiling.** Where the bulk is liquid but the wall is above
:math:`T_\\text{sat}`, real channels boil at the wall. That is not modelled
either. What *is* handled is the property-lookup artefact it used to cause: the
Colburn film temperature :math:`(T_\\text{cool} + T_\\text{cw})/2` would cross
saturation and CoolProp would return *vapour* properties for what is physically
superheated liquid, putting a spurious cliff in ``h_cold``. Film properties are
now clamped to saturated liquid instead (see
:meth:`CoolantState.film_properties`). Still an approximation, but no longer a
discontinuity, and wrong in the conservative direction.

**Two-phase pressure drop.** Density in the dome is the homogeneous (no-slip)
value :math:`1/(x/\\rho_g + (1-x)/\\rho_f)`. This captures the large
acceleration effect but not the two-phase friction multiplier, so channel
:math:`\\Delta p` is under-predicted once quality is appreciable.

Where the enthalpy march is unavailable
---------------------------------------

:func:`probe_coolant` decides per fluid, once, at solver start:

* **Mixtures** (``HEOS::A[..]&B[..]``) fall back to the temperature march. Not
  because the physics is different but because it is computationally
  impossible: a mixture ``(h, p)`` flash costs ~374 ms against ~57 us for a
  pure fluid, a factor of 6500, which would put a 100-node run into the hours.
  Mixtures also have no ``pcrit`` through ``Props1SI`` and no surface tension,
  so most of the two-phase machinery has nothing to stand on. Revisit when
  tabulated backends land.
* **Incompressible backends** (``INCOMP::``) have no saturation data.
* **Supercritical operation** (:math:`p > p_\\text{crit}`) keeps the enthalpy
  march -- it is still the better-conditioned variable -- but there is no dome,
  so quality is undefined and reported as ``nan``. This is the normal state of
  affairs for hydrogen (:math:`p_\\text{crit}` = 13.0 bar) and for kerosene-like
  fluids (n-dodecane, 18.2 bar) at any realistic channel pressure. Heat-transfer
  deterioration near the pseudo-critical line is not modelled.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import CoolProp.CoolProp as CP

__all__ = [
    "PHASE_LIQUID",
    "PHASE_TWO_PHASE",
    "PHASE_VAPOUR",
    "PHASE_SUPERCRITICAL",
    "PHASE_UNKNOWN",
    "CoolantCapability",
    "SaturationState",
    "BulkState",
    "CoolantState",
    "probe_coolant",
]


#: Bulk coolant is subcooled liquid.
PHASE_LIQUID = "liquid"
#: Bulk coolant is inside the boiling dome; ``quality`` is meaningful.
PHASE_TWO_PHASE = "two_phase"
#: Bulk coolant is superheated vapour.
PHASE_VAPOUR = "vapour"
#: Local pressure is above the critical pressure; no dome exists.
PHASE_SUPERCRITICAL = "supercritical"
#: No saturation data available; the march is running single-phase blind.
PHASE_UNKNOWN = "unknown"

#: Temperature step used to back away from CoolProp's saturation guard band [K].
_SATURATION_NUDGE = 1.0e-3


@dataclass(frozen=True)
class CoolantCapability:
    """What two-phase machinery is usable for a given coolant.

    Attributes
    ----------
    enthalpy_march : bool
        Whether the coolant supports the ``(h, p)`` march. ``False`` sends the
        solver back to the temperature march.
    saturation : bool
        Whether saturation properties (and therefore quality) are available.
    is_mixture : bool
        Whether the CoolProp string names more than one component.
    backend : str
        CoolProp backend prefix, e.g. ``'HEOS'`` or ``'INCOMP'``.
    p_crit, T_crit : float or None
        Critical point [Pa], [K]. ``None`` when it could not be determined.
    reason : str or None
        Human-readable explanation when something is unavailable.
    """

    enthalpy_march: bool
    saturation: bool
    is_mixture: bool
    backend: str
    p_crit: Optional[float] = None
    T_crit: Optional[float] = None
    reason: Optional[str] = None


@dataclass(frozen=True)
class SaturationState:
    """Saturation properties at one pressure.

    Built once per axial station and cached, so that the inner wall solve --
    which may evaluate the coolant side a couple of hundred times per node --
    never triggers a flash.

    Attributes
    ----------
    p : float
        Pressure [Pa].
    T_sat : float
        Saturation temperature [K].
    h_f, h_g : float
        Saturated liquid and vapour specific enthalpy [J kg⁻¹].
    rho_f, rho_g : float
        Saturated liquid and vapour density [kg m⁻³].
    mu_f, k_f, cp_f : float
        Saturated liquid viscosity [Pa s], conductivity [W m⁻¹ K⁻¹] and
        specific heat [J kg⁻¹ K⁻¹].
    """

    p: float
    T_sat: float
    h_f: float
    h_g: float
    rho_f: float
    rho_g: float
    mu_f: float
    k_f: float
    cp_f: float

    @property
    def h_fg(self) -> float:
        """Latent heat of vaporisation [J kg⁻¹]."""
        return self.h_g - self.h_f

    def homogeneous_density(self, quality: float) -> float:
        """Homogeneous (no-slip) two-phase density [kg m⁻³].

        Parameters
        ----------
        quality : float
            Vapour mass fraction, clipped into ``[0, 1]``.
        """
        x = min(max(float(quality), 0.0), 1.0)
        return 1.0 / (x / self.rho_g + (1.0 - x) / self.rho_f)


@dataclass(frozen=True)
class BulkState:
    """Resolved bulk coolant state at one axial station.

    Attributes
    ----------
    T : float
        Bulk temperature [K]. Equal to ``T_sat`` inside the dome.
    p : float
        Pressure the state was resolved at [Pa].
    h : float
        Specific enthalpy [J kg⁻¹]. ``nan`` under the temperature march.
    quality : float
        Vapour mass fraction, ``nan`` outside the dome (including subcooled
        liquid, where it is *not* zero but undefined-by-convention here --
        see :attr:`phase` for the distinction).
    rho : float
        Density [kg m⁻³], homogeneous inside the dome.
    phase : str
        One of the ``PHASE_*`` constants.
    T_sat : float
        Local saturation temperature [K], ``nan`` when unavailable.
    """

    T: float
    p: float
    h: float
    quality: float
    rho: float
    phase: str
    T_sat: float

    @property
    def is_two_phase(self) -> bool:
        """True when the bulk sits inside the boiling dome."""
        return self.phase == PHASE_TWO_PHASE


def probe_coolant(fluid: str, p_ref: float) -> CoolantCapability:
    """Determine what two-phase support a CoolProp fluid string has.

    Run once per solve, at the inlet pressure. Cheap: a handful of property
    calls.

    Parameters
    ----------
    fluid : str
        CoolProp fluid string, e.g. ``'HEOS::Ethanol[1.000]'``.
    p_ref : float
        Reference pressure to probe at [Pa], normally the coolant inlet.

    Returns
    -------
    CoolantCapability

    Notes
    -----
    The critical-pressure probe is deliberately tolerant. ``Props1SI`` refuses
    mixtures outright, and ``AbstractState.p_critical()`` can fail even for
    well-behaved binaries (ethanol/water reports "found 4 critical points"), so
    a missing critical point downgrades the capability rather than raising.
    """
    backend = fluid.split("::")[0] if "::" in fluid else "HEOS"
    is_mixture = "&" in fluid

    if is_mixture:
        return CoolantCapability(
            enthalpy_march=False,
            saturation=False,
            is_mixture=True,
            backend=backend,
            reason=(
                "coolant is a mixture; the (h, p) flash CoolProp needs for an "
                "enthalpy march costs ~6500x a pure-fluid call, which is not "
                "usable in a nodal march. Falling back to the temperature "
                "march -- boiling is not modelled"
            ),
        )

    if backend.upper() in ("INCOMP", "IF97"):
        return CoolantCapability(
            enthalpy_march=False,
            saturation=False,
            is_mixture=False,
            backend=backend,
            reason=(
                f"the {backend} backend carries no saturation data. Falling "
                f"back to the temperature march -- boiling is not modelled"
            ),
        )

    p_crit = _try(lambda: CP.PropsSI("pcrit", fluid))
    T_crit = _try(lambda: CP.PropsSI("Tcrit", fluid))
    supercritical_at_ref = p_crit is not None and p_ref >= p_crit

    # Whether the *fluid* has a usable dome is a different question from
    # whether it has one at p_ref, so probe saturation below the critical
    # pressure when the reference pressure is above it.
    p_sat_probe = 0.5 * p_crit if supercritical_at_ref else p_ref
    sat_ok = _try(lambda: CP.PropsSI("T", "P", p_sat_probe, "Q", 0.0, fluid)) is not None

    # Reference enthalpy for the flash round-trip. A (T, p) call is the wrong
    # probe here: at p_ref the interesting temperature is T_sat, and CoolProp
    # refuses a (T, p) pair sitting on the saturation line. Take h_f straight
    # from a P-Q flash instead, and fall back to the critical isotherm when
    # there is no dome at p_ref.
    if supercritical_at_ref:
        h_ref = (
            _try(lambda: CP.PropsSI("Hmass", "T", T_crit, "P", p_ref, fluid))
            if T_crit is not None
            else None
        )
    else:
        h_ref = _try(lambda: CP.PropsSI("Hmass", "P", p_ref, "Q", 0.0, fluid))

    march_ok = (
        h_ref is not None
        and _try(lambda: CP.PropsSI("T", "Hmass", h_ref, "P", p_ref, fluid)) is not None
    )

    if not march_ok:
        reason = (
            f"CoolProp could not run an (h, p) flash for {fluid!r}. Falling "
            f"back to the temperature march -- boiling is not modelled"
        )
    elif not sat_ok:
        reason = (
            f"CoolProp could not produce saturation properties for {fluid!r}; "
            f"the enthalpy march still runs but quality is not reported"
        )
    else:
        reason = None

    return CoolantCapability(
        enthalpy_march=march_ok,
        saturation=sat_ok,
        is_mixture=False,
        backend=backend,
        p_crit=p_crit,
        T_crit=T_crit,
        reason=reason,
    )


def _try(fn):
    """Return ``fn()``, or ``None`` if CoolProp refuses."""
    try:
        return fn()
    except Exception:  # noqa: BLE001
        return None


class CoolantState:
    """Resolve coolant states from ``(h, p)``, with cached saturation data.

    Parameters
    ----------
    coolant_transport : Any
        Object carrying a CoolProp fluid string on ``.fluid``, normally a
        :class:`~pyskyfire.skycea.coolant_transport.CoolantTransport`.
    p_ref : float
        Reference pressure for the capability probe [Pa], normally the coolant
        inlet pressure.

    Attributes
    ----------
    capability : CoolantCapability
        Result of :func:`probe_coolant`. Check ``enthalpy_march`` before using
        :meth:`from_hp`.
    """

    def __init__(self, coolant_transport, p_ref: float):
        self.transport = coolant_transport
        self.fluid = getattr(coolant_transport, "fluid", None)
        if self.fluid is None:
            raise TypeError(
                "coolant_transport does not expose a CoolProp fluid string on "
                "'.fluid'; cannot resolve two-phase coolant states."
            )
        self.capability = probe_coolant(self.fluid, float(p_ref))
        self._sat_cache: dict[float, Optional[SaturationState]] = {}

    # ------------------------------------------------------------------
    # saturation
    # ------------------------------------------------------------------
    def saturation(self, p: float) -> Optional[SaturationState]:
        """Saturation state at ``p`` [Pa], or ``None`` where there is no dome.

        ``None`` means either that the fluid has no usable saturation data or
        that ``p`` is at/above the critical pressure. Results are cached on
        ``p``, so the inner wall solve costs nothing after the first call at a
        station.
        """
        p = float(p)
        if p in self._sat_cache:
            return self._sat_cache[p]

        state = self._build_saturation(p)
        self._sat_cache[p] = state
        return state

    def _build_saturation(self, p: float) -> Optional[SaturationState]:
        cap = self.capability
        if not cap.saturation:
            return None
        if cap.p_crit is not None and p >= cap.p_crit:
            return None

        try:
            T_sat = CP.PropsSI("T", "P", p, "Q", 0, self.fluid)
            h_f = CP.PropsSI("Hmass", "P", p, "Q", 0, self.fluid)
            h_g = CP.PropsSI("Hmass", "P", p, "Q", 1, self.fluid)
            rho_f = CP.PropsSI("Dmass", "P", p, "Q", 0, self.fluid)
            rho_g = CP.PropsSI("Dmass", "P", p, "Q", 1, self.fluid)
        except Exception:  # noqa: BLE001 - above p_crit, or outside the EOS range
            return None

        # Transport properties are not defined *on* the saturation line for
        # every fluid, so fall back a millikelvin into the liquid.
        mu_f = self._liquid_transport("viscosity", T_sat, p)
        k_f = self._liquid_transport("conductivity", T_sat, p)
        cp_f = self._liquid_transport("Cpmass", T_sat, p)
        if None in (mu_f, k_f, cp_f):
            return None

        return SaturationState(
            p=p,
            T_sat=T_sat,
            h_f=h_f,
            h_g=h_g,
            rho_f=rho_f,
            rho_g=rho_g,
            mu_f=mu_f,
            k_f=k_f,
            cp_f=cp_f,
        )

    def _liquid_transport(self, output: str, T_sat: float, p: float):
        value = _try(lambda: CP.PropsSI(output, "P", p, "Q", 0, self.fluid))
        if value is None:
            value = _try(
                lambda: CP.PropsSI(output, "T", T_sat - _SATURATION_NUDGE, "P", p, self.fluid)
            )
        return value

    # ------------------------------------------------------------------
    # state resolution
    # ------------------------------------------------------------------
    def enthalpy(self, T: float, p: float) -> float:
        """Specific enthalpy at ``(T, p)`` [J kg⁻¹].

        Used once, to turn the inlet boundary condition into the march's state
        variable.
        """
        return float(CP.PropsSI("Hmass", "T", float(T), "P", float(p), self.fluid))

    def from_hp(self, h: float, p: float) -> BulkState:
        """Resolve the bulk state from specific enthalpy and pressure.

        Parameters
        ----------
        h : float
            Specific enthalpy [J kg⁻¹].
        p : float
            Static pressure [Pa].

        Returns
        -------
        BulkState

        Notes
        -----
        Inside the dome the temperature is taken as ``T_sat`` and the density as
        the homogeneous mixture value, rather than trusting a flash to land on
        the right branch. Outside it, both come from an ``(h, p)`` flash, which
        is unambiguous on either side.
        """
        h = float(h)
        p = float(p)
        sat = self.saturation(p)

        if sat is None:
            T = float(CP.PropsSI("T", "Hmass", h, "P", p, self.fluid))
            rho = float(CP.PropsSI("Dmass", "Hmass", h, "P", p, self.fluid))
            cap = self.capability
            supercritical = cap.p_crit is not None and p >= cap.p_crit
            return BulkState(
                T=T,
                p=p,
                h=h,
                quality=math.nan,
                rho=rho,
                phase=PHASE_SUPERCRITICAL if supercritical else PHASE_UNKNOWN,
                T_sat=math.nan,
            )

        if sat.h_f <= h <= sat.h_g:
            quality = (h - sat.h_f) / sat.h_fg
            return BulkState(
                T=sat.T_sat,
                p=p,
                h=h,
                quality=quality,
                rho=sat.homogeneous_density(quality),
                phase=PHASE_TWO_PHASE,
                T_sat=sat.T_sat,
            )

        T = float(CP.PropsSI("T", "Hmass", h, "P", p, self.fluid))
        rho = float(CP.PropsSI("Dmass", "Hmass", h, "P", p, self.fluid))
        return BulkState(
            T=T,
            p=p,
            h=h,
            quality=math.nan,
            rho=rho,
            phase=PHASE_LIQUID if h < sat.h_f else PHASE_VAPOUR,
            T_sat=sat.T_sat,
        )

    def from_tp(self, T: float, p: float) -> BulkState:
        """Resolve the bulk state from ``(T, p)``, for the temperature march.

        Cannot represent a two-phase state -- that is the whole limitation of
        the temperature march -- so the phase is reported as
        :data:`PHASE_UNKNOWN` unless the fluid is known to be supercritical.
        """
        T = float(T)
        p = float(p)
        rho = float(CP.PropsSI("Dmass", "T", T, "P", p, self.fluid))
        cap = self.capability
        supercritical = cap.p_crit is not None and p >= cap.p_crit
        sat = self.saturation(p)
        return BulkState(
            T=T,
            p=p,
            h=math.nan,
            quality=math.nan,
            rho=rho,
            phase=PHASE_SUPERCRITICAL if supercritical else PHASE_UNKNOWN,
            T_sat=sat.T_sat if sat is not None else math.nan,
        )

    # ------------------------------------------------------------------
    # properties for the coolant-side correlation
    # ------------------------------------------------------------------
    def film_properties(self, T_film: float, p: float, getters) -> tuple:
        """Coolant film properties for the Colburn correlation.

        Parameters
        ----------
        T_film : float
            Film temperature :math:`(T_\\text{cool} + T_\\text{cw})/2` [K].
        p : float
            Coolant pressure [Pa].
        getters : sequence of callable
            Property getters taking ``(T, p)``, in the order
            ``(conductivity, cp, viscosity)``.

        Returns
        -------
        tuple
            ``(k, cp, mu)`` at the film condition.

        Notes
        -----
        The film temperature routinely exceeds :math:`T_\\text{sat}` while the
        bulk is still subcooled liquid, because the wall is hot. A raw ``(T, p)``
        lookup there returns *vapour* properties, which is physically wrong --
        the fluid at the wall is superheated liquid on the verge of nucleating,
        not vapour -- and puts a cliff in ``h_cold`` that the local wall solve
        can fail to cross. Clamping to saturated liquid removes the cliff. It is
        still not a boiling model: the real coefficient in that region is higher
        than this returns.
        """
        k_get, cp_get, mu_get = getters
        sat = self.saturation(p)

        if sat is not None and T_film >= sat.T_sat:
            return sat.k_f, sat.cp_f, sat.mu_f

        return (
            _off_saturation(k_get, T_film, p),
            _off_saturation(cp_get, T_film, p),
            _off_saturation(mu_get, T_film, p),
        )


def _off_saturation(getter, T, p):
    """Evaluate a property, stepping off CoolProp's saturation guard band.

    CoolProp refuses a ``(T, p)`` pair within 1e-4 % of the saturation line
    because both branches are valid there and it will not choose. Retry a
    millikelvin to the liquid side. This is a crash guard, not a model.
    """
    try:
        return getter(T, p)
    except ValueError:
        return getter(T - _SATURATION_NUDGE, p)
