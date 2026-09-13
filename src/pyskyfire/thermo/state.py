"""CoolProp-backed fluid state resolution, including the saturation dome.

Specific enthalpy and pressure identify states through the boiling dome,
where temperature and pressure alone cannot determine vapour quality.
``FluidStateResolver`` caches saturation data and resolves bulk states using
homogeneous (no-slip) two-phase density. ``probe_fluid`` reports whether the
backend supports enthalpy-based resolution and saturation properties.

This module supplies thermodynamic states, not heat-transfer or pressure-drop
correlations. Mixtures and unsupported backends retain the temperature-based
fallback; above the critical pressure, quality is undefined.
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
    "FluidCapability",
    "SaturationState",
    "FluidState",
    "FluidStateResolver",
    "probe_fluid",
]


#: Bulk fluid is subcooled liquid.
PHASE_LIQUID = "liquid"
#: Bulk fluid is inside the boiling dome; ``quality`` is meaningful.
PHASE_TWO_PHASE = "two_phase"
#: Bulk fluid is superheated vapour.
PHASE_VAPOUR = "vapour"
#: Local pressure is above the critical pressure; no dome exists.
PHASE_SUPERCRITICAL = "supercritical"
#: No saturation data available; the march is running single-phase blind.
PHASE_UNKNOWN = "unknown"

#: Temperature step used to back away from CoolProp's saturation guard band [K].
_SATURATION_NUDGE = 1.0e-3


@dataclass(frozen=True)
class FluidCapability:
    """What two-phase machinery is usable for a given fluid.

    Attributes
    ----------
    enthalpy_march : bool
        Whether the fluid supports the ``(h, p)`` march. ``False`` sends the
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
    which may evaluate the fluid side a couple of hundred times per node --
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
class FluidState:
    """Resolved bulk fluid state at one axial station.

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


def probe_fluid(fluid: str, p_ref: float) -> FluidCapability:
    """Determine what two-phase support a CoolProp fluid string has.

    Run once per solve, at the inlet pressure. Cheap: a handful of property
    calls.

    Parameters
    ----------
    fluid : str
        CoolProp fluid string, e.g. ``'HEOS::Ethanol[1.000]'``.
    p_ref : float
        Reference pressure to probe at [Pa], normally the fluid inlet.

    Returns
    -------
    FluidCapability

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
        return FluidCapability(
            enthalpy_march=False,
            saturation=False,
            is_mixture=True,
            backend=backend,
            reason=(
                "fluid is a mixture; the (h, p) flash CoolProp needs for an "
                "enthalpy march costs ~6500x a pure-fluid call, which is not "
                "usable in a nodal march. Falling back to the temperature "
                "march -- boiling is not modelled"
            ),
        )

    if backend.upper() in ("INCOMP", "IF97"):
        return FluidCapability(
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

    return FluidCapability(
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


class FluidStateResolver:
    """Resolve fluid states from ``(h, p)``, with cached saturation data.

    Parameters
    ----------
    transport : Any
        Object carrying a CoolProp fluid string on ``.fluid``, normally a
        :class:`~pyskyfire.thermo.transport.FluidTransport`.
    p_ref : float
        Reference pressure for the capability probe [Pa], normally the fluid
        inlet pressure.

    Attributes
    ----------
    capability : FluidCapability
        Result of :func:`probe_fluid`. Check ``enthalpy_march`` before using
        :meth:`from_hp`.
    """

    def __init__(self, transport, p_ref: float):
        self.transport = transport
        self.fluid = getattr(transport, "fluid", None)
        if self.fluid is None:
            raise TypeError(
                "transport does not expose a CoolProp fluid string on "
                "'.fluid'; cannot resolve two-phase fluid states."
            )
        self.capability = probe_fluid(self.fluid, float(p_ref))
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

    def from_hp(self, h: float, p: float) -> FluidState:
        """Resolve the bulk state from specific enthalpy and pressure.

        Parameters
        ----------
        h : float
            Specific enthalpy [J kg⁻¹].
        p : float
            Static pressure [Pa].

        Returns
        -------
        FluidState

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
            return FluidState(
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
            return FluidState(
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
        return FluidState(
            T=T,
            p=p,
            h=h,
            quality=math.nan,
            rho=rho,
            phase=PHASE_LIQUID if h < sat.h_f else PHASE_VAPOUR,
            T_sat=sat.T_sat,
        )

    def from_tp(self, T: float, p: float) -> FluidState:
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
        return FluidState(
            T=T,
            p=p,
            h=math.nan,
            quality=math.nan,
            rho=rho,
            phase=PHASE_SUPERCRITICAL if supercritical else PHASE_UNKNOWN,
            T_sat=sat.T_sat if sat is not None else math.nan,
        )
