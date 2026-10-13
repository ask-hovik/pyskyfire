"""Equations of state: the thermodynamic relations a design reads a fluid through.

A turbomachine design needs more of a fluid than its density at one point. It
needs the states at either end of the machine, and the state an isentropic
process reaches at the other end, since that fixes the head. A turbine's
meanline also needs the states between its rows, which it knows by enthalpy
and by the entropy the losses so far have added, and the speed of sound at
each. An :class:`EquationOfState` answers those questions and returns the
package's :class:`~pyskyfire.thermo.state.FluidState` each time, so the code
asking does not care whether the answer came from CoolProp or from a constant.

:class:`Incompressible` is a constant density. Giving a design a density
instead of a medium means this. :class:`CoolPropEOS` is a real fluid, resolved
through the same :class:`~pyskyfire.thermo.state.FluidStateResolver` the
cooling march uses, so the saturation dome is handled the same way here as
there.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from types import SimpleNamespace

import CoolProp.CoolProp as CP

from .fluid import Fluid
from .state import (
    PHASE_SUPERCRITICAL,
    PHASE_UNKNOWN,
    FluidState,
    FluidStateResolver,
    _SATURATION_NUDGE,
)

__all__ = ["EquationOfState", "Incompressible", "CoolPropEOS", "as_equation_of_state"]


class EquationOfState(ABC):
    """The state lookups a machine design asks of a fluid.

    Every pressure is in Pa, temperature in K and specific enthalpy in
    J kg⁻¹. Enthalpy is only ever compared with enthalpy from the same
    equation of state, so its datum does not matter.
    """

    @abstractmethod
    def state_pT(self, p: float, T: float) -> FluidState:
        """The state at pressure ``p`` and temperature ``T``."""

    @abstractmethod
    def state_ph(self, p: float, h: float) -> FluidState:
        """The state at pressure ``p`` and specific enthalpy ``h``."""

    @abstractmethod
    def isentropic(self, state: FluidState, p: float) -> FluidState:
        """The state reached from ``state`` by an isentropic process to ``p``."""

    @abstractmethod
    def entropy(self, state: FluidState) -> float:
        """Specific entropy of ``state`` [J kg⁻¹ K⁻¹], on this EOS's own datum."""

    @abstractmethod
    def state_hs(self, h: float, s: float) -> FluidState:
        """The state with specific enthalpy ``h`` and specific entropy ``s``."""

    @abstractmethod
    def speed_of_sound(self, state: FluidState) -> float:
        """Speed of sound in ``state`` [m s⁻¹]."""


class Incompressible(EquationOfState):
    """A fluid of constant density.

    Enthalpy is ``cp*T + p/rho``. With constant density, an isentropic
    process is also isothermal, so it moves only the ``p/rho`` term. The
    isentropic head between two pressures is then ``Δp/rho`` exactly: the
    familiar pump head.

    ``cp`` is optional. Without it temperature is carried through isentropic
    processes but cannot be recovered from an enthalpy, so :meth:`state_ph`
    reports it as NaN, and there is no entropy: :meth:`entropy` and
    :meth:`state_hs` need ``cp``. Entropy is then ``cp*ln(T)``, datum 1 K.
    Sound is infinitely fast in an incompressible fluid.

    Parameters
    ----------
    rho : float
        Density [kg m⁻³].
    cp : float, optional
        Specific heat [J kg⁻¹ K⁻¹].
    """

    def __init__(self, rho: float, cp: float | None = None):
        if not (math.isfinite(rho) and rho > 0):
            raise ValueError("rho must be a positive finite number")
        if cp is not None and not (math.isfinite(cp) and cp > 0):
            raise ValueError("cp must be a positive finite number")
        self.rho = float(rho)
        self.cp = None if cp is None else float(cp)

    def __repr__(self):
        return f"Incompressible(rho={self.rho:g}" + ("" if self.cp is None else f", cp={self.cp:g}") + ")"

    def _state(self, p, T, h):
        return FluidState(T=T, p=float(p), h=h, quality=math.nan, rho=self.rho,
                          phase=PHASE_UNKNOWN, T_sat=math.nan)

    def state_pT(self, p, T):
        thermal = 0.0 if self.cp is None else self.cp*float(T)
        return self._state(p, float(T), thermal + float(p)/self.rho)

    def state_ph(self, p, h):
        T = math.nan if self.cp is None else (float(h) - float(p)/self.rho)/self.cp
        return self._state(p, T, float(h))

    def isentropic(self, state, p):
        return self._state(p, state.T, state.h + (float(p) - state.p)/self.rho)

    def _need_cp(self, what):
        if self.cp is None:
            raise ValueError(f"Incompressible needs cp for {what}")

    def entropy(self, state):
        self._need_cp("entropy")
        return self.cp*math.log(state.T)

    def state_hs(self, h, s):
        self._need_cp("a state from enthalpy and entropy")
        T = math.exp(float(s)/self.cp)
        return self._state(self.rho*(float(h) - self.cp*T), T, float(h))

    def speed_of_sound(self, state):
        return math.inf


class CoolPropEOS(EquationOfState):
    """A real fluid, from CoolProp.

    Parameters
    ----------
    fluid : str or Fluid or object with ``.fluid``
        A CoolProp fluid string such as ``'Methane'``, a
        :class:`~pyskyfire.thermo.fluid.Fluid`, or a transport object
        carrying the string on ``.fluid`` (a
        :class:`~pyskyfire.thermo.transport.FluidTransport`).

    Notes
    -----
    A ``(p, T)`` pair lying on the saturation line doesn't identify a state:
    it could be liquid, vapour or anything between. :meth:`state_pT` reads
    such a pair as saturated liquid, a millikelvin into the liquid, because
    a pump inlet drawn from a tank at its own vapour pressure is where this
    comes up. Give the enthalpy through :meth:`state_ph` to mean anything
    else.
    """

    def __init__(self, fluid):
        name = fluid_name(fluid)
        if name is None:
            raise TypeError("fluid must be a CoolProp string, a thermo.Fluid, "
                            "or an object carrying one on .fluid")
        self.fluid = name
        self._resolver = None

    def __repr__(self):
        return f"CoolPropEOS({self.fluid!r})"

    def __getstate__(self):
        # The resolver is a cache of CoolProp lookups; rebuild it after loading.
        return {"fluid": self.fluid, "_resolver": None}

    def resolver(self, p_ref: float) -> FluidStateResolver:
        """The state resolver, probed at the first pressure it is asked about."""
        if self._resolver is None:
            self._resolver = FluidStateResolver(SimpleNamespace(fluid=self.fluid), p_ref)
        return self._resolver

    def state_pT(self, p, T):
        p, T = float(p), float(T)
        try:
            h = CP.PropsSI("Hmass", "T", T, "P", p, self.fluid)
        except ValueError:
            h = CP.PropsSI("Hmass", "T", T - _SATURATION_NUDGE, "P", p, self.fluid)
        return self.state_ph(p, h)

    def state_ph(self, p, h):
        p, h = float(p), float(h)
        resolver = self.resolver(p)
        if resolver.capability.enthalpy_march:
            return resolver.from_hp(h, p)
        # Mixtures and backends the resolver will not flash: slower, and
        # without phase bookkeeping, but still the right state.
        cap = resolver.capability
        return FluidState(
            T=float(CP.PropsSI("T", "Hmass", h, "P", p, self.fluid)), p=p, h=h,
            quality=math.nan, rho=float(CP.PropsSI("Dmass", "Hmass", h, "P", p, self.fluid)),
            phase=(PHASE_SUPERCRITICAL if cap.p_crit is not None and p >= cap.p_crit
                   else PHASE_UNKNOWN),
            T_sat=math.nan)

    def entropy(self, state: FluidState) -> float:
        """Specific entropy of ``state`` [J kg⁻¹ K⁻¹]."""
        return float(CP.PropsSI("Smass", "Hmass", state.h, "P", state.p, self.fluid))

    def isentropic(self, state, p):
        h = CP.PropsSI("Hmass", "Smass", self.entropy(state), "P", float(p), self.fluid)
        return self.state_ph(p, h)

    def state_hs(self, h, s, *, p_guess=1e6, tolerance=1e-10):
        """The state at ``(h, s)``.

        CoolProp's own ``(h, s)`` flash is used where it answers. It refuses
        states beyond the equation's temperature range, though its ``(p, h)``
        flash still extrapolates there, so where it refuses the pressure is
        found along the isenthalp instead. At constant enthalpy
        ``ds/dp = -1/(rho*T)``, and a Newton step in ``ln p`` with that
        slope converges in a few steps from far off.
        """
        h, s = float(h), float(s)
        try:
            return self.state_ph(CP.PropsSI("P", "Hmass", h, "Smass", s, self.fluid), h)
        except ValueError:
            pass
        p = float(p_guess)
        for _ in range(50):
            state = self.state_ph(p, h)
            error = self.entropy(state) - s
            step = error*state.rho*state.T/p
            p *= math.exp(max(min(step, 2.0), -2.0))
            if abs(step) < tolerance:
                return self.state_ph(p, h)
        raise ValueError(f"No state found at h={h:g}, s={s:g} for {self.fluid}")

    def speed_of_sound(self, state):
        """Single-phase speed of sound; NaN inside the dome, where it is not defined."""
        try:
            return float(CP.PropsSI("A", "Hmass", state.h, "P", state.p, self.fluid))
        except ValueError:
            return math.nan


def fluid_name(fluid) -> str | None:
    """The CoolProp name ``fluid`` refers to, or None when it names none.

    Reads a CoolProp string, a :class:`~pyskyfire.thermo.fluid.Fluid`, or
    anything carrying one on ``.fluid``, such as a transport or a
    :class:`CoolPropEOS`. A constant-property stand-in names no fluid.
    """
    fluid = getattr(fluid, "fluid", fluid)
    if isinstance(fluid, Fluid):
        fluid = fluid.coolprop_string()
    return fluid if isinstance(fluid, str) else None


def as_equation_of_state(fluid) -> EquationOfState:
    """Read ``fluid`` as an equation of state.

    An :class:`EquationOfState` is returned as it is. Anything else names a
    CoolProp fluid and is handed to :class:`CoolPropEOS`.
    """
    if isinstance(fluid, EquationOfState):
        return fluid
    return CoolPropEOS(fluid)
