"""What a turbomachine is asked to do, before anything is decided about how.

A duty is the fluid's side of the contract: the mass flow, the state it
arrives in, the pressure it must leave at and the shaft speed. Everything a
design sizes itself on follows from those four, and :class:`MachineDuty` holds
the results in the form designs read:

- the isentropic head, i.e. the enthalpy change of an ideal machine across
  the same pressures (J kg⁻¹, positive for pumps and turbines alike);
- the volumetric flow at the inlet, and at the outlet an ideal machine would
  reach.

Only the isentropic outlet is known here. The real outlet depends on the
machine's losses, and a duty says nothing about the machine.

A duty given as head and flow directly, as a pump often is, has no absolute
state: ``inlet``, ``outlet`` and ``eos`` are then ``None``, and ``mdot`` is
only known if a density was. The pump and turbine packages subclass this with
the constructors and derived numbers each needs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .eos import EquationOfState, as_equation_of_state
from .state import FluidState

__all__ = ["MachineDuty"]


def _positive(value, label):
    if value is None or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be a positive finite number, got {value!r}")
    return float(value)


@dataclass(frozen=True)
class MachineDuty:
    """Head, flow and speed asked of one turbomachine.

    Build one through a subclass's constructors rather than field by field;
    they are what keep the fields consistent with one another.

    Attributes
    ----------
    n : float
        Shaft speed [rpm].
    isentropic_head : float
        Magnitude of the isentropic enthalpy change between the inlet state
        and the outlet pressure [J kg⁻¹].
    Q_in : float
        Volumetric flow at the inlet [m³ s⁻¹].
    Q_out : float
        Volumetric flow at the isentropic outlet state [m³ s⁻¹].
    mdot : float or None
        Mass flow [kg s⁻¹].
    inlet : FluidState or None
        Inlet state.
    outlet : FluidState or None
        Isentropic outlet state: the inlet entropy at the outlet pressure.
    eos : EquationOfState or None
        The equation of state the states were resolved with.
    """

    n: float
    isentropic_head: float
    Q_in: float
    Q_out: float
    mdot: float | None = None
    inlet: FluidState | None = None
    outlet: FluidState | None = None
    eos: EquationOfState | None = None

    def __post_init__(self):
        for label in ("n", "isentropic_head", "Q_in", "Q_out"):
            object.__setattr__(self, label, _positive(getattr(self, label), label))
        if self.mdot is not None:
            object.__setattr__(self, "mdot", _positive(self.mdot, "mdot"))
        if (self.inlet is None) != (self.outlet is None):
            raise ValueError("inlet and outlet are given together or not at all")
        if self.inlet is not None:
            if self.mdot is None:
                raise ValueError("a duty with states needs its mass flow")
            for label, state, Q in (("Q_in", self.inlet, self.Q_in),
                                    ("Q_out", self.outlet, self.Q_out)):
                if not math.isclose(self.mdot/state.rho, Q, rel_tol=1e-9):
                    raise ValueError(f"{label} does not match mdot over its state's density")

    @staticmethod
    def _resolve(fluid, mdot, p_in, T_in, p_out):
        """Inlet and isentropic outlet states, and the fields they give."""
        eos = as_equation_of_state(fluid)
        mdot = _positive(mdot, "mdot")
        _positive(p_in, "p_in")
        _positive(p_out, "p_out")
        inlet = eos.state_pT(p_in, T_in)
        outlet = eos.isentropic(inlet, p_out)
        return dict(isentropic_head=abs(outlet.h - inlet.h), Q_in=mdot/inlet.rho,
                    Q_out=mdot/outlet.rho, mdot=mdot, inlet=inlet, outlet=outlet, eos=eos)

    @property
    def omega(self) -> float:
        """Shaft speed [rad s⁻¹]."""
        return self.n*2*math.pi/60

    @property
    def rho_in(self) -> float | None:
        """Inlet density [kg m⁻³], if the mass flow is known."""
        return None if self.mdot is None else self.mdot/self.Q_in

    @property
    def rho_out(self) -> float | None:
        """Isentropic outlet density [kg m⁻³], if the mass flow is known."""
        return None if self.mdot is None else self.mdot/self.Q_out

    @property
    def isentropic_power(self) -> float | None:
        """Power an ideal machine would exchange with the fluid [W]."""
        return None if self.mdot is None else self.mdot*self.isentropic_head
