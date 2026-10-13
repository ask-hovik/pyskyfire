"""The duty a pump is designed for: flow, head and speed.

A pump is usually specified as ``Q``, ``H`` and ``n``, but in an engine those
fall out of a medium, a mass flow and two pressures. :class:`PumpDuty` takes
either and stores the same thing, so a design reads it the same way whichever
the user had to hand.

``H`` is the isentropic head in metres of the fluid pumped: the enthalpy an
ideal pump adds, over ``g``. For a constant density it is ``Δp/(rho*g)``. For
a real fluid it is not quite, since the fluid compresses on the way through,
and the duty says by how much.
"""

from __future__ import annotations

from dataclasses import dataclass

from pyskyfire.thermo.duty import MachineDuty, _positive
from pyskyfire.thermo.eos import Incompressible

from .impeller import G0, specific_speed

__all__ = ["PumpDuty"]


@dataclass(frozen=True)
class PumpDuty(MachineDuty):
    """What a pump is asked to do. See :class:`~pyskyfire.thermo.MachineDuty`.

    Build one with a constructor:

    - :meth:`from_head` for flow, head and speed, as pumps are specified;
    - :meth:`from_pressure_rise` for a mass flow and a pressure rise, at one
      density;
    - :meth:`from_states` for a medium, a mass flow and two pressures;
    - :meth:`from_stations` for two solved network stations.

    Pressures are total pressures, as the network carries them, so the head
    is total-to-total.
    """

    @classmethod
    def from_head(cls, Q, H, n, *, rho=None):
        """Flow [m³ s⁻¹], head [m] and speed [rpm], and a density if known.

        Without a density there is no mass flow and no state, which is all
        the impeller's geometry needs.
        """
        Q = _positive(Q, "Q")
        head = _positive(H, "H")*G0
        if rho is None:
            return cls(n=n, isentropic_head=head, Q_in=Q, Q_out=Q)
        eos = Incompressible(rho)
        return cls(n=n, isentropic_head=head, Q_in=Q, Q_out=Q, mdot=rho*Q, eos=eos)

    @classmethod
    def from_pressure_rise(cls, *, mdot, rho, dp, n):
        """Mass flow [kg s⁻¹], density [kg m⁻³], pressure rise [Pa], speed [rpm]."""
        rho = _positive(rho, "rho")
        return cls.from_head(_positive(mdot, "mdot")/rho, _positive(dp, "dp")/(rho*G0), n, rho=rho)

    @classmethod
    def from_states(cls, *, fluid, mdot, p_in, T_in, p_out, n):
        """A medium, a mass flow and the pressures either side of the pump.

        ``fluid`` is a CoolProp name, a :class:`~pyskyfire.thermo.Fluid`, or
        any :class:`~pyskyfire.thermo.EquationOfState`, including
        :class:`~pyskyfire.thermo.Incompressible`.
        """
        if not p_out > p_in:
            raise ValueError("A pump's outlet pressure must be above its inlet pressure")
        return cls(n=n, **cls._resolve(fluid, mdot, p_in, T_in, p_out))

    @classmethod
    def from_stations(cls, st_in, st_out, *, fluid, n):
        """Two network stations, by their pressures and the inlet's T and mdot.

        The outlet station's temperature is not read: it is what the real
        pump's losses made of the flow, and a duty is about the pressure it
        was asked to reach.
        """
        return cls.from_states(fluid=fluid, mdot=st_in.mdot, p_in=st_in.p, T_in=st_in.T,
                               p_out=st_out.p, n=n)

    @property
    def Q(self) -> float:
        """Volumetric flow at the inlet [m³ s⁻¹], the flow a pump is rated at."""
        return self.Q_in

    @property
    def H(self) -> float:
        """Isentropic head [m of the fluid pumped]."""
        return self.isentropic_head/G0

    @property
    def specific_speed(self) -> float:
        """Gülich's ``n_q``: rpm, m³ s⁻¹, m."""
        return float(specific_speed(self.n, self.Q, self.H))
