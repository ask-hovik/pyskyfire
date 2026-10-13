"""The duty a turbine is designed for: a gas expanding between two pressures.

Unlike a pump's, a turbine's duty cannot be given as head and flow alone.
The gas expands on its way through, so the volumetric flow at the exit can be
a good deal larger than at the inlet, and every blade height downstream is
sized on that. :class:`TurbineDuty` therefore always carries the states: a
medium, a mass flow, the inlet state and the pressure the gas leaves at.

The inlet is a total state. The outlet pressure is whatever the caller means
by it, and ``outlet_basis`` records which: ``'total'`` when the kinetic energy
leaving the last row counts as delivered (as network stations carry it), or
``'static'`` when it is charged to the turbine as a loss. The isentropic head
is computed the same way either way; only what a design does with the leaving
velocity differs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from pyskyfire.thermo.duty import MachineDuty, _positive

__all__ = ["TurbineDuty"]

_BASES = ("total", "static")


@dataclass(frozen=True)
class TurbineDuty(MachineDuty):
    """What a turbine is asked to do. See :class:`~pyskyfire.thermo.MachineDuty`.

    Build one with :meth:`from_states` or :meth:`from_stations`.

    Attributes
    ----------
    outlet_basis : {'total', 'static'}
        Whether the outlet pressure is total or static.
    """

    outlet_basis: str = "total"

    def __post_init__(self):
        super().__post_init__()
        if self.outlet_basis not in _BASES:
            raise ValueError(f"outlet_basis must be one of {_BASES}, got {self.outlet_basis!r}")
        if self.inlet is None:
            raise ValueError("A turbine duty needs its states: build it with "
                             "TurbineDuty.from_states or TurbineDuty.from_stations")

    @classmethod
    def from_states(cls, *, fluid, mdot, p0_in, T0_in, p_out, n, outlet_basis="total"):
        """A medium, a mass flow, the inlet total state and the outlet pressure.

        ``fluid`` is a CoolProp name, a :class:`~pyskyfire.thermo.Fluid`, or
        any :class:`~pyskyfire.thermo.EquationOfState`.
        """
        if not p_out < p0_in:
            raise ValueError("A turbine's outlet pressure must be below its inlet pressure")
        return cls(n=n, outlet_basis=outlet_basis, **cls._resolve(fluid, mdot, p0_in, T0_in, p_out))

    @classmethod
    def from_stations(cls, st_in, st_out, *, fluid, n):
        """Two network stations, by their pressures and the inlet's T and mdot.

        Network stations are total states, so the outlet basis is total. The
        outlet station's temperature is not read: it is what the real
        turbine's losses made of the flow.
        """
        return cls.from_states(fluid=fluid, mdot=st_in.mdot, p0_in=st_in.p, T0_in=st_in.T,
                               p_out=st_out.p, n=n, outlet_basis="total")

    @property
    def pressure_ratio(self) -> float:
        """Inlet total pressure over outlet pressure."""
        return self.inlet.p/self.outlet.p

    @property
    def spouting_velocity(self) -> float:
        """``c0 = sqrt(2*Δh_s)`` [m s⁻¹], the speed the whole head would make."""
        return math.sqrt(2*self.isentropic_head)

    @property
    def specific_speed(self) -> float:
        """Dimensionless ``Ns = omega*sqrt(Q_out)/Δh_s**0.75``.

        The flow is taken at the isentropic outlet, the convention of Balje's
        charts and of turboflow.
        """
        return self.omega*math.sqrt(self.Q_out)/self.isentropic_head**0.75

    def specific_diameter(self, diameter) -> float:
        """Dimensionless ``Ds = D*Δh_s**0.25/sqrt(Q_out)`` for a rotor of ``diameter`` [m]."""
        return _positive(diameter, "diameter")*self.isentropic_head**0.25/math.sqrt(self.Q_out)

    def blade_speed(self, diameter) -> float:
        """Blade speed [m s⁻¹] at ``diameter`` [m]."""
        return self.omega*_positive(diameter, "diameter")/2

    def blade_jet_ratio(self, diameter) -> float:
        """``u/c0`` at ``diameter`` [m]."""
        return self.blade_speed(diameter)/self.spouting_velocity
