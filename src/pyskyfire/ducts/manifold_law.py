"""Rules that set how a manifold's passage grows as it gathers flow.

A volute's area law answers to swirl: the passage is sized so the angular
momentum it collects is conserved, and departing from that is what costs the
tangential loss. A chamber manifold has no swirl to conserve. Its circumferential
velocity is a nuisance rather than a design variable -- it exists only because
fluid has to travel round to reach the far channels -- so the laws here answer a
different question: how much dynamic head does the passage carry, and is it the
same everywhere the channels tap into it?

Two are provided.

:class:`UniformArea`
    A passage of constant cross-section. Cheapest to make, and perfectly valid
    where the manifold is oversized enough that its dynamic head is small
    against the channel pressure drop, which is the usual design rule.

:class:`ConstantVelocity`
    Area in proportion to the flow gathered so far, so the through-flow
    velocity is the same at every station. The point is not the velocity but
    the head: a uniform dynamic head leaves every branch seeing the same static
    pressure, which is what makes the split even.

Not provided, and worth knowing about: holding *static* pressure uniform is not
the same as holding velocity uniform. A dividing header decelerates, so
momentum recovery lifts static pressure back and partly cancels friction, and
the taper that truly flattens static pressure follows from a momentum balance
that includes the branch flow. It needs the branch pressure drop, which couples
the manifold to the circuit it feeds, so it waits until that coupling exists.

The laws are stated as geometry rather than as an operating point, the way
:class:`~pyskyfire.ducts.area_law.ConstantArea` is: a duty is turned into an
area once, through :meth:`ConstantVelocity.for_duty`, and the law itself then
carries no fluid state.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

__all__ = ["ManifoldLaw", "UniformArea", "ConstantVelocity"]


class ManifoldLaw(ABC):
    """A rule relating flow gathered to passage size."""

    @abstractmethod
    def residual(self, collected, area):
        """Return the signed miss at one station.

        Negative means the trial section is too small. Must increase with
        section area so that the solver can bracket a root.

        Parameters
        ----------
        collected : float
            Fraction of the manifold's whole mass flow that has passed through
            this station. Runs from zero at the dead point to ``1/branches`` at
            the port, so a law never has to know how the ring is divided.
        area : float
            Trial section area [m^2].
        """


@dataclass(frozen=True)
class UniformArea(ManifoldLaw):
    """Hold the passage area constant all the way round.

    The plain toroidal manifold. It carries more dynamic head near the port
    than at the dead point, since the same area passes a growing share of the
    flow, so it relies on being oversized enough that the variation does not
    reach the channels. That is the usual rule for a coolant manifold, and it
    is why the parameter here is simply the area the designer chose.

    Parameters
    ----------
    area : float
        Passage cross-sectional area [m^2].
    """

    area: float

    def __post_init__(self):
        area = float(self.area)
        if not np.isfinite(area) or area <= 0.0:
            raise ValueError("area must be finite and positive")
        object.__setattr__(self, "area", area)

    def residual(self, collected, area):
        return area - self.area


@dataclass(frozen=True)
class ConstantVelocity(ManifoldLaw):
    """Grow the passage in step with the flow it has gathered.

    ``A = collected * full_area`` holds the through-flow velocity uniform, and
    with it the dynamic head, so every branch tapping the manifold sees the
    same static pressure. That evenness is the point; the velocity is how it is
    obtained.

    The area goes to zero at the dead point, which no section family can build.
    Stations where the law asks for less than the smallest section are flagged
    through ``geometry.clipped``, exactly as a volute's tongue is.

    Parameters
    ----------
    full_area : float
        Area that would carry the manifold's whole mass flow at the design
        velocity [m^2]. A ring split into two branches reaches half of it at
        the port, since each branch only ever gathers half the flow.
    """

    full_area: float

    def __post_init__(self):
        full_area = float(self.full_area)
        if not np.isfinite(full_area) or full_area <= 0.0:
            raise ValueError("full_area must be finite and positive")
        object.__setattr__(self, "full_area", full_area)

    @classmethod
    def for_duty(cls, *, velocity, mass_flow, density):
        """Size the law from a duty and the velocity it should be carried at."""
        for name, value in (("velocity", velocity), ("mass_flow", mass_flow),
                            ("density", density)):
            value = float(value)
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        return cls(float(mass_flow)/(float(density)*float(velocity)))

    def velocity(self, mass_flow, density):
        """The through-flow velocity this law holds, for a given duty [m/s]."""
        return float(mass_flow)/(float(density)*self.full_area)

    def residual(self, collected, area):
        return area - collected*self.full_area
