"""Rules that set how a volute cross-section grows with azimuth.

An area law answers one question: at polar angle ``theta`` measured from the
tongue, how big should the passage be? It answers it as a residual rather than
as an area, because the shaped laws constrain the area *and* the centroid
radius together, and both follow from the section size. The volute solver
drives that residual to zero.

Two laws are provided.

:class:`ConstantAngularMomentum`
    The shaped volute. Aungier calls this the method of choice, since departing
    from conserved angular momentum is what drives the tangential-velocity
    loss.

:class:`ConstantArea`
    The unshaped collector: a passage of uniform cross-section all the way
    around. This is the plain toroidal manifold used on cooling-circuit inlets
    where manufacturing cost outweighs hydraulic efficiency.

References
----------
Aungier, R. H., *Centrifugal Compressors: A Strategy for Aerodynamic Design and
Analysis*, ASME Press, 2000, Chapter 11.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

__all__ = [
    "VoluteInlet",
    "AreaLaw",
    "ConstantAngularMomentum",
    "ConstantArea",
]


@dataclass(frozen=True)
class VoluteInlet:
    """Flow conditions and placement at the volute inlet, Aungier's station 5.

    Attributes
    ----------
    radius : float
        Inlet radius ``r5`` [m], where the collected flow enters the volute.
    C_m : float
        Meridional (through-the-slot) velocity component [m/s]. For a volute
        around a pump discharge this is radial; the tilt of the inlet plane
        rotates it rigidly and does not change its magnitude.
    C_u : float
        Tangential velocity component [m/s]. This is the swirl the volute
        collects, and it must be positive: a volute with no incoming swirl has
        no angular momentum to conserve, and neither the design law nor the
        loss model applies to it.
    """

    radius: float
    C_m: float
    C_u: float

    def __post_init__(self):
        for name in ("radius", "C_m", "C_u"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"VoluteInlet.{name} must be finite and positive")
            object.__setattr__(self, name, value)

    @property
    def C(self) -> float:
        """Inlet velocity magnitude [m/s]. Losses are referred to its head."""
        return float(np.hypot(self.C_m, self.C_u))

    @property
    def tan_alpha(self) -> float:
        """Tangent of the inlet flow angle measured from tangential."""
        return self.C_m / self.C_u

    @property
    def alpha(self) -> float:
        """Inlet flow angle from tangential [rad]."""
        return float(np.arctan2(self.C_m, self.C_u))

    @property
    def swirl(self) -> float:
        """Inlet angular momentum ``r5 * C_u5`` [m^2/s]."""
        return self.radius * self.C_u


class AreaLaw(ABC):
    """A rule relating azimuth to passage size."""

    #: Whether the passage diffuses the flow toward rest as collection begins,
    #: which a uniform collector does and a shaped volute does not. Sets
    #: whether Aungier's Eq. 5-108 correction applies to the tangential loss.
    collector = False

    @abstractmethod
    def residual(self, theta, area, centroid_radius, inlet_width, inlet):
        """Return the signed miss at one azimuth.

        Negative means the trial section is too small. Must increase with
        section size so that the solver can bracket a root.

        Parameters
        ----------
        theta : float
            Polar angle from the tongue [rad], increasing in the direction of
            collection.
        area : float
            Trial section area [m^2].
        centroid_radius : float
            Radius of the trial section's area centroid [m].
        inlet_width : float
            Inlet slot width ``b5`` [m].
        inlet : VoluteInlet
            Inlet flow conditions.
        """


@dataclass(frozen=True)
class ConstantAngularMomentum(AreaLaw):
    """Size the passage to conserve the angular momentum it collects.

    Imposing ``(r*C_u)_c = (r*C_u)_5`` together with mass conservation gives
    Aungier's Eqs. 11-1 to 11-3, generalised by the sizing parameter in
    Eq. 11-5::

        A = theta * SP * r_c * b5 * tan(alpha5)

    so the area divided by the centroid radius grows linearly with azimuth.
    Because both ``A`` and ``r_c`` follow from the section size, this is
    implicit in the size and is solved rather than evaluated.

    Parameters
    ----------
    SP : float, optional
        Sizing parameter: the ratio of actual passage area to the area that
        exactly conserves angular momentum. Defaults to ``1.0``, which is the
        minimum-loss design -- Aungier's tangential-velocity loss (Eq. 5-107)
        vanishes on both of its branches there. Values above one oversize the
        volute deliberately, which costs a little loss at the design point in
        exchange for tolerance to off-design operation; Aungier notes 1.0 to
        1.2 as common practice, biased upward because undersizing is penalised
        far more steeply than oversizing.
    """

    SP: float = 1.0

    def __post_init__(self):
        SP = float(self.SP)
        if not np.isfinite(SP) or SP <= 0.0:
            raise ValueError("SP must be finite and positive")
        object.__setattr__(self, "SP", SP)

    def residual(self, theta, area, centroid_radius, inlet_width, inlet):
        target = theta * self.SP * centroid_radius * inlet_width * inlet.tan_alpha
        return area - target


@dataclass(frozen=True)
class ConstantArea(AreaLaw):
    """Hold the passage area uniform all the way around.

    This is Aungier's "collector" rather than a volute. It collects the flow
    without conserving angular momentum, so its sizing parameter is whatever
    the chosen area happens to produce and its tangential-velocity loss is
    correspondingly larger. Aungier's Eq. 5-108 additionally averages that loss
    against a total loss of tangential head, because a uniform passage tries to
    diffuse the entering flow to rest where collection begins.

    Parameters
    ----------
    area : float
        Passage cross-sectional area [m^2], constant in azimuth.
    """

    area: float
    collector = True

    def __post_init__(self):
        area = float(self.area)
        if not np.isfinite(area) or area <= 0.0:
            raise ValueError("area must be finite and positive")
        object.__setattr__(self, "area", area)

    def residual(self, theta, area, centroid_radius, inlet_width, inlet):
        return area - self.area
