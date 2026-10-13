"""The wetted flow path of a turbomachine or annular duct, seen in meridional section.

A passage is two walls of revolution about the shaft (+x), each given as an
``(x, r)`` curve: the ``hub`` nearer the axis and the ``shroud`` outside it.
Both curves have the same number of points, and point ``i`` on one pairs with
point ``i`` on the other. Each pair is the end of a span line, and positions
along the passage are fractions of the hub's arc length. That pairing is what
lets blade sections be placed at a span fraction, and what an annular duct
reads its area from.

The walls also record what an efficiency evaluator needs and the curves
alone do not say:

- :class:`WallMotion`: which stretches of each wall turn with the shaft. A
  rotor's hub is its disk; the casing over its tips stands still.
- :class:`Clearance`: where a blade edge stops short of a wall, and by how
  much. The walls stay where they are drawn: a clearance is a gap between
  the blade and the wall, not a moved wall.

Nothing here knows what kind of machine the passage belongs to. Checks that
only hold for one kind, such as a radial impeller's monotone walls, belong
to that machine.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.interpolate import PchipInterpolator

from .blade_surface import frozen_array

__all__ = ["WallMotion", "Clearance", "AnnulusEnd", "MeridionalPassage", "WALLS"]

#: The two walls of a passage.
WALLS = ("hub", "shroud")


def _fraction_span(start, end, label):
    start, end = float(start), float(end)
    if not (0.0 <= start < end <= 1.0):
        raise ValueError(f"{label} must satisfy 0 <= start < end <= 1, got ({start:g}, {end:g})")
    return start, end


def _wall(name):
    if name not in WALLS:
        raise ValueError(f"wall must be one of {WALLS}, got {name!r}")
    return name


@dataclass(frozen=True)
class WallMotion:
    """Which stretches of each wall rotate with the shaft.

    Each wall lists the ``(start, end)`` station fractions over which it
    turns; everywhere else it stands still. The default is a passage that
    stands still throughout, as a duct does.

    Parameters
    ----------
    hub, shroud : sequence of (start, end)
        Rotating intervals, as fractions of the passage's stations.
    """

    hub: tuple = ()
    shroud: tuple = ()

    def __post_init__(self):
        for name in WALLS:
            intervals = sorted(_fraction_span(a, b, f"{name} rotating interval")
                               for a, b in getattr(self, name))
            for (_, end), (start, _) in zip(intervals, intervals[1:]):
                if start < end:
                    raise ValueError(f"{name} rotating intervals overlap")
            object.__setattr__(self, name, tuple(intervals))

    @classmethod
    def whole(cls, *walls):
        """The named walls rotate over their whole length, the rest stand still."""
        return cls(**{_wall(name): ((0.0, 1.0),) for name in walls})

    def rotating(self, wall, fraction):
        """Whether ``wall`` turns at each station ``fraction``. Ends count as inside."""
        f = np.asarray(fraction, dtype=float)
        result = np.zeros(f.shape, dtype=bool)
        for start, end in getattr(self, _wall(wall)):
            result |= (f >= start) & (f <= end)
        return result


@dataclass(frozen=True)
class Clearance:
    """A blade edge that stops short of a passage wall.

    Parameters
    ----------
    blade : str
        Which blade this is, in the owning machine's own naming, such as
        ``'main'`` or ``'R1'``.
    wall : {'hub', 'shroud'}
        The wall the edge faces.
    start, end : float
        Station fractions the edge runs between.
    gap : float
        Distance from the edge to the wall, along the span line [m].
    """

    blade: str
    wall: str
    start: float
    end: float
    gap: float

    def __post_init__(self):
        _wall(self.wall)
        start, end = _fraction_span(self.start, self.end, f"clearance {self.blade!r}")
        object.__setattr__(self, "start", start)
        object.__setattr__(self, "end", end)
        if not (np.isfinite(self.gap) and self.gap > 0):
            raise ValueError(f"clearance {self.blade!r} needs a positive finite gap")
        object.__setattr__(self, "gap", float(self.gap))


@dataclass(frozen=True)
class AnnulusEnd:
    """One end of a passage: what anything joining it there has to match.

    Attributes
    ----------
    x : float
        Axial position of the hub end of the span line [m].
    hub, shroud : ndarray, shape (2,)
        ``(x, r)`` of each wall's end point.
    hub_tangent, shroud_tangent : ndarray, shape (2,)
        Unit ``(x, r)`` direction of each wall, pointing downstream.
    """

    hub: np.ndarray
    shroud: np.ndarray
    hub_tangent: np.ndarray
    shroud_tangent: np.ndarray

    @property
    def x(self) -> float:
        return float(self.hub[0])

    @property
    def r_hub(self) -> float:
        return float(self.hub[1])

    @property
    def r_shroud(self) -> float:
        return float(self.shroud[1])

    @property
    def r_mean(self) -> float:
        return 0.5*(self.r_hub + self.r_shroud)

    @property
    def height(self) -> float:
        return float(np.linalg.norm(self.shroud - self.hub))

    @property
    def area(self) -> float:
        """Area the span line sweeps, as :meth:`MeridionalPassage.area`."""
        return float(np.pi*(self.r_hub + self.r_shroud)*self.height)


@dataclass(frozen=True)
class MeridionalPassage:
    """Paired hub and shroud curves, in ``(x, r)`` metres, inlet to outlet.

    Parameters
    ----------
    hub, shroud : array_like, shape (N, 2)
        Wall curves. Point ``i`` on each is one end of span line ``i``.
    motion : WallMotion, optional
        Which wall stretches rotate. Default: none.
    clearances : sequence of Clearance, optional
        Blade edges standing off a wall.
    """

    hub: np.ndarray
    shroud: np.ndarray
    # A plain default rather than a factory, so it lives on the class and a
    # passage pickled before the field existed still reads as standing still.
    motion: WallMotion = WallMotion()
    clearances: tuple = ()

    def __post_init__(self):
        for name in WALLS:
            p = frozen_array(getattr(self, name))
            if p.ndim != 2 or p.shape[1] != 2 or len(p) < 2:
                raise ValueError(f"{name} must have shape (N, 2) with N >= 2")
            if np.any(p[:, 1] < 0):
                raise ValueError(f"{name} radii must not be negative")
            if np.any(np.linalg.norm(np.diff(p, axis=0), axis=1) <= 0):
                raise ValueError(f"{name} needs distinct consecutive stations")
            object.__setattr__(self, name, p)
        if self.hub.shape != self.shroud.shape:
            raise ValueError("Hub and shroud stations must be paired")
        if np.any(self.shroud[:, 1] < self.hub[:, 1] - 1e-12):
            raise ValueError("Shroud radius must not be below hub radius")
        if np.any(self.height() <= 1e-12):
            raise ValueError("Hub and shroud must remain separated, including the outlet")
        if not isinstance(self.motion, WallMotion):
            raise TypeError("motion must be a WallMotion")
        clearances = tuple(self.clearances)
        for c in clearances:
            if not isinstance(c, Clearance):
                raise TypeError("clearances must be Clearance objects")
            f = np.linspace(c.start, c.end, 33)
            if np.any(self.height(f) <= c.gap):
                raise ValueError(f"clearance {c.blade!r} is not smaller than the passage height")
        object.__setattr__(self, "clearances", clearances)

    # -- stations ------------------------------------------------------------
    @property
    def stations(self):
        """Station fractions: the hub's arc length, normalised to 0..1."""
        s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(self.hub, axis=0), axis=1))]
        return s/s[-1]

    def _interpolate(self, curve):
        # Pchip needs at least two points and keeps a two-point wall straight.
        return PchipInterpolator(self.stations, curve)

    def sample(self, fraction, span_fraction):
        """``(x, r)`` at station ``fraction`` and ``span_fraction`` (0 hub, 1 shroud)."""
        h = self._interpolate(self.hub)(fraction)
        s = self._interpolate(self.shroud)(fraction)
        return h + np.asarray(span_fraction)[..., None]*(s - h)

    def wall_tangents(self, fraction):
        """Unit hub and shroud tangents at the ends of one span line."""
        tangents = [self._interpolate(wall).derivative()(fraction) for wall in (self.hub, self.shroud)]
        return np.array([t/np.linalg.norm(t) for t in tangents])

    # -- section properties --------------------------------------------------
    def _span_ends(self, fraction):
        if fraction is None:
            return self.hub, self.shroud
        return self.sample(fraction, 0.0), self.sample(fraction, 1.0)

    def height(self, fraction=None):
        """Span-line length [m], at the stored stations or at ``fraction``."""
        hub, shroud = self._span_ends(fraction)
        return np.linalg.norm(shroud - hub, axis=-1)

    def area(self, fraction=None):
        """Area of the surface the span line sweeps round the axis [m²].

        That is the frustum ``pi*(r_hub + r_shroud)*height``: an annulus
        ``pi*(r_s² - r_h²)`` where the span line is radial, a cylinder
        ``2*pi*r*b`` where it is axial. It is the flow area where the flow
        crosses the span line square, which is what a span line is drawn for.
        """
        hub, shroud = self._span_ends(fraction)
        return np.pi*(hub[..., 1] + shroud[..., 1])*np.linalg.norm(shroud - hub, axis=-1)

    def hydraulic_diameter(self, fraction=None):
        """``4A/P`` over both walls [m], which for an annulus is twice its height."""
        return 2*self.height(fraction)

    # -- walls and blades ----------------------------------------------------
    def rotating(self, wall, fraction=None):
        """Whether ``wall`` turns with the shaft, at the stored stations or at ``fraction``."""
        return self.motion.rotating(wall, self.stations if fraction is None else fraction)

    def clearance(self, blade):
        """The clearances recorded for ``blade``, as a tuple."""
        return tuple(c for c in self.clearances if c.blade == blade)

    def with_clearance(self, clearance, wall="shroud"):
        """This passage with ``wall`` moved a distance ``clearance`` towards the other.

        This is the boundary a blade is built to, not the wetted wall: the gap
        between the two is the clearance. The pairing, motion and recorded
        clearances are kept.
        """
        _wall(wall)
        gap = self.height()
        if clearance < 0 or np.any(gap <= clearance):
            raise ValueError("tip_clearance must be nonnegative and smaller than the passage height")
        moved, other = (self.shroud, self.hub) if wall == "shroud" else (self.hub, self.shroud)
        return replace(self, **{wall: moved + clearance/gap[:, None]*(other - moved)})

    def end(self, which):
        """The ``'inlet'`` or ``'outlet'`` end, as an :class:`AnnulusEnd`."""
        if which not in ("inlet", "outlet"):
            raise ValueError(f"which must be 'inlet' or 'outlet', got {which!r}")
        i, f = (0, 0.0) if which == "inlet" else (-1, 1.0)
        hub_t, shroud_t = self.wall_tangents(f)
        return AnnulusEnd(self.hub[i].copy(), self.shroud[i].copy(), hub_t, shroud_t)

    def surface(self, wall, n_theta=160):
        """Cartesian points of ``wall`` revolved about +x, shape ``(N, n_theta, 3)``."""
        curve = getattr(self, _wall(wall))
        theta = np.linspace(0, 2*np.pi, n_theta)
        x, angle = np.meshgrid(curve[:, 0], theta, indexing="ij")
        radius = curve[:, 1, None]
        return np.stack((x, radius*np.cos(angle), radius*np.sin(angle)), axis=-1)
