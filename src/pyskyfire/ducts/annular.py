"""Annular ducts about the machine's shaft.

:class:`AnnularDuct` is a passage between a hub and a shroud, both walls of
revolution about the shaft (+x). Its walls are curves in the meridional
``(x, r)`` plane, so it can turn between axial and radial: a turbine's inlet
or exhaust annulus, or the crossover that carries a pump stage's radial
discharge round to the next eye. It is built on
:class:`~pyskyfire.common.MeridionalPassage`, with its span lines drawn normal
to the mean line so the area they sweep is the flow area.

An annulus about a centreline that bends anywhere in space, such as a pipe
round a centrebody or coaxial lines, is a :class:`~pyskyfire.ducts.Pipe` with
a ``core_diameter``.

The duct reduces to station data -- flow area, hydraulic diameter, mean
radius, length and curvature -- and is solved by the march every duct shares
(see :mod:`pyskyfire.ducts.pipe`). It is adiabatic and carries the stagnation
pressure, less the irreversible losses, and the swirl ``c_u`` about the shaft
as the angular momentum ``r*c_u``: unchanged by a change of radius, which is
reversible, and reduced only by wall friction. A meridional turn presses the
flow from one wall to the other equally all the way round, so there is no
secondary flow round the circumference; it is charged as friction on a curved
duct. Every loss coefficient is referred to the through-flow head
``rho*c_m**2/2``.
"""

from __future__ import annotations

from functools import cached_property
import math
from typing import Callable
import warnings

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from pyskyfire.common.friction import annulus_laminar_factor
from pyskyfire.common.meridional import MeridionalPassage
from pyskyfire.thermo.state import FluidState

from .path import _GEOMETRY_TOL
from .pipe import (DuctResult, PipeCorrelationWarning, _area_change_coefficients, _march,
                   _roughness)

__all__ = ["AnnularDuct"]

#: Meridional turns tighter than this many hydraulic diameters get a warning.
#: Roughly where a pipe bend's inner-wall flow separates, and well past the
#: gentle curvature a helical-coil friction multiplier describes.
SHARP_TURN_RATIO = 1.0


def _profile(value, s, name):
    """A per-station array from a scalar, a sequence, or a callable of ``s``."""
    if callable(value):
        values = [value(float(x)) for x in s]
    else:
        values = value
    try:
        array = np.broadcast_to(np.asarray(values, float), (len(s),)).copy()
    except ValueError:
        raise ValueError(f"{name} must be a scalar, one value per station, or a callable of s") from None
    if not np.isfinite(array).all() or np.any(array <= 0.0):
        raise ValueError(f"{name} must be finite and positive at every station")
    array.setflags(write=False)
    return array


def _arc_length(points):
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]


def _unit_tangents(points):
    """Station tangents bisecting the neighbouring chords, as :class:`Path` does.

    Unlike a finite-difference derivative, a bisector never points outside the
    chords either side, so a straight run joining a turn stays straight.
    """
    chord = np.diff(points, axis=0)
    chord /= np.linalg.norm(chord, axis=1)[:, None]
    t = np.vstack((chord[0], chord[:-1] + chord[1:], chord[-1]))
    norm = np.linalg.norm(t, axis=1)
    if np.any(norm < 1.0e-10):
        raise ValueError("the mean line reverses at a station")
    return t/norm[:, None]


def _normals(tangents):
    """Unit normals turned a quarter turn from the tangent, toward the shroud."""
    return np.column_stack((-tangents[:, 1], tangents[:, 0]))


class _Wall:
    """A wall curve as PCHIP of its normalised arc length, with a dense sampling to search."""

    def __init__(self, points, count=2001):
        s = _arc_length(points)
        self.curve = PchipInterpolator(s/s[-1], points)
        self.u = np.linspace(0.0, 1.0, count)
        self.points = self.curve(self.u)

    def crossing(self, origin, tangent):
        """Where the wall crosses the line through ``origin`` normal to ``tangent``, or None."""
        g = (self.points - origin) @ tangent
        hits = np.flatnonzero(np.sign(g[:-1])*np.sign(g[1:]) <= 0)
        if not len(hits):
            return None
        points = []
        for k in hits:
            a, b = self.u[k], self.u[k + 1]
            if g[k] == 0.0 or g[k + 1] == 0.0:
                u = a if g[k] == 0.0 else b
            else:
                u = brentq(lambda v: float((self.curve(v) - origin) @ tangent), a, b, xtol=1e-14)
            points.append(self.curve(u))
        points = np.array(points)
        return points[np.argmin(np.linalg.norm(points - origin, axis=1))]


class AnnularDuct:
    """A stationary annulus about the shaft (+x), turning in the meridional plane.

    The duct is its mean line, an ``(x, r)`` curve from inlet to outlet, and
    the height of the passage across it. Each span line is drawn through a
    mean-line station normal to the mean line, so the frustum it sweeps round
    the axis is the area the through-flow crosses. The walls are the span
    lines' ends. Use :meth:`from_walls` to start from given walls instead.

    Parameters
    ----------
    mean_line : array_like, shape (N, 2)
        Mean-line stations, ``(x, r)`` [m]. They are also the march stations,
        so give enough of them to resolve every turn.
    height : float, sequence or callable
        Span-line length [m]: one value, one per station, or a function of
        mean-line arc length.
    roughness : float or callable, optional
        Absolute wall roughness [m], or a function of arc length. Default 0.

    Notes
    -----
    The passage may not run backward along the shaft: the walls are named
    hub and shroud by radius, and a mean line heading toward -x would swap
    them. Its walls may not fold either, so every turn's radius must exceed
    half the height.
    """

    def __init__(self, mean_line, height, *, roughness: float | Callable[[float], float] = 0.0):
        mean_line = np.array(mean_line, dtype=float)
        if mean_line.ndim != 2 or mean_line.shape[1] != 2 or len(mean_line) < 2:
            raise ValueError("mean_line must have shape (N, 2) with N >= 2")
        if not np.isfinite(mean_line).all():
            raise ValueError("mean_line must be finite")
        if np.any(np.linalg.norm(np.diff(mean_line, axis=0), axis=1) <= _GEOMETRY_TOL):
            raise ValueError("mean_line needs distinct consecutive stations")
        height = _profile(height, _arc_length(mean_line), "height")
        tangents = _unit_tangents(mean_line)
        normal = _normals(tangents)
        self._setup(mean_line - 0.5*height[:, None]*normal,
                    mean_line + 0.5*height[:, None]*normal, roughness)

    @classmethod
    def from_walls(cls, hub, shroud, *, stations: int = 41, roughness: float | Callable[[float], float] = 0.0,
                   max_iterations: int = 100):
        """Re-pair two given walls with span lines normal to their mean line.

        The mean line is found by iteration. Starting from the midpoints of
        the walls sampled at equal fractions of their own lengths, each
        interior station is moved to the midpoint of the chord that the
        normal to the current mean line cuts between the walls, and the
        stations are respaced evenly along it, until they stop moving. The
        end span lines are the walls' own ends, left where they are given:
        they are where the duct joins its neighbours.

        ``hub`` and ``shroud`` are ``(x, r)`` curves such as a
        :class:`~pyskyfire.common.MeridionalPassage`'s; any pairing between
        their points is discarded.
        """
        walls = []
        for name, wall in (("hub", hub), ("shroud", shroud)):
            wall = np.array(wall, dtype=float)
            if wall.ndim != 2 or wall.shape[1] != 2 or len(wall) < 2 or not np.isfinite(wall).all():
                raise ValueError(f"{name} must be finite with shape (N, 2), N >= 2")
            if np.any(np.linalg.norm(np.diff(wall, axis=0), axis=1) <= _GEOMETRY_TOL):
                raise ValueError(f"{name} needs distinct consecutive points")
            walls.append(_Wall(wall))
        if isinstance(stations, bool) or not isinstance(stations, int) or stations < 2:
            raise ValueError("stations must be an integer of at least 2")
        hub_wall, shroud_wall = walls
        u = np.linspace(0.0, 1.0, stations)
        inner, outer = hub_wall.curve(u), shroud_wall.curve(u)
        scale = float(np.ptp(np.vstack((hub_wall.points, shroud_wall.points)), axis=0).max())
        for _ in range(max_iterations):
            mean = 0.5*(inner + outer)
            s = _arc_length(mean)
            mean = np.column_stack([np.interp(u*s[-1], s, mean[:, k]) for k in (0, 1)])
            tangents = _unit_tangents(mean)
            new_inner, new_outer = inner.copy(), outer.copy()
            for j in range(1, stations - 1):
                a = hub_wall.crossing(mean[j], tangents[j])
                b = shroud_wall.crossing(mean[j], tangents[j])
                if a is None or b is None:
                    raise ValueError(f"the span line normal to the mean line at station {j} "
                                     f"does not reach both walls")
                new_inner[j], new_outer[j] = a, b
            shift = max(np.abs(new_inner - inner).max(), np.abs(new_outer - outer).max())
            inner, outer = new_inner, new_outer
            if shift < 1.0e-10*scale:
                break
        else:
            raise ValueError("the mean line between the walls did not converge")
        duct = cls.__new__(cls)
        duct._setup(inner, outer, roughness)
        return duct

    def _setup(self, hub, shroud, roughness):
        mean = 0.5*(hub + shroud)
        tangents = _unit_tangents(mean)
        # A tolerance for discretisation noise at a junction; the passage
        # itself insists on the shroud's radius not falling below the hub's.
        if np.any(tangents[:, 0] < -1.0e-3):
            raise ValueError("an annular duct must not run backward along the shaft (-x)")
        # A radial stretch has its walls at equal radii, where the same noise
        # can leave the shroud a hair inside the hub; level such span lines.
        hub, shroud = hub.copy(), shroud.copy()
        dip = shroud[:, 1] - hub[:, 1]
        level = (dip < 0.0) & (dip > -1.0e-3*np.linalg.norm(shroud - hub, axis=1))
        hub[level, 1] = shroud[level, 1] = 0.5*(hub[level, 1] + shroud[level, 1])
        passage = MeridionalPassage(hub, shroud)
        s = _arc_length(mean)
        angle = np.unwrap(np.arctan2(tangents[:, 1], tangents[:, 0]))
        turning = np.abs(np.diff(angle))
        with np.errstate(divide="ignore"):
            radius = np.where(turning > 1.0e-12, np.diff(s)/np.maximum(turning, 1e-300), np.inf)
            station_radius = 1.0/np.abs(np.gradient(angle, s)) if len(s) > 2 else np.full(len(s), np.inf)
        if np.any(station_radius <= 0.5*passage.height()):
            raise ValueError("a turn's radius must exceed half the passage height, or its inner wall folds")
        self.passage = passage
        self.mean_line = mean
        self.mean_line.setflags(write=False)
        self.s = s
        self.s.setflags(write=False)
        self.curvature_radius = radius
        self.curvature_radius.setflags(write=False)
        self._tangents = tangents
        self.roughness = _roughness(roughness)

    # -- geometry ------------------------------------------------------------
    @property
    def hub(self) -> np.ndarray:
        """Inner wall, ``(x, r)`` at each station [m]."""
        return self.passage.hub

    @property
    def shroud(self) -> np.ndarray:
        """Outer wall, ``(x, r)`` at each station [m]."""
        return self.passage.shroud

    @property
    def height(self) -> np.ndarray:
        """Span-line length at each station [m]."""
        return self.passage.height()

    @property
    def area(self) -> np.ndarray:
        """Flow area at each station [m²]."""
        return self.passage.area()

    @property
    def hydraulic_diameter(self) -> np.ndarray:
        """``4A/P`` at each station, twice the height [m]."""
        return self.passage.hydraulic_diameter()

    @property
    def r_mean(self) -> np.ndarray:
        """Mean-line radius at each station [m]."""
        return self.mean_line[:, 1]

    @property
    def lengths(self) -> np.ndarray:
        """Mean-line length of each interval [m]."""
        return np.diff(self.s)

    @property
    def cone_angle(self) -> np.ndarray:
        """Full cone angle of each interval [rad], positive where the area grows.

        The angle of the cone whose area grows as the duct's does,
        ``2*atan(d(sqrt(A/pi))/ds)``, as :attr:`Pipe.cone_angle
        <pyskyfire.ducts.Pipe.cone_angle>` measures it. An annulus is only
        roughly comparable with a cone, and a radial one carrying swirl
        diffuses without its through-flow alone setting the limit.
        """
        return 2.0*np.arctan(np.diff(np.sqrt(self.area/np.pi))/self.lengths)

    @cached_property
    def area_change_loss_coefficient(self) -> np.ndarray:
        """Loss coefficient of each interval for its change of area [-].

        As :attr:`Pipe.area_change_loss_coefficient
        <pyskyfire.ducts.Pipe.area_change_loss_coefficient>`: Rennels'
        conical-diffuser or conical-contraction term for each run over which
        the area keeps growing or shrinking, read as the cone of the same
        area, on the through-flow head. An annular diffuser is only roughly a
        cone, since how soon it separates depends on both walls' angles, and
        a radial passage is closer to a vaneless diffuser than to either.
        """
        D = self.hydraulic_diameter
        K = _area_change_coefficients(self.area, self._mid_area(), self.lengths,
                                      0.5*(D[:-1] + D[1:]))
        K.setflags(write=False)
        return K

    def end(self, which):
        """The ``'inlet'`` or ``'outlet'`` end, as an :class:`~pyskyfire.common.AnnulusEnd`."""
        return self.passage.end(which)

    def surface(self, wall, n_theta=160):
        """Cartesian points of ``wall`` revolved about +x, shape ``(N, n_theta, 3)``."""
        return self.passage.surface(wall, n_theta)

    def _laminar_factors(self):
        # A wall of revolution curves round the circumference with radius
        # r/cos(beta), beta its angle to the axis. That is the radius of the
        # concentric annulus the section resembles: r_hub for an axial duct,
        # infinite -- parallel plates -- for a radial one.
        t = 0.5*(self._tangents[:-1] + self._tangents[1:])
        axial = np.abs(t[:, 0])/np.linalg.norm(t, axis=1)
        r_hub = 0.5*(self.hub[:-1, 1] + self.hub[1:, 1])
        h = 0.5*(self.height[:-1] + self.height[1:])
        factors = []
        for a, r, b in zip(axial, r_hub, h):
            ratio = 1.0 if a < 1.0e-12 else (r/a)/(r/a + b)
            factors.append(annulus_laminar_factor(min(ratio, 1.0)))
        return np.array(factors)

    def _mid_area(self):
        hub, shroud = 0.5*(self.hub[:-1] + self.hub[1:]), 0.5*(self.shroud[:-1] + self.shroud[1:])
        return np.pi*(hub[:, 1] + shroud[:, 1])*np.linalg.norm(shroud - hub, axis=1)

    # -- solver --------------------------------------------------------------
    def solve(self, *, inlet_state: FluidState, mass_flow: float, transport, inlet_swirl: float = 0.0,
              friction_model: str = "fluids") -> DuctResult:
        """March irreversible losses from inlet to outlet.

        Each interval charges wall friction on the hydraulic diameter, times
        the curvature multiplier where the mean line turns, along the helix a
        swirling flow follows, plus the loss of any change of area
        (:attr:`area_change_loss_coefficient`) on the through-flow head.

        Parameters
        ----------
        inlet_state : FluidState
            Inlet stagnation pressure and temperature.
        mass_flow : float
            [kg/s]
        transport
            Provides ``rho``, ``mu`` and ``cp`` as ``get_<name>(T, p)`` or
            ``<name>(T, p)``.
        inlet_swirl : float, optional
            Circumferential velocity at the inlet mean radius [m/s], positive
            in either sense as long as it is used consistently. Default 0.
        friction_model : {'fluids', 'native'}, optional
            Backend for the friction factor and the curvature multiplier.

        Warns
        -----
        PipeCorrelationWarning
            Where the mean line turns tighter than :data:`SHARP_TURN_RATIO`
            hydraulic diameters, since the distributed curvature model then
            misses the separation a sharp turn causes.

        Notes
        -----
        The area-change loss reads the annulus as a cone of the same area.
        Whether a diffusing annulus reaches its recovery without separating
        is a question for an annular diffuser map (Sovran & Klomp, Japikse),
        which is not applied.
        """
        D = self.hydraulic_diameter
        mid_D = 0.5*(D[:-1] + D[1:])
        if np.any(self.curvature_radius < SHARP_TURN_RATIO*mid_D):
            warnings.warn(
                "annular duct turns tighter than one hydraulic diameter; distributed curvature "
                "friction does not capture separation in a sharp turn and underestimates its loss",
                PipeCorrelationWarning, stacklevel=2,
            )
        count = len(self.s) - 1
        return _march(
            lengths=self.lengths, area=self.area, mid_area=self._mid_area(), hydraulic_diameter=mid_D,
            r_mean=self.r_mean, mid_r_mean=0.5*(self.r_mean[:-1] + self.r_mean[1:]),
            laminar_factor=self._laminar_factors(), curvature_radius=self.curvature_radius,
            bends=[None]*count, bend_diameter=mid_D,
            area_change=self.area_change_loss_coefficient, roughness=self.roughness,
            inlet_state=inlet_state, mass_flow=mass_flow, inlet_swirl=inlet_swirl,
            transport=transport, friction_model=friction_model,
        )
