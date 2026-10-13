"""Built turbine blades, and what can only be measured on them.

A :class:`TurbineGeometry` is what :meth:`TurbineDefinition.build_geometry`
returns: one :class:`RowGeometry` per row, and the passage they sit in. Each
row carries one blade's surfaces, which the scene repeats round the wheel,
and the numbers the drawing left to the built section:

- the throat: the shortest line from one blade's trailing edge to its
  neighbour, at every span station. With the pitch it gives the gauging
  angle ``acos(o/s)``, which is roughly where the flow leaves the row.
  Setting that beside the angle the design asked for shows how well the
  section does what was meant;
- the narrowest gap anywhere between neighbouring blades. Set it beside the
  width the flow comes in through (``inlet_width``, the pitch times the
  cosine of the inlet metal angle) and the throat it leaves through. A gap
  narrower than both is a pinch: the flow speeds up into it and slows again
  after, which costs loss. It chokes only if the flow there reaches sonic
  speed. An impulse blade takes its flow in as steeply as it lets it out, so
  its narrowest gap is usually at the inlet, and the throat alone is the
  wrong yardstick for it;
- the maximum thickness the section actually has.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyskyfire.common.blade_surface import BladeSurface, frozen_array
from pyskyfire.common.meridional import MeridionalPassage
from pyskyfire.common.rotation import COUNTERCLOCKWISE, parse_rotation, rotation_sign

from .definition import BladeRow

__all__ = ["RowGeometry", "TurbineGeometry"]


@dataclass(frozen=True)
class RowGeometry:
    """One row's blade and its measured section properties.

    Attributes
    ----------
    row : BladeRow
        The row as drawn.
    blade : BladeSurface
        One blade, ``(span, chord, 3)`` in ``(x, r, theta)``, with its section
        centroids at ``theta = 0``.
    span_fractions : ndarray
        Span station of each section, 0 at the hub end of the blade.
    section_radius : ndarray
        Radius each section was drawn at and its pitch taken at [m].
    inlet_metal_angle, outlet_metal_angle : ndarray
        Metal angles of each section [rad].
    throat : ndarray
        Shortest distance from the trailing edge to the neighbouring blade,
        per section [m].
    max_thickness : ndarray
        Largest thickness of each section [m].
    narrowest : ndarray
        Narrowest gap anywhere between neighbouring blades, per section [m].
    narrowest_at : ndarray
        Where that gap is, as a fraction of the axial chord from the leading
        edge, per section.
    """

    row: BladeRow
    blade: BladeSurface
    span_fractions: np.ndarray
    section_radius: np.ndarray
    inlet_metal_angle: np.ndarray
    outlet_metal_angle: np.ndarray
    throat: np.ndarray
    max_thickness: np.ndarray
    narrowest: np.ndarray
    narrowest_at: np.ndarray

    def __post_init__(self):
        for name in ("span_fractions", "section_radius", "inlet_metal_angle", "outlet_metal_angle",
                     "throat", "max_thickness", "narrowest", "narrowest_at"):
            object.__setattr__(self, name, frozen_array(getattr(self, name)))

    @property
    def pitch(self) -> np.ndarray:
        """Blade spacing at each section's radius [m]."""
        return 2*np.pi*self.section_radius/self.row.blade_count

    @property
    def inlet_width(self) -> np.ndarray:
        """Width the flow enters through, ``pitch*cos(inlet metal angle)`` [m]."""
        return self.pitch*np.cos(self.inlet_metal_angle)

    @property
    def pinch(self) -> np.ndarray:
        """Narrowest gap over the smaller of the inlet width and the throat.

        One where the passage narrows no further than its ends; below one, the
        flow accelerates into the pinch and diffuses out of it.
        """
        return self.narrowest/np.minimum(self.inlet_width, self.throat)

    @property
    def narrowest_height(self) -> np.ndarray:
        """Blade height where the narrowest gap is, between the row's edge heights [m].

        With it, ``narrowest*narrowest_height`` is the passage's smallest flow
        area, to set beside ``inlet_width*row.height_in``.
        """
        r = self.row
        return r.height_in + self.narrowest_at*(r.height_out - r.height_in)

    @property
    def gauging_angle(self) -> np.ndarray:
        """``acos(throat/pitch)``, signed as the outlet metal angle [rad]."""
        ratio = np.clip(self.throat/self.pitch, 0.0, 1.0)
        return np.sign(self.outlet_metal_angle)*np.arccos(ratio)

    def at_mid_span(self, values) -> float:
        """``values`` (one per section) interpolated to mid-span."""
        return float(np.interp(0.5, self.span_fractions, values))


@dataclass(frozen=True)
class TurbineGeometry:
    """Every row's built blade, and the passage they stand in."""

    passage: MeridionalPassage
    rows: tuple
    rotation: str = COUNTERCLOCKWISE
    backend_version: str = ""

    def __post_init__(self):
        object.__setattr__(self, "rows", tuple(self.rows))
        object.__setattr__(self, "rotation", parse_rotation(self.rotation))

    @property
    def rotation_sign(self) -> float:
        return rotation_sign(self.rotation)

    def row(self, name) -> RowGeometry:
        for r in self.rows:
            if r.row.name == name:
                return r
        raise KeyError(name)


def _max_thickness(suction, pressure):
    """Largest distance from a suction-side point to the pressure side."""
    from scipy.spatial import cKDTree
    return float(cKDTree(pressure).query(suction)[0].max())


def cascade_gaps(points, trailing_edge, pitch):
    """The throat and the narrowest gap of a cascade, from one section's outline.

    ``points`` outline one section in the unrolled cascade plane, as
    ``(axial, pitchwise)``, and ``trailing_edge`` are the outline points
    around its trailing edge. The neighbours are the same outline one pitch
    either side. The throat is the shortest distance from the trailing edge
    to a neighbour; the narrowest gap is the shortest from anywhere on the
    blade. Returns the throat, the narrowest gap, and the axial position of
    the gap's midpoint.
    """
    from scipy.spatial import cKDTree
    shift = np.array([0.0, pitch])
    neighbours_points = np.vstack((points + shift, points - shift))
    neighbours = cKDTree(neighbours_points)
    gaps, nearest = neighbours.query(points)
    k = int(np.argmin(gaps))
    x = 0.5*(points[k, 0] + neighbours_points[nearest[k], 0])
    return float(neighbours.query(trailing_edge)[0].min()), float(gaps[k]), float(x)
