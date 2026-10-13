"""The drawing of an axial turbine's bladed flow path.

A :class:`TurbineDefinition` is hardware and nothing else: an ordered tuple of
:class:`BladeRow` objects and the annular passage they sit in. It holds no
velocity triangles, no duty and no efficiency, so it does not care how it was
designed. A meanline design produces one, a user can write one out by hand,
and evaluators and the geometry builder read one.

Conventions
-----------
- The shaft is +x and the flow runs along +x.
- Angles are radians from axial, positive in the direction the rotor turns,
  which is the direction of blade speed. The machine's ``rotation`` then says
  which way that is about +x, so a clockwise turbine is the mirror image of a
  counterclockwise one and nothing else changes.
- Metal angles are the direction of the camber line at each edge. The flow
  meets and leaves at slightly different angles: incidence is metal minus
  flow at the inlet, and deviation is how far the leaving flow falls short of
  the outlet metal angle.
- Section properties (radii, metal angles, stagger) are given at mean radius.
  ``twist`` says how the metal angles change from hub to tip; without it the
  section is the same at every span.

What is *not* here is anything the built section decides: the throat, and
the realised maximum thickness. :meth:`TurbineDefinition.build_geometry`
measures those from the blades it builds.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np

from pyskyfire.common.meridional import Clearance, MeridionalPassage, WallMotion
from pyskyfire.common.rotation import COUNTERCLOCKWISE, parse_rotation, rotation_sign

__all__ = ["BladeRow", "FreeVortexTwist", "ThicknessDistribution", "TurbineDefinition",
           "metal_angles"]

KINDS = ("stator", "rotor")
_MAX_ANGLE = math.radians(85)


def _positive(value, label):
    if not (np.isfinite(value) and value > 0):
        raise ValueError(f"{label} must be a positive finite number, got {value!r}")
    return float(value)


def metal_angles(flow_in, flow_out, incidence=0.0, deviation=0.0):
    """Metal angles for a row turning the flow from ``flow_in`` to ``flow_out``.

    Incidence is added to the inlet flow angle. Deviation is added to the
    outlet in the direction the row turns the flow, since the leaving flow
    falls short of the metal and the blade has to turn it further.
    """
    turning = np.sign(flow_out - flow_in) or 1.0
    return flow_in + incidence, flow_out + turning*deviation


@dataclass(frozen=True)
class FreeVortexTwist:
    """Metal angles from hub to tip for a free-vortex row.

    In a free vortex the axial velocity is the same at every radius and the
    swirl falls off as ``1/r``, so ``r*c_u`` is constant across the span. A
    rotor's blade speed grows with ``r``, so its relative angles twist more
    than a stator's absolute ones.

    Parameters
    ----------
    r_mean : float
        Radius the swirl and blade speed are given at [m].
    c_x_in, c_x_out : float
        Axial velocity at the row's inlet and outlet [m s⁻¹].
    c_u_in, c_u_out : float
        Absolute swirl at mean radius, positive in the rotation direction [m s⁻¹].
    u_mean : float
        Blade speed at mean radius: zero for a stator [m s⁻¹].
    incidence, deviation : float
        As :func:`metal_angles` [rad].
    """

    r_mean: float
    c_x_in: float
    c_x_out: float
    c_u_in: float
    c_u_out: float
    u_mean: float = 0.0
    incidence: float = 0.0
    deviation: float = 0.0

    def flow_angles(self, r):
        """Inlet and outlet flow angles, in the row's own frame, at radius ``r``."""
        ratio = self.r_mean/r
        u = self.u_mean/ratio
        return (math.atan((self.c_u_in*ratio - u)/self.c_x_in),
                math.atan((self.c_u_out*ratio - u)/self.c_x_out))

    def __call__(self, span_fraction, r):
        return metal_angles(*self.flow_angles(r), self.incidence, self.deviation)


@dataclass(frozen=True)
class ThicknessDistribution:
    """How a section's thickness is laid along its camber line.

    The section is drawn as offsets from the camber line: control points on
    the suction side and on the pressure side, each a fraction of the row's
    ``max_thickness``, spaced from leading edge towards trailing edge. The
    built section's actual thickness is measured, not assumed from these.

    Parameters
    ----------
    suction, pressure : sequence of float
        Offsets of each side's control points, as fractions of the maximum
        thickness. A negative pressure-side offset bends that side back
        across the camber line.
    suction_extent, pressure_extent : float
        Fraction of the camber line the control points spread over.
    spacing : float
        Growth ratio between control-point spacings: 1 spaces them evenly,
        larger values bunch them towards the leading edge.

    The default is front-loaded, thickest near the leading edge, as a
    reaction blade or a nozzle vane is. :meth:`symmetric` peaks at
    mid-chord and thins towards both edges, as an impulse blade does.
    """

    suction: tuple = (1.0, 0.9, 0.7, 0.5)
    pressure: tuple = (0.4, 0.3, 0.15)
    suction_extent: float = 0.8
    pressure_extent: float = 0.95
    spacing: float = 1.2

    def __post_init__(self):
        for name in ("suction", "pressure"):
            values = tuple(float(v) for v in getattr(self, name))
            if not values or not all(np.isfinite(values)):
                raise ValueError(f"{name} needs at least one finite offset")
            object.__setattr__(self, name, values)
        for name in ("suction_extent", "pressure_extent"):
            if not 0 < getattr(self, name) <= 1:
                raise ValueError(f"{name} must lie in (0, 1]")
        _positive(self.spacing, "spacing")

    @classmethod
    def symmetric(cls):
        """Thickest at mid-chord, thin at both edges: an impulse blade."""
        return cls(suction=(0.5, 0.8, 1.0, 0.8, 0.5), pressure=(0.2, 0.3, 0.2),
                   suction_extent=0.9, pressure_extent=0.9, spacing=1.0)


@dataclass(frozen=True)
class BladeRow:
    """One row of blades, as hardware.

    Parameters
    ----------
    name : str
        The row's name, such as ``'S1'`` or ``'R1'``.
    kind : {'stator', 'rotor'}
        Whether the row turns with the shaft.
    blade_count : int
        Blades around the full circumference.
    x_le : float
        Axial position of the leading edge [m].
    axial_chord : float
        Axial extent of the blade [m].
    r_hub_in, r_tip_in, r_hub_out, r_tip_out : float
        Passage radii at the leading and trailing edges [m].
    inlet_metal_angle, outlet_metal_angle : float
        Camber-line angles at mean radius [rad].
    stagger : float
        Angle of the chord line from axial at mean radius [rad].
    leading_edge_diameter : float
        Leading-edge thickness [m].
    trailing_edge_thickness : float
        Trailing-edge thickness [m].
    max_thickness : float
        The thickness the section is drawn with [m]. The built section's
        actual maximum is measured and reported alongside.
    trailing_edge_wedge_angles : (float, float)
        Suction- and pressure-side wedge angles at the trailing edge [rad].
    thickness_distribution : ThicknessDistribution
        How the thickness is laid along the camber line.
    tip_clearance : float
        Gap from a rotor's tips to the casing [m]. Zero is a shrouded rotor,
        whose shroud turns with it and forms the outer wall over the row.
    hub_clearance : float
        Gap from a cantilevered stator's hub end to the rotating drum under
        it [m]. Zero is a stator carried on a stationary hub.
    admission : float
        Fraction of the circumference fed with flow.
    twist : callable, optional
        ``twist(span_fraction, r) -> (inlet_metal, outlet_metal)``. None gives
        the mean section at every span.
    """

    name: str
    kind: str
    blade_count: int
    x_le: float
    axial_chord: float
    r_hub_in: float
    r_tip_in: float
    r_hub_out: float
    r_tip_out: float
    inlet_metal_angle: float
    outlet_metal_angle: float
    stagger: float
    leading_edge_diameter: float
    trailing_edge_thickness: float
    max_thickness: float
    trailing_edge_wedge_angles: tuple = (math.radians(5), math.radians(2))
    thickness_distribution: ThicknessDistribution = ThicknessDistribution()
    tip_clearance: float = 0.0
    hub_clearance: float = 0.0
    admission: float = 1.0
    twist: object = None

    def __post_init__(self):
        label = f"row {self.name!r}"
        if self.kind not in KINDS:
            raise ValueError(f"{label}: kind must be one of {KINDS}")
        if isinstance(self.blade_count, bool) or int(self.blade_count) != self.blade_count \
                or self.blade_count < 2:
            raise ValueError(f"{label}: blade_count must be an integer >= 2")
        object.__setattr__(self, "blade_count", int(self.blade_count))
        if not np.isfinite(self.x_le):
            raise ValueError(f"{label}: x_le must be finite")
        for name in ("axial_chord", "r_hub_in", "r_tip_in", "r_hub_out", "r_tip_out",
                     "leading_edge_diameter", "trailing_edge_thickness", "max_thickness"):
            _positive(getattr(self, name), f"{label}: {name}")
        if self.r_tip_in <= self.r_hub_in or self.r_tip_out <= self.r_hub_out:
            raise ValueError(f"{label}: tip radius must be above hub radius")
        for name in ("inlet_metal_angle", "outlet_metal_angle", "stagger"):
            if not abs(getattr(self, name)) < _MAX_ANGLE:
                raise ValueError(f"{label}: {name} must lie within ±85° of axial")
        wedges = tuple(float(w) for w in self.trailing_edge_wedge_angles)
        if len(wedges) != 2 or not all(0 <= w < math.pi/2 for w in wedges):
            raise ValueError(f"{label}: trailing_edge_wedge_angles needs two angles in [0, 90°)")
        object.__setattr__(self, "trailing_edge_wedge_angles", wedges)
        if not isinstance(self.thickness_distribution, ThicknessDistribution):
            raise TypeError(f"{label}: thickness_distribution must be a ThicknessDistribution")
        for name in ("tip_clearance", "hub_clearance"):
            if not (np.isfinite(getattr(self, name)) and getattr(self, name) >= 0):
                raise ValueError(f"{label}: {name} must be nonnegative")
        if self.kind == "rotor" and self.hub_clearance:
            raise ValueError(f"{label}: a rotor is carried on its hub and has no hub clearance")
        if self.kind == "stator" and self.tip_clearance:
            raise ValueError(f"{label}: a stator is carried by the casing and has no tip clearance")
        if max(self.tip_clearance, self.hub_clearance) >= min(self.height_in, self.height_out):
            raise ValueError(f"{label}: clearance must be smaller than the blade height")
        if not 0 < self.admission <= 1:
            raise ValueError(f"{label}: admission must lie in (0, 1]")
        if self.twist is not None and not callable(self.twist):
            raise TypeError(f"{label}: twist must be callable or None")

    # -- derived hardware numbers ---------------------------------------------
    @property
    def x_te(self) -> float:
        return self.x_le + self.axial_chord

    @property
    def is_rotor(self) -> bool:
        return self.kind == "rotor"

    @property
    def chord(self) -> float:
        """True chord at mean radius [m]."""
        return self.axial_chord/math.cos(self.stagger)

    @property
    def height_in(self) -> float:
        return self.r_tip_in - self.r_hub_in

    @property
    def height_out(self) -> float:
        return self.r_tip_out - self.r_hub_out

    @property
    def r_mean(self) -> float:
        return 0.25*(self.r_hub_in + self.r_tip_in + self.r_hub_out + self.r_tip_out)

    @property
    def pitch(self) -> float:
        """Blade spacing at mean radius [m]."""
        return 2*math.pi*self.r_mean/self.blade_count

    @property
    def aspect_ratio(self) -> float:
        """Mean blade height over chord."""
        return 0.5*(self.height_in + self.height_out)/self.chord

    def section_angles(self, span_fraction, r):
        """Inlet metal, outlet metal and stagger at one span station [rad].

        The stagger keeps its offset from the mean of the two metal angles,
        so a twisted section leans with its camber.
        """
        if self.twist is None:
            return self.inlet_metal_angle, self.outlet_metal_angle, self.stagger
        inlet, outlet = (float(a) for a in self.twist(span_fraction, r))
        offset = self.stagger - 0.5*(self.inlet_metal_angle + self.outlet_metal_angle)
        return inlet, outlet, 0.5*(inlet + outlet) + offset


@dataclass(frozen=True)
class TurbineDefinition:
    """An axial turbine's rows and the annular passage through them.

    Parameters
    ----------
    rows : sequence of BladeRow
        In flow order, without overlapping.
    rotation : str
        Which way the rotor turns about +x. See :mod:`pyskyfire.common.rotation`.
    inlet_length, outlet_length : float
        Straight annulus ahead of the first row and behind the last [m]. The
        passage ends there, where a duct would join it.

    Attributes
    ----------
    passage : MeridionalPassage
        Hub and casing through all rows. The hub turns under each rotor, the
        casing turns over a shrouded one, and every clearance is recorded
        against the row it belongs to.
    """

    rows: tuple
    rotation: str = COUNTERCLOCKWISE
    inlet_length: float = 0.0
    outlet_length: float = 0.0
    passage: MeridionalPassage = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        rows = tuple(self.rows)
        if not rows or not all(isinstance(r, BladeRow) for r in rows):
            raise ValueError("rows must be a non-empty sequence of BladeRow")
        names = [r.name for r in rows]
        if len(set(names)) != len(names):
            raise ValueError("row names must be unique")
        for a, b in zip(rows, rows[1:]):
            if b.x_le < a.x_te - 1e-12:
                raise ValueError(f"rows {a.name!r} and {b.name!r} overlap axially")
        for name in ("inlet_length", "outlet_length"):
            if not (np.isfinite(getattr(self, name)) and getattr(self, name) >= 0):
                raise ValueError(f"{name} must be nonnegative")
        object.__setattr__(self, "rows", rows)
        object.__setattr__(self, "rotation", parse_rotation(self.rotation))
        object.__setattr__(self, "passage", self._build_passage())

    def _build_passage(self, points_per_segment=8):
        # Corner stations: every edge, plus the ends of the straight lengths.
        # Each row's walls run straight from its leading to its trailing edge,
        # and a gap between rows runs straight from one to the next.
        corners = []
        first, last = self.rows[0], self.rows[-1]
        if self.inlet_length > 0:
            corners.append((first.x_le - self.inlet_length, first.r_hub_in, first.r_tip_in))
        for row in self.rows:
            corners.append((row.x_le, row.r_hub_in, row.r_tip_in))
            corners.append((row.x_te, row.r_hub_out, row.r_tip_out))
        if self.outlet_length > 0:
            corners.append((last.x_te + self.outlet_length, last.r_hub_out, last.r_tip_out))
        points = [corners[0]]
        for a, b in zip(corners, corners[1:]):
            if b[0] - a[0] <= 1e-12:
                if not np.allclose(a[1:], b[1:]):
                    raise ValueError("rows with no gap between them must meet at the same radii")
                continue
            for t in np.linspace(0, 1, points_per_segment + 1)[1:]:
                points.append(tuple(np.add(a, t*np.subtract(b, a))))
        points = np.array(points)
        hub, shroud = points[:, [0, 1]], points[:, [0, 2]]
        base = MeridionalPassage(hub, shroud)

        def fraction(x):
            return float(np.interp(x, hub[:, 0], base.stations))

        hub_turns, shroud_turns, clearances = [], [], []
        for row in self.rows:
            span = (fraction(row.x_le), fraction(row.x_te))
            if row.is_rotor:
                hub_turns.append(span)
                if row.tip_clearance > 0:
                    clearances.append(Clearance(row.name, "shroud", *span, row.tip_clearance))
                else:
                    shroud_turns.append(span)
            elif row.hub_clearance > 0:
                hub_turns.append(span)
                clearances.append(Clearance(row.name, "hub", *span, row.hub_clearance))
        return replace(base, motion=WallMotion(hub=hub_turns, shroud=shroud_turns),
                       clearances=tuple(clearances))

    # -- lookups ------------------------------------------------------------
    @property
    def rotation_sign(self) -> float:
        """+1 counterclockwise about +x, -1 clockwise."""
        return rotation_sign(self.rotation)

    def row(self, name) -> BladeRow:
        for r in self.rows:
            if r.name == name:
                return r
        raise KeyError(name)

    @property
    def axial_gaps(self) -> tuple:
        """Axial spacing between consecutive rows [m]."""
        return tuple(b.x_le - a.x_te for a, b in zip(self.rows, self.rows[1:]))

    @property
    def length(self) -> float:
        """Axial length of the passage, including the straight ends [m]."""
        return float(self.passage.hub[-1, 0] - self.passage.hub[0, 0])

    @property
    def inlet_annulus(self):
        return self.passage.end("inlet")

    @property
    def exit_annulus(self):
        return self.passage.end("outlet")

    def build_geometry(self, *, span_points=9, chord_points=129):
        """Blade surfaces for every row, and each row's throat, measured from them.

        Needs the optional ``pyturbo-aero`` package.
        """
        from ._pyturbo import build_geometry
        return build_geometry(self, span_points=span_points, chord_points=chord_points)
