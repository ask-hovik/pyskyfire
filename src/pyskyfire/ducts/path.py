"""Shared centerline construction, rigid placement, and transported frames.

Paths contain geometry only. Flow areas and fluid states belong to components.
All vectors use Cartesian coordinates and all angles are radians.
"""

from __future__ import annotations
from dataclasses import dataclass, replace, field
import math
from typing import Sequence
import numpy as np

__all__ = [
    "PathNode",
    "PathStation",
    "Path",
    "Straight",
    "Continuous",
    "RoundedBend",
    "PathPlacement",
    "rotation_about_axis",
    "generate_path_stations",
]

_GEOMETRY_TOL = 1.0e-12
_ROUNDED_METHODS = (
    "Rennels",
    "Miller",
    "Crane",
    "Crane standard",
    "Ito",
    "Swamee",
)


def _position(value, *, name: str) -> tuple[float, float, float]:
    vector = np.asarray(value, dtype=float)
    if vector.shape != (3,) or not np.isfinite(vector).all():
        raise ValueError(f"{name} must be a finite three-dimensional position")
    return tuple(float(component) for component in vector)


def _positive(value, *, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


@dataclass(frozen=True)
class Straight:
    """Distributed straight-pipe friction between two adjacent stations."""


@dataclass(frozen=True)
class Continuous:
    """A locally resolved continuous-curvature connection.

    ``radius`` is the exact centerline radius when known. If omitted on a
    manually supplied layout, :class:`Path` estimates it transparently from
    neighboring station triples using Menger curvature. Attached to a
    :class:`PathNode`, it rounds that corner with a circular arc of ``radius``,
    which the generator emits as one connection carrying its exact arc
    ``length``, ``center`` and ``normal``. A constant-section arc needs no
    further stations; :class:`Pipe` adds them where its section varies.
    """

    radius: float | None = None
    length: float | None = None
    center: tuple[float, float, float] | None = None
    normal: tuple[float, float, float] | None = None

    def __post_init__(self):
        if self.radius is not None:
            object.__setattr__(self, "radius", _positive(self.radius, name="radius"))
        if self.length is not None:
            object.__setattr__(self, "length", _positive(self.length, name="length"))
        if (self.center is None) != (self.normal is None):
            raise ValueError("continuous arc geometry requires both center and normal")
        if self.center is not None:
            if self.radius is None or self.length is None:
                raise ValueError("continuous arc geometry requires radius and length")
            object.__setattr__(self, "center", _position(self.center, name="center"))
            normal = np.asarray(_position(self.normal, name="normal"))
            norm = np.linalg.norm(normal)
            if norm <= _GEOMETRY_TOL:
                raise ValueError("normal must have nonzero magnitude")
            object.__setattr__(self, "normal", tuple(normal / norm))


@dataclass(frozen=True)
class RoundedBend:
    """A complete lumped rounded-bend loss between adjacent stations.

    Parameters
    ----------
    radius : float
        Bend centerline radius [m].
    method : str, optional
        A method accepted by :func:`fluids.fittings.bend_rounded`.
    angle : float, optional
        Resolved turn angle [rad]. Required by :class:`Pipe`; the station
        generator calculates it from the neighboring control legs.
    center, normal : tuple, optional
        Resolved arc geometry retained for inspection and visualization.
    """

    radius: float
    method: str = "Rennels"
    angle: float | None = None
    center: tuple[float, float, float] | None = None
    normal: tuple[float, float, float] | None = None

    def __post_init__(self):
        object.__setattr__(self, "radius", _positive(self.radius, name="radius"))
        if self.method not in _ROUNDED_METHODS:
            choices = ", ".join(_ROUNDED_METHODS)
            raise ValueError(
                f"unknown rounded-bend method {self.method!r}; choose from {choices}"
            )
        if self.angle is not None:
            angle = float(self.angle)
            if not np.isfinite(angle) or not 0.0 < angle < math.pi:
                raise ValueError("rounded-bend angle must be between 0 and pi radians")
            object.__setattr__(self, "angle", angle)
        if self.center is not None:
            object.__setattr__(self, "center", _position(self.center, name="center"))
        if self.normal is not None:
            normal = np.asarray(_position(self.normal, name="normal"))
            magnitude = np.linalg.norm(normal)
            if magnitude <= _GEOMETRY_TOL:
                raise ValueError("normal must have nonzero magnitude")
            object.__setattr__(self, "normal", tuple(normal / magnitude))


@dataclass(frozen=True)
class PathNode:
    """Human-friendly route vertex consumed by :func:`generate_path_stations`."""

    position: tuple[float, float, float]
    bend: RoundedBend | Continuous | None = None

    def __post_init__(self):
        object.__setattr__(
            self, "position", _position(self.position, name="PathNode.position")
        )
        if self.bend is not None and not isinstance(
            self.bend, (RoundedBend, Continuous)
        ):
            raise TypeError("PathNode.bend must be RoundedBend, Continuous, or None")
        if isinstance(self.bend, Continuous) and self.bend.radius is None:
            raise ValueError("a Continuous bend on a PathNode requires a radius")


@dataclass(frozen=True)
class PathStation:
    """A physical centerline station shared by flow components."""

    position: tuple[float, float, float]

    def __post_init__(self):
        object.__setattr__(
            self, "position", _position(self.position, name="PathStation.position")
        )


@dataclass(frozen=True)
class _CornerGeometry:
    incoming: np.ndarray
    outgoing: np.ndarray
    tangent_in: np.ndarray
    tangent_out: np.ndarray
    center: np.ndarray
    normal: np.ndarray
    angle: float
    tangent_distance: float


def _corner_geometry(
    previous: PathNode, node: PathNode, following: PathNode
) -> _CornerGeometry:
    vertex = np.asarray(node.position)
    incoming_leg = vertex - np.asarray(previous.position)
    outgoing_leg = np.asarray(following.position) - vertex
    incoming_length = np.linalg.norm(incoming_leg)
    outgoing_length = np.linalg.norm(outgoing_leg)
    if incoming_length <= _GEOMETRY_TOL or outgoing_length <= _GEOMETRY_TOL:
        raise ValueError("path nodes must not contain consecutive duplicate positions")

    incoming = incoming_leg / incoming_length
    outgoing = outgoing_leg / outgoing_length
    cosine = float(np.clip(np.dot(incoming, outgoing), -1.0, 1.0))
    angle = math.acos(cosine)
    if angle <= 1.0e-10:
        raise ValueError("a bend cannot be attached to a straight-through path node")
    if math.pi - angle <= 1.0e-8:
        raise ValueError(
            "a rounded or continuous bend cannot resolve a reversing path node"
        )

    radius = node.bend.radius
    tangent_distance = radius * math.tan(0.5 * angle)
    if tangent_distance >= incoming_length - _GEOMETRY_TOL:
        raise ValueError(
            f"bend radius {radius:g} m at {node.position} requires tangent distance "
            f"{tangent_distance:g} m, but the incoming leg is only {incoming_length:g} m"
        )
    if tangent_distance >= outgoing_length - _GEOMETRY_TOL:
        raise ValueError(
            f"bend radius {radius:g} m at {node.position} requires tangent distance "
            f"{tangent_distance:g} m, but the outgoing leg is only {outgoing_length:g} m"
        )

    sine = math.sin(angle)
    radial_inward = (outgoing - cosine * incoming) / sine
    tangent_in = vertex - tangent_distance * incoming
    tangent_out = vertex + tangent_distance * outgoing
    center = tangent_in + radius * radial_inward
    normal = np.cross(incoming, outgoing) / sine

    if np.dot(tangent_in - np.asarray(previous.position), incoming) <= 0.0:
        raise ValueError("bend incoming tangent lies backwards on its control leg")
    if np.dot(np.asarray(following.position) - tangent_out, outgoing) <= 0.0:
        raise ValueError("bend outgoing tangent lies backwards on its control leg")

    return _CornerGeometry(
        incoming,
        outgoing,
        tangent_in,
        tangent_out,
        center,
        normal,
        angle,
        tangent_distance,
    )


def _rotate(vector: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotate ``vector`` around unit ``axis`` using Rodrigues' formula."""
    return (
        vector * math.cos(angle)
        + np.cross(axis, vector) * math.sin(angle)
        + axis * np.dot(axis, vector) * (1.0 - math.cos(angle))
    )


def _arc_angle(connection) -> float:
    """Turn angle [rad] of a resolved arc connection."""
    if isinstance(connection, RoundedBend):
        return connection.angle
    return connection.length / connection.radius


def _bend_arc(start, end, bend, segments: int) -> list[np.ndarray]:
    """Points along a resolved arc connection, excluding its starting station."""
    center = np.asarray(bend.center)
    normal = np.asarray(bend.normal)
    radial = np.asarray(start, dtype=float) - center
    angle = _arc_angle(bend)
    arc = [
        center + _rotate(radial, normal, angle * step / segments)
        for step in range(1, segments + 1)
    ]
    tolerance = 1.0e-6 * max(bend.radius, float(np.linalg.norm(radial)))
    if (
        abs(float(np.linalg.norm(radial)) - bend.radius) > tolerance
        or float(np.linalg.norm(arc[-1] - np.asarray(end, dtype=float))) > tolerance
    ):
        raise ValueError(
            "rounded-bend center and normal are inconsistent with its tangent stations"
        )
    arc[-1] = np.asarray(end, dtype=float)
    return arc


def generate_path_stations(
    nodes: Sequence[PathNode],
) -> tuple[
    tuple[PathStation, ...],
    tuple[Straight | Continuous | RoundedBend, ...],
]:
    """Resolve control vertices into physical stations and connections.

    Interior direction changes must state their interpretation explicitly with
    :class:`RoundedBend` or :class:`Continuous`. Rounded bends become one
    lumped connection between tangent points. Continuous bends become one
    connection carrying the exact arc length, center and normal.
    """
    nodes = tuple(nodes)
    if len(nodes) < 2:
        raise ValueError("at least two path nodes are required")
    if not all(isinstance(node, PathNode) for node in nodes):
        raise TypeError("nodes must contain only PathNode instances")
    if nodes[0].bend is not None or nodes[-1].bend is not None:
        raise ValueError("the first and last path nodes cannot define a bend")

    for first, second in zip(nodes, nodes[1:]):
        if (
            np.linalg.norm(np.asarray(second.position) - np.asarray(first.position))
            <= _GEOMETRY_TOL
        ):
            raise ValueError(
                "path nodes must not contain consecutive duplicate positions"
            )

    corners: dict[int, _CornerGeometry] = {}
    for index in range(1, len(nodes) - 1):
        node = nodes[index]
        before = np.asarray(node.position) - np.asarray(nodes[index - 1].position)
        after = np.asarray(nodes[index + 1].position) - np.asarray(node.position)
        before /= np.linalg.norm(before)
        after /= np.linalg.norm(after)
        turn = math.acos(float(np.clip(np.dot(before, after), -1.0, 1.0)))
        if node.bend is None:
            if turn > 1.0e-10:
                if math.pi - turn <= 1.0e-8:
                    raise ValueError("the path reverses at an unresolved interior node")
                raise ValueError(
                    f"non-collinear path node {index} must define RoundedBend or Continuous"
                )
        else:
            corners[index] = _corner_geometry(nodes[index - 1], node, nodes[index + 1])

    for index in range(len(nodes) - 1):
        leg_length = np.linalg.norm(
            np.asarray(nodes[index + 1].position) - np.asarray(nodes[index].position)
        )
        trim_start = corners[index].tangent_distance if index in corners else 0.0
        trim_end = corners[index + 1].tangent_distance if index + 1 in corners else 0.0
        if trim_start + trim_end >= leg_length - _GEOMETRY_TOL:
            raise ValueError(
                f"bends adjacent to control leg {index} overlap or touch; required trims "
                f"sum to {trim_start + trim_end:g} m on a {leg_length:g} m leg"
            )

    stations: list[PathStation] = [PathStation(nodes[0].position)]
    connections: list[Straight | Continuous | RoundedBend] = []

    def append_station(position, connection):
        station = PathStation(tuple(position))
        if (
            np.linalg.norm(
                np.asarray(station.position) - np.asarray(stations[-1].position)
            )
            <= _GEOMETRY_TOL
        ):
            raise ValueError("resolved pipe contains coincident adjacent stations")
        connections.append(connection)
        stations.append(station)

    for index in range(1, len(nodes) - 1):
        node = nodes[index]
        geometry = corners.get(index)
        if geometry is None:
            append_station(node.position, Straight())
            continue

        append_station(geometry.tangent_in, Straight())
        if isinstance(node.bend, RoundedBend):
            resolved = replace(
                node.bend,
                angle=geometry.angle,
                center=tuple(geometry.center),
                normal=tuple(geometry.normal),
            )
            append_station(geometry.tangent_out, resolved)
        else:
            connection = Continuous(
                radius=node.bend.radius,
                length=node.bend.radius * geometry.angle,
                center=tuple(geometry.center),
                normal=tuple(geometry.normal),
            )
            append_station(geometry.tangent_out, connection)

    append_station(nodes[-1].position, Straight())
    return tuple(stations), tuple(connections)


# Uniform parameter intervals sampled before adaptive halving begins. A feature
# that turns and returns within one of them can be missed.
_CURVE_INITIAL_INTERVALS = 16


def _evaluate_curve(curve, t: np.ndarray) -> np.ndarray:
    """Evaluate a curve-like object at parameters ``t``, returning ``(len(t), 3)``.

    The curve is first called with the whole array, as SciPy splines and
    vectorized functions expect, and otherwise once per parameter value. Planar
    ``(len(t), 2)`` results are taken to lie in the x-y plane.
    """
    t = np.asarray(t, dtype=float)

    def checked(points):
        points = np.asarray(points, dtype=float)
        if points.shape == (len(t), 2):
            points = np.column_stack((points, np.zeros(len(t))))
        return points if points.shape == (len(t), 3) else None

    error = None
    try:
        points = checked(curve(t))
    except Exception as exc:  # noqa: BLE001 -- retried point by point below
        points, error = None, exc
    if points is None:
        try:
            points = checked([np.asarray(curve(float(x)), dtype=float).ravel() for x in t])
        except Exception as exc:
            raise TypeError(
                "curve must be callable on a parameter value or an array of them"
            ) from (error or exc)
    if points is None:
        raise ValueError(
            "curve must return points of shape (len(t), 3) or (len(t), 2) for an "
            "array t, or one 2D or 3D point for a scalar t; wrap other forms, "
            "such as splev output, as lambda t: np.array(splev(t, tck)).T"
        )
    if not np.isfinite(points).all():
        raise ValueError("curve returned a non-finite point")
    return points


def _half_turn(start, middle, end) -> np.ndarray:
    """Angle [rad] between the two half-chords of each interval."""
    u, v = middle - start, end - middle
    return np.arctan2(np.linalg.norm(np.cross(u, v), axis=1), np.einsum("ij,ij->i", u, v))


def _circle_connection(start, middle, end):
    """The circular arc through three points, as one connection.

    Straight when the points are collinear. Otherwise a :class:`Continuous`
    carrying the circle's exact radius, arc length, center, and normal.
    """
    a, b = start - end, middle - end
    normal = np.cross(a, b)
    twice_area = np.linalg.norm(normal)
    if twice_area <= 1.0e-12 * np.linalg.norm(a) * np.linalg.norm(b):
        return Straight()
    center = end + np.cross(np.dot(a, a) * b - np.dot(b, b) * a, normal) / (2.0 * twice_area**2)
    radius = float(np.linalg.norm(start - center))
    chord = float(np.linalg.norm(end - start))
    angle = 2.0 * math.asin(min(1.0, chord / (2.0 * radius)))
    turn = np.cross(middle - start, end - middle)
    return Continuous(radius=radius, length=radius * angle, center=tuple(center),
                      normal=tuple(turn / np.linalg.norm(turn)))


def _menger_radius(a, b, c) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    ab = np.linalg.norm(b - a)
    bc = np.linalg.norm(c - b)
    ac = np.linalg.norm(c - a)
    twice_area = np.linalg.norm(np.cross(b - a, c - a))
    if twice_area <= _GEOMETRY_TOL * max(ab * bc, 1.0):
        return math.inf
    return ab * bc * ac / (2.0 * twice_area)


def _estimated_connection_radius(
    stations: tuple[PathStation, ...], index: int
) -> float:
    candidates = []
    if index > 0:
        candidates.append(
            _menger_radius(
                stations[index - 1].position,
                stations[index].position,
                stations[index + 1].position,
            )
        )
    if index + 2 < len(stations):
        candidates.append(
            _menger_radius(
                stations[index].position,
                stations[index + 1].position,
                stations[index + 2].position,
            )
        )
    finite = [radius for radius in candidates if np.isfinite(radius)]
    return float(np.mean(finite)) if finite else math.inf


@dataclass(frozen=True)
class Path:
    """Resolved stations and connections in order from inlet to outlet.

    ``s`` is true arc length [m]. Geometry carries no bore diameter or state.
    Use :meth:`refined` to choose the thermal/flow integration resolution.
    """

    stations: tuple[PathStation, ...]
    connections: tuple[Straight | Continuous | RoundedBend, ...]

    def __post_init__(self):
        object.__setattr__(self, "stations", tuple(self.stations))
        object.__setattr__(self, "connections", tuple(self.connections))
        if len(self.stations) < 2 or not all(
            isinstance(s, PathStation) for s in self.stations
        ):
            raise ValueError("Path requires at least two PathStation instances")
        if len(self.connections) != len(self.stations) - 1:
            raise ValueError(
                "Path requires exactly one connection between each station pair"
            )
        if not all(
            isinstance(c, (Straight, Continuous, RoundedBend)) for c in self.connections
        ):
            raise TypeError("unsupported path connection type")
        if np.any(
            np.linalg.norm(np.diff(self.centerline, axis=0), axis=1) <= _GEOMETRY_TOL
        ):
            raise ValueError(
                "Path stations must not contain consecutive duplicate positions"
            )
        for c in self.connections:
            if isinstance(c, RoundedBend) and c.angle is None:
                raise ValueError("RoundedBend connections require a resolved angle")

    @classmethod
    def from_nodes(cls, nodes):
        return cls(*generate_path_stations(nodes))

    @classmethod
    def from_curve(cls, curve, t_span, *, max_turn: float = math.radians(5.0)):
        """Sample a curve-like object into stations, as finely as its bends need.

        ``curve`` is anything that maps a parameter ``t`` to a point: called
        with an array of parameters it returns points of shape ``(len(t), 3)``,
        or it returns one point per scalar ``t``. SciPy's ``CubicSpline``,
        ``make_interp_spline`` and ``PchipInterpolator`` qualify as they are, and
        so does a plain function. Planar points ``(len(t), 2)`` lie in the x-y
        plane. ``t_span`` is the ``(start, end)`` parameter range to follow.

        Sampling is adaptive. Any interval whose curve turns by more than
        ``max_turn`` [rad] is halved, so stations crowd into tight bends and
        stay sparse along straight stretches. Each interval then becomes the
        circular arc through its two stations and the curve's own midpoint: a
        :class:`Continuous` carrying that arc's radius, length, center, and
        normal, or a :class:`Straight` where the three are collinear.

        Raises
        ------
        ValueError
            If the curve has a corner, a sharp turn no amount of halving
            resolves, or does not advance between two parameter values.
        """
        start_t, end_t = (float(x) for x in t_span)
        if not np.isfinite([start_t, end_t]).all() or start_t == end_t:
            raise ValueError("t_span must be two distinct finite parameter values")
        max_turn = float(max_turn)
        if not 0.0 < max_turn <= math.pi / 2:
            raise ValueError("max_turn must be between 0 and pi/2 radians")
        smallest = 1.0e-9 * abs(end_t - start_t)

        t = np.linspace(start_t, end_t, _CURVE_INITIAL_INTERVALS + 1)
        points = _evaluate_curve(curve, t)
        while True:
            mid_t = 0.5 * (t[:-1] + t[1:])
            middle = _evaluate_curve(curve, mid_t)
            # Turn inside each interval, and across each station between the
            # half-chords on either side, which is where a corner that falls
            # exactly on a sample would otherwise hide.
            split = 2.0 * _half_turn(points[:-1], middle, points[1:]) > max_turn
            across = 2.0 * _half_turn(middle[:-1], points[1:-1], middle[1:]) > max_turn
            split[:-1] |= across
            split[1:] |= across
            if not split.any():
                break
            if np.any(np.abs(np.diff(t))[split] < smallest):
                where = mid_t[split][0]
                raise ValueError(
                    f"curve has a corner near t={where:g} that a pipe cannot follow; "
                    "round it with a bend or smooth the curve"
                )
            order = np.argsort(np.concatenate((t, mid_t[split])), kind="stable")
            t = np.concatenate((t, mid_t[split]))[order]
            points = np.concatenate((points, middle[split]))[order]

        chords = np.linalg.norm(np.diff(points, axis=0), axis=1)
        if np.any(chords <= _GEOMETRY_TOL):
            where = mid_t[np.argmax(chords <= _GEOMETRY_TOL)]
            raise ValueError(f"curve does not advance near t={where:g}")
        connections = [
            _circle_connection(a, m, b)
            for a, m, b in zip(points[:-1], middle, points[1:])
        ]
        return cls(tuple(PathStation(tuple(point)) for point in points), connections)

    @property
    def centerline(self):
        return np.asarray([s.position for s in self.stations], float)

    @property
    def lengths(self):
        chord = np.linalg.norm(np.diff(self.centerline, axis=0), axis=1)
        return np.asarray(
            [
                c.radius * c.angle
                if isinstance(c, RoundedBend)
                else c.length
                if isinstance(c, Continuous) and c.length is not None
                else chord[i]
                for i, c in enumerate(self.connections)
            ]
        )

    @property
    def s(self):
        return np.r_[0.0, np.cumsum(self.lengths)]

    @property
    def curvature_radius(self):
        return np.asarray(
            [
                np.inf
                if isinstance(c, Straight)
                else c.radius
                if c.radius is not None
                else _estimated_connection_radius(self.stations, i)
                for i, c in enumerate(self.connections)
            ]
        )

    @property
    def tangents(self):
        points = self.centerline
        chord = np.diff(points, axis=0)
        chord /= np.linalg.norm(chord, axis=1)[:, None]
        start, end = chord.copy(), chord.copy()
        for i, c in enumerate(self.connections):
            if (
                getattr(c, "center", None) is not None
                and getattr(c, "normal", None) is not None
            ):
                start[i] = np.cross(c.normal, points[i] - c.center)
                end[i] = np.cross(c.normal, points[i + 1] - c.center)
                start[i] /= np.linalg.norm(start[i])
                end[i] /= np.linalg.norm(end[i])
        tangent = np.vstack((start[0], end[:-1] + start[1:], end[-1]))
        norm = np.linalg.norm(tangent, axis=1)
        if np.any(norm < 1e-10):
            raise ValueError("a reversing path cannot define transported frames")
        return tangent / norm[:, None]

    def refined(self, max_step):
        """Split distributed connections, retaining exact generated arc geometry.

        Manually supplied Continuous connections without a center and normal
        must already be resolved finely enough; their centerline is interpolated
        linearly while their supplied arc length and curvature are retained.
        """
        max_step = _positive(max_step, name="max_step")
        stations = [self.stations[0]]
        connections = []
        points = self.centerline
        for i, (c, length) in enumerate(zip(self.connections, self.lengths)):
            if isinstance(c, RoundedBend):
                connections.append(c)
                stations.append(self.stations[i + 1])
                continue
            count = max(1, math.ceil(length / max_step))
            for j in range(1, count + 1):
                fraction = j / count
                if isinstance(c, Continuous) and c.center is not None:
                    point = np.asarray(c.center) + _rotate(
                        points[i] - c.center,
                        np.asarray(c.normal),
                        fraction * length / c.radius,
                    )
                else:
                    point = points[i] + fraction * (points[i + 1] - points[i])
                stations.append(
                    PathStation(tuple(points[i + 1] if j == count else point))
                )
                connections.append(
                    replace(c, length=length / count)
                    if isinstance(c, Continuous)
                    else c
                )
        return Path(stations, connections)


def rotation_about_axis(axis, angle):
    """Return a proper 3x3 rotation; ``angle`` is in radians."""
    axis = np.asarray(_position(axis, name="axis"))
    norm = np.linalg.norm(axis)
    if norm <= _GEOMETRY_TOL or not np.isfinite(angle):
        raise ValueError("rotation requires a nonzero axis and finite angle")
    axis /= norm
    return np.column_stack([_rotate(e, axis, angle) for e in np.eye(3)])


@dataclass(frozen=True)
class PathPlacement:
    """Place the first path station and establish its section roll.

    The rotation maps template coordinates into assembly coordinates. The
    initial floor direction is expressed in assembly coordinates, points from
    the passage centroid toward its floor, and is projected perpendicular to
    the starting tangent. Frame columns are (floor-to-roof, lateral, tangent).
    All subsequent section orientations use minimal-rotation transport.
    """

    start_position: tuple[float, float, float] = (0.0, 0.0, 0.0)
    rotation: np.ndarray = field(default_factory=lambda: np.eye(3))
    initial_floor_direction: tuple[float, float, float] = (0.0, -1.0, 0.0)

    def __post_init__(self):
        object.__setattr__(
            self,
            "start_position",
            _position(self.start_position, name="start_position"),
        )
        object.__setattr__(
            self,
            "initial_floor_direction",
            _position(self.initial_floor_direction, name="initial_floor_direction"),
        )
        r = np.array(self.rotation, float, copy=True)
        if (
            r.shape != (3, 3)
            or not np.isfinite(r).all()
            or not np.allclose(r.T @ r, np.eye(3), atol=1e-10, rtol=0)
            or not np.isclose(np.linalg.det(r), 1, atol=1e-10, rtol=0)
        ):
            raise ValueError("rotation must be a proper orthonormal 3x3 matrix")
        r.setflags(write=False)
        object.__setattr__(self, "rotation", r)

    def place(self, path):
        points = (
            path.centerline - path.centerline[0]
        ) @ self.rotation.T + self.start_position
        connections = []
        for c in path.connections:
            if getattr(c, "center", None) is not None:
                c = replace(
                    c,
                    center=tuple(
                        self.rotation @ (np.asarray(c.center) - path.centerline[0])
                        + self.start_position
                    ),
                    normal=tuple(self.rotation @ c.normal),
                )
            connections.append(c)
        return Path(tuple(PathStation(tuple(p)) for p in points), connections)

    def frames(self, path):
        tangent = self.place(path).tangents
        d = np.asarray(self.initial_floor_direction)
        floor = d - np.dot(d, tangent[0]) * tangent[0]
        if np.linalg.norm(floor) <= 1e-10 * max(1.0, np.linalg.norm(d)):
            raise ValueError(
                "initial_floor_direction must not be zero or parallel to the path tangent"
            )
        up = -floor / np.linalg.norm(floor)
        frames = []
        for i, t in enumerate(tangent):
            if i:
                axis = np.cross(tangent[i - 1], t)
                sine = np.linalg.norm(axis)
                cosine = np.clip(np.dot(tangent[i - 1], t), -1.0, 1.0)
                if cosine < -1 + 1e-10:
                    raise ValueError(
                        "a reversing path cannot define transported frames"
                    )
                if sine > 1e-12:
                    up = _rotate(up, axis / sine, math.atan2(sine, cosine))
                up -= np.dot(up, t) * t
                up /= np.linalg.norm(up)
            frames.append(np.column_stack((up, np.cross(t, up), t)))
        return np.asarray(frames)
