"""Explicit pipe geometry construction and pressure-loss marching.

A :class:`Pipe` contains the same :class:`PipeStation` and connection objects
that its solver uses. :func:`generate_pipe_stations` is only a convenience for
turning human-friendly corner vertices into tangent stations and connections;
users may construct the resolved objects directly.

The first solver version marches irreversible wall and bend losses. It is
adiabatic, so station temperature is unchanged, and it deliberately does not
apply the reversible acceleration term associated with changing area or
compressibility. Consequently ``PipeResult.p`` is the loss-marched pressure
used by the engine network, not a complete variable-area static-pressure
solution.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Callable, Sequence
import warnings

import fluids
import numpy as np

from pyskyfire.thermo.state import FluidState, PHASE_UNKNOWN

from .friction import curvature_multiplier, darcy_friction_factor

__all__ = [
    "PathNode",
    "PipeStation",
    "Straight",
    "Continuous",
    "RoundedBend",
    "Pipe",
    "PipeResult",
    "PipeCorrelationWarning",
    "DisplayPath",
    "generate_pipe_stations",
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


class PipeCorrelationWarning(RuntimeWarning):
    """A pipe correlation is being evaluated outside a documented range."""


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
    manually supplied layout, :class:`Pipe` estimates it transparently from
    neighboring station triples using Menger curvature. ``segments`` controls
    circular-arc generation when this object is attached to a :class:`PathNode`.
    Generated local connections carry their exact arc ``length``.
    """

    radius: float | None = None
    segments: int = 16
    length: float | None = None

    def __post_init__(self):
        if self.radius is not None:
            object.__setattr__(self, "radius", _positive(self.radius, name="radius"))
        if isinstance(self.segments, bool) or not isinstance(self.segments, int) or self.segments < 1:
            raise ValueError("segments must be a positive integer")
        if self.length is not None:
            object.__setattr__(self, "length", _positive(self.length, name="length"))


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
            raise ValueError(f"unknown rounded-bend method {self.method!r}; choose from {choices}")
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
    """Human-friendly route vertex consumed by :func:`generate_pipe_stations`."""

    position: tuple[float, float, float]
    diameter: float
    bend: RoundedBend | Continuous | None = None

    def __post_init__(self):
        object.__setattr__(self, "position", _position(self.position, name="PathNode.position"))
        object.__setattr__(self, "diameter", _positive(self.diameter, name="PathNode.diameter"))
        if self.bend is not None and not isinstance(self.bend, (RoundedBend, Continuous)):
            raise TypeError("PathNode.bend must be RoundedBend, Continuous, or None")
        if isinstance(self.bend, Continuous) and self.bend.radius is None:
            raise ValueError("a Continuous bend on a PathNode requires a radius")


@dataclass(frozen=True)
class PipeStation:
    """A physical centerline station used directly by the pipe solver."""

    position: tuple[float, float, float]
    diameter: float

    def __post_init__(self):
        object.__setattr__(self, "position", _position(self.position, name="PipeStation.position"))
        object.__setattr__(self, "diameter", _positive(self.diameter, name="PipeStation.diameter"))


@dataclass(frozen=True)
class PipeResult:
    """Station states and per-connection diagnostics from a pipe loss march."""

    stations: tuple[PipeStation, ...]
    connections: tuple[Straight | Continuous | RoundedBend, ...]
    states: tuple[FluidState, ...]
    s: np.ndarray
    p: np.ndarray
    T: np.ndarray
    density: np.ndarray
    velocity: np.ndarray
    pressure_loss: np.ndarray
    reynolds: np.ndarray
    darcy_factor: np.ndarray
    loss_coefficient: np.ndarray

    @property
    def total_pressure_loss(self) -> float:
        """Total irreversible pressure loss [Pa]."""
        return float(self.p[0] - self.p[-1])


@dataclass(frozen=True)
class DisplayPath:
    """A refined centerline for drawing, never for pressure-loss marching.

    ``points`` [m] and ``radius`` [m] describe the swept circular tube. The
    path carries its own provenance so results can follow it: ``station`` is
    the fractional :class:`PipeStation` index of every point, and
    ``connection`` is the originating connection index of every segment
    between consecutive points.
    """

    points: np.ndarray
    radius: np.ndarray
    station: np.ndarray
    connection: np.ndarray

    def __post_init__(self):
        points = np.array(self.points, dtype=float, copy=True)
        radius = np.array(self.radius, dtype=float, copy=True)
        station = np.array(self.station, dtype=float, copy=True)
        connection = np.array(self.connection, dtype=np.int64, copy=True)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("DisplayPath points must be finite with shape (M,3)")
        if radius.shape != (len(points),) or station.shape != (len(points),):
            raise ValueError("DisplayPath radius and station must have one value per point")
        if connection.shape != (max(len(points) - 1, 0),):
            raise ValueError("DisplayPath connection must have one value per segment")
        for array in (points, radius, station, connection):
            array.setflags(write=False)
        object.__setattr__(self, "points", points)
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "station", station)
        object.__setattr__(self, "connection", connection)

    @property
    def connection_count(self) -> int:
        """Number of source connections the path was refined from."""
        return int(self.connection[-1]) + 1 if len(self.connection) else 0

    @property
    def station_count(self) -> int:
        """Number of source :class:`PipeStation` objects the path came from."""
        return self.connection_count + 1 if len(self.connection) else len(self.points)


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


def _corner_geometry(previous: PathNode, node: PathNode, following: PathNode) -> _CornerGeometry:
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
        raise ValueError("a rounded or continuous bend cannot resolve a reversing path node")

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
        incoming, outgoing, tangent_in, tangent_out, center, normal, angle,
        tangent_distance,
    )


def _rotate(vector: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotate ``vector`` around unit ``axis`` using Rodrigues' formula."""
    return (
        vector * math.cos(angle)
        + np.cross(axis, vector) * math.sin(angle)
        + axis * np.dot(axis, vector) * (1.0 - math.cos(angle))
    )


def _bend_arc(start, end, bend: RoundedBend, segments: int) -> list[np.ndarray]:
    """Points along a resolved rounded bend, excluding its starting station."""
    center = np.asarray(bend.center)
    normal = np.asarray(bend.normal)
    radial = np.asarray(start, dtype=float) - center
    arc = [
        center + _rotate(radial, normal, bend.angle * step / segments)
        for step in range(1, segments + 1)
    ]
    tolerance = 1.0e-6 * max(bend.radius, float(np.linalg.norm(radial)))
    if (abs(float(np.linalg.norm(radial)) - bend.radius) > tolerance
            or float(np.linalg.norm(arc[-1] - np.asarray(end, dtype=float))) > tolerance):
        raise ValueError(
            "rounded-bend center and normal are inconsistent with its tangent stations"
        )
    arc[-1] = np.asarray(end, dtype=float)
    return arc


def generate_pipe_stations(
    nodes: Sequence[PathNode],
) -> tuple[
    tuple[PipeStation, ...],
    tuple[Straight | Continuous | RoundedBend, ...],
]:
    """Resolve control vertices into physical stations and connections.

    Interior direction changes must state their interpretation explicitly with
    :class:`RoundedBend` or :class:`Continuous`. Rounded bends become one
    lumped connection between tangent points. Continuous bends become a
    circular arc discretized according to ``Continuous.segments``.
    """
    nodes = tuple(nodes)
    if len(nodes) < 2:
        raise ValueError("at least two path nodes are required")
    if not all(isinstance(node, PathNode) for node in nodes):
        raise TypeError("nodes must contain only PathNode instances")
    if nodes[0].bend is not None or nodes[-1].bend is not None:
        raise ValueError("the first and last path nodes cannot define a bend")

    for first, second in zip(nodes, nodes[1:]):
        if np.linalg.norm(np.asarray(second.position) - np.asarray(first.position)) <= _GEOMETRY_TOL:
            raise ValueError("path nodes must not contain consecutive duplicate positions")

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

    stations: list[PipeStation] = [PipeStation(nodes[0].position, nodes[0].diameter)]
    connections: list[Straight | Continuous | RoundedBend] = []

    def append_station(position, diameter, connection):
        station = PipeStation(tuple(position), diameter)
        if np.linalg.norm(np.asarray(station.position) - np.asarray(stations[-1].position)) <= _GEOMETRY_TOL:
            raise ValueError("resolved pipe contains coincident adjacent stations")
        connections.append(connection)
        stations.append(station)

    for index in range(1, len(nodes) - 1):
        node = nodes[index]
        geometry = corners.get(index)
        if geometry is None:
            append_station(node.position, node.diameter, Straight())
            continue

        append_station(geometry.tangent_in, node.diameter, Straight())
        if isinstance(node.bend, RoundedBend):
            resolved = replace(
                node.bend,
                angle=geometry.angle,
                center=tuple(geometry.center),
                normal=tuple(geometry.normal),
            )
            append_station(geometry.tangent_out, node.diameter, resolved)
        else:
            radial_start = geometry.tangent_in - geometry.center
            segment_angle = geometry.angle / node.bend.segments
            segment_length = node.bend.radius * segment_angle
            for segment in range(1, node.bend.segments + 1):
                radial = _rotate(radial_start, geometry.normal, segment * segment_angle)
                point = geometry.center + radial
                connection = Continuous(
                    radius=node.bend.radius,
                    segments=1,
                    length=segment_length,
                )
                append_station(point, node.diameter, connection)

    append_station(nodes[-1].position, nodes[-1].diameter, Straight())
    return tuple(stations), tuple(connections)


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


def _estimated_connection_radius(stations: tuple[PipeStation, ...], index: int) -> float:
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


def _property(transport, name: str, T: float, p: float) -> float:
    getter = getattr(transport, f"get_{name}", None)
    if getter is None:
        getter = getattr(transport, name, None)
    if getter is None or not callable(getter):
        raise TypeError(f"transport must provide get_{name}(T, p) or {name}(T, p)")
    value = float(getter(T, p))
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"transport returned invalid {name}={value!r}")
    return value


def _state(T: float, p: float, rho: float) -> FluidState:
    return FluidState(
        T=float(T),
        p=float(p),
        h=math.nan,
        quality=math.nan,
        rho=float(rho),
        phase=PHASE_UNKNOWN,
        T_sat=math.nan,
    )


class Pipe:
    """An explicitly resolved circular pipe and its pressure-loss solver."""

    def __init__(self, stations, connections, *, roughness: float | Callable[[float], float] = 0.0):
        self.stations = tuple(stations)
        self.connections = tuple(connections)
        if len(self.stations) < 2:
            raise ValueError("Pipe requires at least two stations")
        if not all(isinstance(station, PipeStation) for station in self.stations):
            raise TypeError("stations must contain only PipeStation instances")
        if len(self.connections) != len(self.stations) - 1:
            raise ValueError("Pipe requires exactly one connection between each station pair")
        if not all(isinstance(item, (Straight, Continuous, RoundedBend)) for item in self.connections):
            raise TypeError("unsupported pipe connection type")

        lengths = np.linalg.norm(np.diff(self.centerline, axis=0), axis=1)
        if np.any(lengths <= _GEOMETRY_TOL):
            raise ValueError("Pipe stations must not contain consecutive duplicate positions")

        for index, connection in enumerate(self.connections):
            if isinstance(connection, RoundedBend):
                if connection.angle is None:
                    raise ValueError("RoundedBend connections require a resolved angle")
                first = self.stations[index].diameter
                second = self.stations[index + 1].diameter
                if not math.isclose(first, second, rel_tol=1.0e-9, abs_tol=1.0e-12):
                    raise ValueError("RoundedBend correlations require equal endpoint diameters")

        if callable(roughness):
            self.roughness = roughness
        else:
            roughness = float(roughness)
            if not np.isfinite(roughness) or roughness < 0.0:
                raise ValueError("roughness must be finite and non-negative")
            self.roughness = roughness

    @property
    def centerline(self) -> np.ndarray:
        """Resolved station positions with shape ``(N, 3)`` [m]."""
        return np.asarray([station.position for station in self.stations], dtype=float)

    @property
    def radius(self) -> np.ndarray:
        """Circular bore radius at each station [m]."""
        return np.asarray([station.diameter / 2.0 for station in self.stations])

    def display_path(self, *, segments: int = 16) -> DisplayPath:
        """Resolve renderable geometry, opening lumped bends into their arcs.

        A :class:`RoundedBend` is one connection between tangent points, so
        :attr:`centerline` crosses its corner as a chord. Here each such bend
        is redrawn as the ``segments``-sided arc that its resolved ``center``
        and ``normal`` describe. A bend left without that geometry keeps its
        chord, because the bend plane is not recoverable from the loss model
        alone. Distributed connections are untouched: :class:`Continuous`
        resolution stays under the user's control through its own ``segments``.

        This changes drawn shape only. Lengths, losses, and every solver
        quantity come from the stations and connections, which are unchanged.
        """
        if isinstance(segments, bool) or not isinstance(segments, int) or segments < 1:
            raise ValueError("segments must be a positive integer")
        centerline = self.centerline
        radii = self.radius
        points = [centerline[0]]
        radius = [radii[0]]
        station = [0.0]
        connection = []
        for index, item in enumerate(self.connections):
            arc = None
            if isinstance(item, RoundedBend) and item.center is not None and item.normal is not None:
                arc = _bend_arc(centerline[index], centerline[index + 1], item, segments)
            if arc is None:
                points.append(centerline[index + 1])
                radius.append(radii[index + 1])
                station.append(index + 1.0)
                connection.append(index)
                continue
            for step, point in enumerate(arc, start=1):
                fraction = step / segments
                points.append(point)
                radius.append(radii[index] + fraction * (radii[index + 1] - radii[index]))
                station.append(index + fraction)
                connection.append(index)
        return DisplayPath(points, radius, station, connection)

    def _roughness(self, s: float) -> float:
        value = self.roughness(s) if callable(self.roughness) else self.roughness
        value = float(value)
        if not np.isfinite(value) or value < 0.0:
            raise ValueError("roughness must evaluate to a finite non-negative value")
        return value

    def _length(self, index: int) -> float:
        connection = self.connections[index]
        if isinstance(connection, RoundedBend):
            return connection.radius * connection.angle
        if isinstance(connection, Continuous) and connection.length is not None:
            return connection.length
        return float(np.linalg.norm(
            np.asarray(self.stations[index + 1].position)
            - np.asarray(self.stations[index].position)
        ))

    def _rounded_K(self, connection, diameter, reynolds, roughness, fd):
        angle_degrees = math.degrees(connection.angle)
        ratio = connection.radius / diameter
        if ratio < 0.5:
            warnings.warn(
                "rounded-bend radius/diameter is below 0.5; fluids correlations are extrapolated",
                PipeCorrelationWarning,
                stacklevel=3,
            )
        if not 0.0 < angle_degrees <= 180.0:
            warnings.warn(
                "rounded-bend angle is outside 0-180 degrees",
                PipeCorrelationWarning,
                stacklevel=3,
            )

        kwargs = dict(Di=diameter, angle=angle_degrees, method=connection.method)
        if connection.method == "Rennels":
            kwargs.update(rc=connection.radius, fd=fd)
        elif connection.method == "Miller":
            kwargs.update(rc=connection.radius, Re=reynolds, roughness=roughness)
        elif connection.method == "Crane":
            kwargs.update(bend_diameters=ratio)
        elif connection.method == "Crane standard":
            pass
        elif connection.method == "Ito":
            kwargs.update(rc=connection.radius, Re=reynolds, roughness=roughness)
        elif connection.method == "Swamee":
            kwargs.update(rc=connection.radius)
        return float(fluids.fittings.bend_rounded(**kwargs))

    def solve(
        self,
        *,
        inlet_state: FluidState,
        mass_flow: float,
        transport,
        friction_model: str = "fluids",
    ) -> PipeResult:
        """March irreversible pressure losses from inlet to outlet.

        Acceleration and heat transfer are intentionally absent in this first
        piping model. Diameter changes alter local velocity, Reynolds number,
        and loss magnitude but do not add a reversible static-pressure term.
        """
        if not isinstance(inlet_state, FluidState):
            raise TypeError("inlet_state must be a thermo.FluidState")
        mass_flow = float(mass_flow)
        if not np.isfinite(mass_flow) or mass_flow <= 0.0:
            raise ValueError("mass_flow must be finite and positive")
        T_in = float(inlet_state.T)
        p_in = float(inlet_state.p)
        if not np.isfinite(T_in) or T_in <= 0.0 or not np.isfinite(p_in) or p_in <= 0.0:
            raise ValueError("inlet state must have finite positive T and p")

        count = len(self.stations)
        lengths = np.asarray([self._length(index) for index in range(count - 1)])
        s = np.concatenate(([0.0], np.cumsum(lengths)))
        p = np.empty(count)
        T = np.full(count, T_in)
        density = np.empty(count)
        velocity = np.empty(count)
        pressure_loss = np.empty(count - 1)
        reynolds_values = np.empty(count - 1)
        darcy_values = np.full(count - 1, np.nan)
        K_values = np.full(count - 1, np.nan)

        p[0] = p_in
        density[0] = _property(transport, "rho", T_in, p_in)
        inlet_area = math.pi * self.stations[0].diameter**2 / 4.0
        velocity[0] = mass_flow / (density[0] * inlet_area)
        states = [_state(T_in, p_in, density[0])]

        for index, connection in enumerate(self.connections):
            D0 = self.stations[index].diameter
            D1 = self.stations[index + 1].diameter
            D = 0.5 * (D0 + D1)
            area = math.pi * D**2 / 4.0
            roughness = self._roughness(0.5 * (s[index] + s[index + 1]))
            p_out = p[index]

            for _ in range(8):
                p_mid = 0.5 * (p[index] + p_out)
                rho = _property(transport, "rho", T_in, p_mid)
                mu = _property(transport, "mu", T_in, p_mid)
                local_velocity = mass_flow / (rho * area)
                reynolds = rho * local_velocity * D / mu
                fd = darcy_friction_factor(reynolds, D, roughness, model=friction_model)

                if isinstance(connection, RoundedBend):
                    K = self._rounded_K(connection, D, reynolds, roughness, fd)
                    loss = K * rho * local_velocity**2 / 2.0
                    factor = math.nan
                else:
                    factor = fd
                    if isinstance(connection, Continuous):
                        radius = connection.radius
                        if radius is None:
                            radius = _estimated_connection_radius(self.stations, index)
                        factor *= curvature_multiplier(
                            reynolds,
                            D,
                            radius,
                            roughness=roughness,
                            model=friction_model,
                        )
                    K = factor * lengths[index] / D
                    loss = K * rho * local_velocity**2 / 2.0

                next_pressure = p[index] - loss
                if next_pressure <= 0.0:
                    raise RuntimeError(
                        f"pipe pressure fell to {next_pressure:g} Pa after connection {index}"
                    )
                if abs(next_pressure - p_out) < 1.0e-6:
                    p_out = next_pressure
                    break
                p_out = next_pressure

            p[index + 1] = p_out
            pressure_loss[index] = loss
            reynolds_values[index] = reynolds
            darcy_values[index] = factor
            K_values[index] = K
            density[index + 1] = _property(transport, "rho", T_in, p_out)
            outlet_area = math.pi * D1**2 / 4.0
            velocity[index + 1] = mass_flow / (density[index + 1] * outlet_area)
            states.append(_state(T_in, p_out, density[index + 1]))

        return PipeResult(
            stations=self.stations,
            connections=self.connections,
            states=tuple(states),
            s=s,
            p=p,
            T=T,
            density=density,
            velocity=velocity,
            pressure_loss=pressure_loss,
            reynolds=reynolds_values,
            darcy_factor=darcy_values,
            loss_coefficient=K_values,
        )
