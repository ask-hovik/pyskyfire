"""Explicit pipe geometry construction and pressure-loss marching.

A :class:`Pipe` contains the same :class:`PathStation` and connection objects
that its solver uses. :func:`generate_path_stations` is only a convenience for
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

from pyskyfire.common.friction import curvature_multiplier, darcy_friction_factor

from .path import Path, PathStation, Straight, Continuous, RoundedBend, _GEOMETRY_TOL, _bend_arc, _estimated_connection_radius

__all__ = ["Pipe", "PipeResult", "DisplayPath", "PipeCorrelationWarning"]

class PipeCorrelationWarning(RuntimeWarning):
    """A pipe correlation is being evaluated outside its documented range."""

@dataclass(frozen=True)
class PipeResult:
    """Station states and per-connection diagnostics from a pipe loss march."""

    stations: tuple[PathStation, ...]
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
    the fractional :class:`PathStation` index of every point, and
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
        """Number of source :class:`PathStation` objects the path came from."""
        return self.connection_count + 1 if len(self.connection) else len(self.points)


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

    def __init__(self, path, connections=None, *, diameter, roughness: float | Callable[[float], float] = 0.0):
        self.path = path if isinstance(path, Path) else Path(path, connections)
        if isinstance(path, Path) and connections is not None:
            raise ValueError("connections are already supplied by Path")
        self.stations = self.path.stations
        self.connections = self.path.connections
        if callable(diameter):
            values = [diameter(float(s)) for s in self.path.s]
        else:
            values = diameter
        self.diameter = np.broadcast_to(np.asarray(values, float), (len(self.stations),)).copy()
        if not np.isfinite(self.diameter).all() or np.any(self.diameter <= 0):
            raise ValueError("diameter must be finite and positive at every station")
        self.diameter.setflags(write=False)
        if len(self.stations) < 2:
            raise ValueError("Pipe requires at least two stations")
        if not all(isinstance(station, PathStation) for station in self.stations):
            raise TypeError("stations must contain only PathStation instances")
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
                first = self.diameter[index]
                second = self.diameter[index + 1]
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
        return self.diameter / 2.0

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
        inlet_area = math.pi * self.diameter[0]**2 / 4.0
        velocity[0] = mass_flow / (density[0] * inlet_area)
        states = [_state(T_in, p_in, density[0])]

        for index, connection in enumerate(self.connections):
            D0 = self.diameter[index]
            D1 = self.diameter[index + 1]
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
