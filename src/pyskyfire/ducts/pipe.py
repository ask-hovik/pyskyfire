"""Explicit pipe geometry construction and pressure-loss marching.

A :class:`Pipe` is a duct swept along a :class:`Path`: a circular bore,
optionally with a concentric core inside it, such as a centrebody or the inner
line of coaxial tubing. It contains the same :class:`PathStation` and
connection objects that its solver uses. :func:`generate_path_stations` is
only a convenience for turning human-friendly corner vertices into tangent
stations and connections; users may construct the resolved objects directly.
Where the section varies, the pipe adds the stations its loss march needs by
itself; see :class:`Pipe`.

The march, which :class:`~pyskyfire.ducts.AnnularDuct` shares, charges wall
friction, bend losses and the losses of a changing flow area, and carries
swirl as angular momentum. It is adiabatic, so station temperature is
unchanged, and it deliberately does not apply the reversible acceleration term
associated with changing area or compressibility. The pressure it marches is
therefore a stagnation pressure: the inlet's, less the irreversible losses,
which is what the engine network carries between blocks.
``DuctResult.p_static`` subtracts the local dynamic pressure from it; it is
not a complete variable-area static-pressure solution.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import cached_property
import math
from typing import Callable
import warnings

import fluids
import numpy as np

from pyskyfire.thermo.eos import fluid_name
from pyskyfire.thermo.state import FluidState, PHASE_UNKNOWN

from pyskyfire.common.friction import (
    annulus_laminar_factor,
    curvature_multiplier,
    darcy_friction_factor,
)

from .path import (Path, PathStation, Straight, Continuous, RoundedBend, _GEOMETRY_TOL,
                   _arc_angle, _bend_arc, _estimated_connection_radius, _rotate)

__all__ = ["Pipe", "DuctResult", "DisplayPath", "PipeCorrelationWarning"]


class PipeCorrelationWarning(RuntimeWarning):
    """A duct loss correlation is being evaluated outside what it describes."""


@dataclass(frozen=True)
class DuctResult:
    """Station states and per-interval diagnostics from a duct loss march.

    :meth:`Pipe.solve` and :meth:`~pyskyfire.ducts.AnnularDuct.solve` both
    return one. ``p_stagnation`` is the marched pressure: the inlet pressure
    less the irreversible losses up to each station. It includes the swirl's
    dynamic head. The march is adiabatic, so ``T_stagnation`` is the inlet
    temperature at every station, and ``states`` and ``density`` are evaluated
    at the two.

    ``meridional_velocity`` is the through-flow velocity ``c_m`` normal to the
    flow area, ``swirl_velocity`` the circumferential ``c_u`` at the mean
    radius, and ``velocity`` their magnitude, which without swirl is the
    through-flow velocity. ``T_static`` is the stagnation temperature less
    ``velocity**2/(2*cp)``, with ``cp`` taken at the same state as the
    density. Like ``p_static`` it is a first-order correction, not a
    static-state solution.

    Per-interval arrays have one value fewer than the stations.
    ``loss_coefficient`` is referred to the through-flow head
    ``rho*c_m**2/2`` on the interval's mean section. ``darcy_factor`` is the
    factor applied, including any curvature multiplier, and NaN over a lumped
    bend.

    ``stations`` and ``connections`` are a :class:`Pipe`'s resolved path, and
    None for an :class:`~pyskyfire.ducts.AnnularDuct`. ``fluid`` is the
    CoolProp name the transport carried, or None when it named none.
    """

    states: tuple[FluidState, ...]
    s: np.ndarray
    area: np.ndarray
    p_stagnation: np.ndarray
    T_stagnation: np.ndarray
    T_static: np.ndarray
    density: np.ndarray
    meridional_velocity: np.ndarray
    swirl_velocity: np.ndarray
    velocity: np.ndarray
    pressure_loss: np.ndarray
    reynolds: np.ndarray
    darcy_factor: np.ndarray
    loss_coefficient: np.ndarray
    stations: tuple[PathStation, ...] | None = None
    connections: tuple[Straight | Continuous | RoundedBend, ...] | None = None
    fluid: str | None = None

    @property
    def total_pressure_loss(self) -> float:
        """Total irreversible pressure loss [Pa]."""
        return float(self.p_stagnation[0] - self.p_stagnation[-1])

    @property
    def p_static(self) -> np.ndarray:
        """Static pressure at each station [Pa].

        The stagnation pressure less the incompressible dynamic pressure
        ``density*velocity**2/2``, swirl included. At Mach ``M`` that
        understates the compressible dynamic pressure by about ``M**2/4`` of
        itself. With swirl it is the mean across the section: swirl holds the
        static pressure higher at the outer wall than at the inner.
        """
        return self.p_stagnation - 0.5*self.density*self.velocity**2

    @property
    def swirl_angle(self) -> np.ndarray:
        """Flow angle from the through-flow direction, ``atan(c_u/c_m)`` [rad]."""
        return np.arctan2(self.swirl_velocity, self.meridional_velocity)


@dataclass(frozen=True)
class DisplayPath:
    """A refined centerline for drawing, never for pressure-loss marching.

    ``points`` [m] and ``radius`` [m] describe the swept circular bore, and
    ``core_radius`` [m] the core inside it, zero where there is none. The
    path carries its own provenance so results can follow it: ``station`` is
    the fractional :class:`PathStation` index of every point, and
    ``connection`` is the originating connection index of every segment
    between consecutive points.
    """

    points: np.ndarray
    radius: np.ndarray
    station: np.ndarray
    connection: np.ndarray
    core_radius: np.ndarray | None = None

    def __post_init__(self):
        points = np.array(self.points, dtype=float, copy=True)
        radius = np.array(self.radius, dtype=float, copy=True)
        station = np.array(self.station, dtype=float, copy=True)
        connection = np.array(self.connection, dtype=np.int64, copy=True)
        core = (np.zeros_like(radius) if self.core_radius is None
                else np.array(self.core_radius, dtype=float, copy=True))
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("DisplayPath points must be finite with shape (M,3)")
        if radius.shape != (len(points),) or station.shape != (len(points),):
            raise ValueError("DisplayPath radius and station must have one value per point")
        if core.shape != (len(points),) or not np.isfinite(core).all() or np.any(core < 0.0):
            raise ValueError("DisplayPath core_radius must be finite and non-negative, one per point")
        if connection.shape != (max(len(points) - 1, 0),):
            raise ValueError("DisplayPath connection must have one value per segment")
        for array in (points, radius, station, connection, core):
            array.setflags(write=False)
        object.__setattr__(self, "points", points)
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "station", station)
        object.__setattr__(self, "connection", connection)
        object.__setattr__(self, "core_radius", core)

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


def _rounded_bend_K(connection, diameter, reynolds, roughness, fd):
    """Complete loss coefficient of a rounded bend in a circular pipe of ``diameter``."""
    angle_degrees = math.degrees(connection.angle)
    ratio = connection.radius / diameter
    if ratio < 0.5:
        warnings.warn(
            "rounded-bend radius/diameter is below 0.5; fluids correlations are extrapolated",
            PipeCorrelationWarning,
            stacklevel=4,
        )
    if not 0.0 < angle_degrees <= 180.0:
        warnings.warn(
            "rounded-bend angle is outside 0-180 degrees",
            PipeCorrelationWarning,
            stacklevel=4,
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


# ---------------------------------------------------------------------------
# Changing flow area
# ---------------------------------------------------------------------------

# Relative change of equivalent diameter below which an interval counts as a
# constant section.
_CONSTANT_SECTION_TOL = 1.0e-9

# A constant-section stretch shorter than this many hydraulic diameters is
# taken to leave a decelerating or accelerating flow undeveloped, so the area
# change on either side of it is one diffuser or one contraction. A modelling
# choice, not a correlation limit.
_RECOVERY_DIAMETERS = 4.0


def _area_change_runs(change: np.ndarray, lengths: np.ndarray, scale: np.ndarray,
                      sign: int) -> list[tuple[int, int]]:
    """Interval ranges ``(start, stop)`` over which the area keeps changing one way.

    ``sign`` is +1 for widening runs, diffusers, and -1 for narrowing ones,
    contractions. A run ends where the area changes the other way, and at a
    constant-section stretch longer than :data:`_RECOVERY_DIAMETERS` of
    ``scale``, the interval's hydraulic diameter. A shorter one, a
    :class:`RoundedBend` included, stays inside the run. Each run starts and
    ends on an interval that changes the run's way.
    """
    runs = []
    start = last = None
    gap = 0.0
    for index, step in enumerate(sign * change):
        if step > _CONSTANT_SECTION_TOL:
            if start is None:
                start = index
            last, gap = index, 0.0
            continue
        if start is None:
            continue
        if step >= -_CONSTANT_SECTION_TOL:
            gap += lengths[index]
            if gap <= _RECOVERY_DIAMETERS * scale[index]:
                continue
        runs.append((start, last + 1))
        start, gap = None, 0.0
    if start is not None:
        runs.append((start, last + 1))
    return runs


def _area_change_coefficients(area, mid_area, lengths, hydraulic_diameter) -> np.ndarray:
    """Loss coefficient of every interval for its change of flow area, beyond friction.

    The section is read as the circle of the same area, of diameter
    ``sqrt(4*A/pi)``. Each run of widening intervals is a diffuser and loses
    Rennels' conical-diffuser expansion term
    (``fluids.fittings.diffuser_conical``) for its end diameters; each run of
    narrowing ones is a contraction and loses his conical-contraction term
    (``contraction_conical``). Both are taken at the local cone angle,
    averaged over the run with each interval's share of the run's change in
    ``1/A**2`` as weight, which is where the flow slows or speeds up. That is
    Rennels' value for a true cone and for a step alike, and neither the
    station spacing nor a nearly parallel stretch at either end of a run
    changes it. The run's loss is shared out in the same proportion.

    The correlations are called with a zero friction factor, leaving only
    their area-change terms, since the march charges the wall friction.
    Coefficients are on each interval's through-flow head at ``mid_area``.
    """
    area = np.asarray(area, dtype=float)
    mid_area = np.asarray(mid_area, dtype=float)
    diameter = np.sqrt(4.0*area/math.pi)
    change = diameter[1:]/diameter[:-1] - 1.0
    angle = np.degrees(2.0*np.arctan(np.abs(np.diff(diameter))/(2.0*lengths)))
    K = np.zeros(len(lengths))
    for sign, correlation in ((1, fluids.fittings.diffuser_conical),
                              (-1, fluids.fittings.contraction_conical)):
        for start, stop in _area_change_runs(change, lengths, hydraulic_diameter, sign):
            A = area[start:stop + 1]
            weight = sign*(A[:-1]**-2.0 - A[1:]**-2.0)
            active = weight > 0.0
            weight = np.where(active, weight, 0.0)/weight[active].sum()
            local = np.array([
                correlation(diameter[start], diameter[stop], angle=a, fd=0.0) if w > 0.0 else 0.0
                for a, w in zip(angle[start:stop], weight)
            ])
            # Both correlations refer their coefficient to the narrow end's head.
            narrow = min(area[start], area[stop])
            K[start:stop] += float(np.dot(weight, local))*weight*(mid_area[start:stop]/narrow)**2
    return K


# ---------------------------------------------------------------------------
# The march
# ---------------------------------------------------------------------------

def _roughness(roughness):
    if callable(roughness):
        return roughness
    roughness = float(roughness)
    if not np.isfinite(roughness) or roughness < 0.0:
        raise ValueError("roughness must be finite and non-negative")
    return roughness


def _roughness_at(roughness, s):
    value = float(roughness(s) if callable(roughness) else roughness)
    if not np.isfinite(value) or value < 0.0:
        raise ValueError("roughness must evaluate to a finite non-negative value")
    return value


def _march(*, lengths, area, mid_area, hydraulic_diameter, r_mean, mid_r_mean, laminar_factor,
           curvature_radius, bends, bend_diameter, area_change, roughness, inlet_state,
           mass_flow, inlet_swirl, transport, friction_model, stations=None,
           connections=None) -> DuctResult:
    """The march every duct shares, over ``len(lengths)`` intervals.

    Station arrays have one more value than the interval ones. ``bends`` holds
    a :class:`RoundedBend` for an interval charged as one lumped bend, and
    None for distributed friction on ``curvature_radius``. ``area_change`` is
    each interval's area-change loss coefficient, from
    :func:`_area_change_coefficients`.
    """
    if not isinstance(inlet_state, FluidState):
        raise TypeError("inlet_state must be a thermo.FluidState")
    mass_flow = float(mass_flow)
    if not np.isfinite(mass_flow) or mass_flow <= 0.0:
        raise ValueError("mass_flow must be finite and positive")
    inlet_swirl = float(inlet_swirl)
    if not np.isfinite(inlet_swirl):
        raise ValueError("inlet_swirl must be finite")
    T_in = float(inlet_state.T)
    p_in = float(inlet_state.p)
    if not np.isfinite(T_in) or T_in <= 0.0 or not np.isfinite(p_in) or p_in <= 0.0:
        raise ValueError("inlet state must have finite positive T and p")

    count = len(lengths) + 1
    s = np.r_[0.0, np.cumsum(lengths)]
    p = np.empty(count)
    density = np.empty(count)
    c_m = np.empty(count)
    angular_momentum = np.empty(count)
    pressure_loss = np.empty(count - 1)
    reynolds_values = np.empty(count - 1)
    darcy_values = np.full(count - 1, np.nan)
    K_values = np.empty(count - 1)

    p[0] = p_in
    density[0] = _property(transport, "rho", T_in, p_in)
    c_m[0] = mass_flow/(density[0]*area[0])
    angular_momentum[0] = r_mean[0]*inlet_swirl
    states = [_state(T_in, p_in, density[0])]

    for i in range(count - 1):
        D = hydraulic_diameter[i]
        L = lengths[i]
        eps = _roughness_at(roughness, 0.5*(s[i] + s[i + 1]))
        p_out, rcu_out = p[i], angular_momentum[i]
        for _ in range(12):
            p_mid = 0.5*(p[i] + p_out)
            rho = _property(transport, "rho", T_in, p_mid)
            mu = _property(transport, "mu", T_in, p_mid)
            cm = mass_flow/(rho*mid_area[i])
            cu = 0.5*(angular_momentum[i] + rcu_out)/mid_r_mean[i]
            c = math.hypot(cm, cu)
            # The wall shear lies along the total velocity, so a swirling
            # flow rubs over a helix longer than the duct by c/c_m.
            helix = c/cm
            reynolds = rho*c*D/mu
            fd = darcy_friction_factor(reynolds, D, eps, model=friction_model,
                                       laminar_factor=laminar_factor[i])
            bend = bends[i]
            if bend is not None:
                Db = bend_diameter[i]
                reynolds_b = rho*cm*Db/mu
                fd_b = darcy_friction_factor(reynolds_b, Db, eps, model=friction_model)
                K = _rounded_bend_K(bend, Db, reynolds_b, eps, fd_b)
                factor = math.nan
                decay = fd*L/(2.0*D)*helix
            else:
                factor = fd
                if np.isfinite(curvature_radius[i]):
                    factor *= curvature_multiplier(reynolds, D, curvature_radius[i],
                                                   roughness=eps, model=friction_model)
                K = factor*L/D*helix**3 + area_change[i]
                decay = factor*L/(2.0*D)*helix
            loss = K*rho*cm**2/2.0
            next_rcu = angular_momentum[i]*math.exp(-decay)
            next_pressure = p[i] - loss
            if next_pressure <= 0.0:
                raise RuntimeError(f"duct pressure fell to {next_pressure:g} Pa after interval {i}")
            converged = abs(next_pressure - p_out) < 1.0e-6 and abs(next_rcu - rcu_out) < 1.0e-12
            p_out, rcu_out = next_pressure, next_rcu
            if converged:
                break

        p[i + 1] = p_out
        angular_momentum[i + 1] = rcu_out
        pressure_loss[i] = loss
        reynolds_values[i] = reynolds
        darcy_values[i] = factor
        K_values[i] = K
        density[i + 1] = _property(transport, "rho", T_in, p_out)
        c_m[i + 1] = mass_flow/(density[i + 1]*area[i + 1])
        states.append(_state(T_in, p_out, density[i + 1]))

    c_u = angular_momentum/r_mean
    velocity = np.hypot(c_m, c_u)
    T = np.full(count, T_in)
    cp = np.array([_property(transport, "cp", T_in, value) for value in p])
    return DuctResult(
        states=tuple(states), s=s, area=np.asarray(area, float).copy(),
        p_stagnation=p, T_stagnation=T, T_static=T - velocity**2/(2.0*cp),
        density=density, meridional_velocity=c_m, swirl_velocity=c_u, velocity=velocity,
        pressure_loss=pressure_loss, reynolds=reynolds_values, darcy_factor=darcy_values,
        loss_coefficient=K_values, stations=stations, connections=connections,
        fluid=fluid_name(transport),
    )


# ---------------------------------------------------------------------------
# Station refinement
# ---------------------------------------------------------------------------

# Diameter samples taken along each connection when deciding how finely to
# split it. A section feature shorter than 1/_DIAMETER_SAMPLES of its
# connection can be missed.
_DIAMETER_SAMPLES = 64


def _diameter_law(diameter, path: Path, name: str) -> tuple[Callable[[np.ndarray], np.ndarray], bool]:
    """Return a diameter as a function of arc length, and whether it is constant.

    One value per station is read as a law that is linear in arc length between
    stations, which is also what the march assumes within one connection.
    """
    if callable(diameter):
        return (lambda s: np.asarray([float(diameter(float(x))) for x in np.atleast_1d(s)]),
                False)
    values = np.asarray(diameter, dtype=float)
    if values.ndim == 0 or values.size == 1:
        value = float(values.ravel()[0])
        return (lambda s: np.full(np.shape(np.atleast_1d(s)), value)), True
    if values.shape != (len(path.stations),):
        raise ValueError(f"{name} must be a scalar, a function of s, or one value per path station")
    s = path.s
    return (lambda x: np.interp(x, s, values)), bool(np.all(values == values[0]))


def _split_fractions(bore, core, s0: float, s1: float, max_change: float) -> np.ndarray:
    """Interior split points, as fractions of one connection, from its diameter laws.

    The march evaluates a connection at its mean section, through its area
    and hydraulic diameter. The splits equidistribute the larger of the
    variations of ``ln`` of the two, the area taken as the diameter
    ``sqrt(D**2 - d**2)`` of the circle of the same area, so neither varies by
    more than ``max_change`` across a piece and stations crowd where the
    section changes fastest. Both stay finite as a core shrinks to nothing.
    Variation is measured inside the connection, not only at its ends, so a
    bulge between equal end diameters is also resolved.
    """
    x = np.linspace(s0, s1, _DIAMETER_SAMPLES + 1)
    outer, inner = bore(x), core(x)
    if (not np.isfinite(outer).all() or not np.isfinite(inner).all() or np.any(outer <= 0)
            or np.any(inner < 0) or np.any(inner >= outer)):
        raise ValueError("diameter must be finite and positive along the whole pipe, "
                         "and any core smaller than it")
    step = np.maximum(np.abs(np.diff(0.5*np.log(outer**2 - inner**2))),
                      np.abs(np.diff(np.log(outer - inner))))
    variation = np.concatenate(([0.0], np.cumsum(step)))
    count = max(1, math.ceil(variation[-1] / math.log1p(max_change) - 1.0e-9))
    if count == 1:
        return np.empty(0)
    targets = variation[-1] * np.arange(1, count) / count
    return (np.interp(targets, variation, x) - s0) / (s1 - s0)


def _split_connection(path: Path, index: int, fractions: np.ndarray):
    """Interior stations and the connections joining them, for one split connection."""
    connection = path.connections[index]
    start, end = path.centerline[index], path.centerline[index + 1]
    length = path.lengths[index]
    bounds = np.concatenate(([0.0], fractions, [1.0]))
    if isinstance(connection, Continuous) and connection.center is not None:
        center, normal = np.asarray(connection.center), np.asarray(connection.normal)
        angle = _arc_angle(connection)
        points = [center + _rotate(start - center, normal, f * angle) for f in fractions]
    else:
        points = [start + f * (end - start) for f in fractions]
    if isinstance(connection, Continuous):
        radius = connection.radius
        if radius is None:
            # Freeze the estimate from the original stations: the new ones lie on
            # the chord, where their own triples would see no curvature at all.
            estimate = _estimated_connection_radius(path.stations, index)
            radius = estimate if np.isfinite(estimate) else None
        pieces = [replace(connection, radius=radius, length=length * (b - a))
                  for a, b in zip(bounds[:-1], bounds[1:])]
    else:
        pieces = [connection] * (len(bounds) - 1)
    return [PathStation(tuple(point)) for point in points], pieces


def _resolve_path(path: Path, bore, core, max_change: float) -> Path:
    """Split distributed connections until the section varies little across each."""
    s = path.s
    stations, connections = [path.stations[0]], []
    for index, connection in enumerate(path.connections):
        fractions = np.empty(0)
        if not isinstance(connection, RoundedBend):
            fractions = _split_fractions(bore, core, s[index], s[index + 1], max_change)
        if fractions.size:
            interior, pieces = _split_connection(path, index, fractions)
            stations.extend(interior)
            connections.extend(pieces)
        else:
            connections.append(connection)
        stations.append(path.stations[index + 1])
    if len(stations) == len(path.stations):
        return path
    return Path(stations, connections)


_BEND_LENGTH_SCALES = ("hydraulic", "outer")


class Pipe:
    """A circular duct swept along a path, with an optional core, and its loss solver.

    The section normal to the path is a circular bore of ``diameter``, with a
    concentric core of ``core_diameter`` inside it where that is not zero. A
    plain pipe has no core. A core may run the whole length, as in coaxial
    lines, or only part of it: it may start from a point, as a turbine's nose
    does, and end, bluntly or tapering.

    The march evaluates each connection at one section, the mean of its two
    ends. Where the section varies, the pipe therefore splits
    :class:`Straight` and :class:`Continuous` connections itself until its
    area and hydraulic diameter each change by at most
    ``max_diameter_change`` (relative, the area as the diameter of the circle
    of the same area) across each piece, placing the extra stations where the
    section changes fastest. A constant section is left exactly as supplied,
    and so is every :class:`RoundedBend`. :attr:`path`, :attr:`stations` and
    every result refer to the resolved stations, which include all of the
    supplied ones.

    Parameters
    ----------
    path : Path or sequence of PathStation
        Resolved stations, with ``connections`` if a plain sequence is given.
    diameter : float, sequence or callable
        Bore diameter [m]: a constant, one value per supplied station (linear
        in arc length between them), or a function of arc length ``s`` [m].
    core_diameter : float, sequence or callable, optional
        Core diameter [m], given the same ways, smaller than the bore
        everywhere. Zero, the default, where there is no core.
    roughness : float or callable, optional
        Absolute roughness of every wall [m], or a function of ``s``.
    max_diameter_change : float, optional
        Largest relative change of the section across one connection. The
        default, 2 %, keeps the error of the single-section evaluation below
        0.05 %. ``math.inf`` keeps the supplied stations exactly.
    """

    def __init__(self, path, connections=None, *, diameter, core_diameter=0.0,
                 roughness: float | Callable[[float], float] = 0.0,
                 max_diameter_change: float = 0.02):
        if isinstance(path, Path) and connections is not None:
            raise ValueError("connections are already supplied by Path")
        path = path if isinstance(path, Path) else Path(path, connections)
        max_diameter_change = float(max_diameter_change)
        if np.isnan(max_diameter_change) or max_diameter_change <= 0.0:
            raise ValueError("max_diameter_change must be positive")
        self.max_diameter_change = max_diameter_change
        bore, bore_constant = _diameter_law(diameter, path, "diameter")
        core, core_constant = _diameter_law(core_diameter, path, "core_diameter")
        if not (bore_constant and core_constant) and np.isfinite(max_diameter_change):
            path = _resolve_path(path, bore, core, max_diameter_change)
        self.path = path
        self.stations = self.path.stations
        self.connections = self.path.connections
        if len(self.stations) < 2:
            raise ValueError("Pipe requires at least two stations")
        if not all(isinstance(station, PathStation) for station in self.stations):
            raise TypeError("stations must contain only PathStation instances")
        if len(self.connections) != len(self.stations) - 1:
            raise ValueError("Pipe requires exactly one connection between each station pair")
        if not all(isinstance(item, (Straight, Continuous, RoundedBend)) for item in self.connections):
            raise TypeError("unsupported pipe connection type")

        self.diameter = np.array(bore(self.path.s), dtype=float)
        self.core_diameter = np.array(core(self.path.s), dtype=float)
        if not np.isfinite(self.diameter).all() or np.any(self.diameter <= 0):
            raise ValueError("diameter must be finite and positive at every station")
        if not np.isfinite(self.core_diameter).all() or np.any(self.core_diameter < 0):
            raise ValueError("core_diameter must be finite and non-negative at every station")
        if np.any(self.core_diameter >= self.diameter):
            raise ValueError("core_diameter must be smaller than diameter at every station")
        self.diameter.setflags(write=False)
        self.core_diameter.setflags(write=False)

        lengths = np.linalg.norm(np.diff(self.centerline, axis=0), axis=1)
        if np.any(lengths <= _GEOMETRY_TOL):
            raise ValueError("Pipe stations must not contain consecutive duplicate positions")

        curvature = self.path.curvature_radius
        for index, connection in enumerate(self.connections):
            if isinstance(connection, RoundedBend):
                if connection.angle is None:
                    raise ValueError("RoundedBend connections require a resolved angle")
                for wall in (self.diameter, self.core_diameter):
                    if not math.isclose(wall[index], wall[index + 1], rel_tol=1.0e-9, abs_tol=1.0e-12):
                        raise ValueError("RoundedBend correlations require equal endpoint diameters")
            if isinstance(connection, (Continuous, RoundedBend)):
                if curvature[index] <= 0.5*max(self.diameter[index], self.diameter[index + 1]):
                    raise ValueError("a bend's centreline radius must exceed the bore radius, "
                                     "or its inner side folds")

        self.roughness = _roughness(roughness)

    # -- geometry ------------------------------------------------------------
    @property
    def centerline(self) -> np.ndarray:
        """Resolved station positions with shape ``(N, 3)`` [m]."""
        return np.asarray([station.position for station in self.stations], dtype=float)

    @property
    def s(self) -> np.ndarray:
        """Arc length at each station [m]."""
        return self.path.s

    @property
    def radius(self) -> np.ndarray:
        """Bore radius at each station [m]."""
        return self.diameter / 2.0

    @property
    def core_radius(self) -> np.ndarray:
        """Core radius at each station [m], zero where there is no core."""
        return self.core_diameter / 2.0

    @property
    def has_core(self) -> bool:
        """Whether a core runs along any of the pipe."""
        return bool(np.any(self.core_diameter > 0.0))

    @property
    def area(self) -> np.ndarray:
        """Flow area at each station [m²]."""
        return math.pi/4.0*(self.diameter**2 - self.core_diameter**2)

    @property
    def hydraulic_diameter(self) -> np.ndarray:
        """``4A/P`` at each station, ``D - d`` [m]: the bore where there is no core."""
        return self.diameter - self.core_diameter

    @property
    def r_mean(self) -> np.ndarray:
        """Mean radius of the section, ``(D + d)/4`` [m], where swirl is measured."""
        return 0.25*(self.diameter + self.core_diameter)

    @property
    def cone_angle(self) -> np.ndarray:
        """Full cone angle of each connection [rad], positive where the area grows.

        The angle of the cone whose area grows as the section's does, from
        the diameter ``sqrt(4A/pi)`` at the connection's ends over its length
        along the centerline. For a pipe without a core it is the bore's own
        wall angle. It is local wherever :class:`Pipe` has split a varying
        section; with ``max_diameter_change=math.inf`` it only sees the
        supplied stations.
        """
        diameter = np.sqrt(4.0*self.area/math.pi)
        return 2.0 * np.arctan(np.diff(diameter) / (2.0 * self.path.lengths))

    @cached_property
    def area_change_loss_coefficient(self) -> np.ndarray:
        """Loss coefficient of each connection for its change of area [-].

        Where the flow area keeps growing, the flow loses to its slowing as a
        diffuser does, and where it keeps shrinking, to its speeding up as a
        contraction does. :meth:`solve` adds this to each connection's wall
        friction. A run of connections changing area the same way is one
        diffuser or contraction: it ends where the area changes the other way,
        or at a constant-section stretch longer than four hydraulic
        diameters. Each run loses Rennels' conical-diffuser or
        conical-contraction term for its end sections, at the cone angle
        averaged over the run, weighted by where the flow slows or speeds up.
        That is his value for a true cone and for a step alike, and it is
        unchanged by how finely the pipe is split or by a nearly parallel
        stretch at either end. The run's loss is shared out in the same
        proportion. Coefficients are on each connection's through-flow head
        at its mean section; zero where the area is constant.
        """
        Do, Di = self.diameter, self.core_diameter
        mid_Do, mid_Di = 0.5*(Do[:-1] + Do[1:]), 0.5*(Di[:-1] + Di[1:])
        K = _area_change_coefficients(self.area, math.pi/4.0*(mid_Do**2 - mid_Di**2),
                                      self.path.lengths, mid_Do - mid_Di)
        K.setflags(write=False)
        return K

    def display_path(self, *, segments: int = 16) -> DisplayPath:
        """Resolve renderable geometry, opening arc connections into their arcs.

        A :class:`RoundedBend` is one connection between tangent points, and a
        generated :class:`Continuous` arc may be one connection too, so
        :attr:`centerline` crosses such a corner as a chord. Here every arc
        connection with a resolved ``center`` and ``normal`` is redrawn as an
        arc of ``segments`` sides per quarter turn, at least one. A connection
        without that geometry keeps its chord, because its plane is not
        recoverable from the loss model alone.

        This changes drawn shape only. Lengths, losses, and every solver
        quantity come from the stations and connections, which are unchanged.
        """
        if isinstance(segments, bool) or not isinstance(segments, int) or segments < 1:
            raise ValueError("segments must be a positive integer")
        centerline = self.centerline
        radii, cores = self.radius, self.core_radius
        points = [centerline[0]]
        radius = [radii[0]]
        core = [cores[0]]
        station = [0.0]
        connection = []
        for index, item in enumerate(self.connections):
            arc = None
            if (isinstance(item, (RoundedBend, Continuous))
                    and item.center is not None and item.normal is not None):
                count = max(1, math.ceil(segments * _arc_angle(item) / (math.pi / 2) - 1.0e-9))
                arc = _bend_arc(centerline[index], centerline[index + 1], item, count)
            if arc is None:
                arc = [centerline[index + 1]]
            for step, point in enumerate(arc, start=1):
                fraction = step / len(arc)
                points.append(point)
                radius.append(radii[index] + fraction * (radii[index + 1] - radii[index]))
                core.append(cores[index] + fraction * (cores[index + 1] - cores[index]))
                station.append(index + fraction)
                connection.append(index)
        return DisplayPath(points, radius, station, connection, core)

    # -- solver --------------------------------------------------------------
    def solve(
        self,
        *,
        inlet_state: FluidState,
        mass_flow: float,
        transport,
        inlet_swirl: float = 0.0,
        friction_model: str = "fluids",
        bend_length_scale: str = "hydraulic",
    ) -> DuctResult:
        """March irreversible pressure losses from inlet to outlet.

        :class:`Straight` and :class:`Continuous` connections charge wall
        friction on the hydraulic diameter, the latter times the curvature
        multiplier, along the helix a swirling flow follows, plus the loss of
        any change of area (:attr:`area_change_loss_coefficient`). A
        :class:`RoundedBend` charges its complete loss coefficient. Every loss
        is on the through-flow head ``rho*c_m**2/2``. Acceleration and heat
        transfer are intentionally absent: a change of area alters the
        velocities, and so the static pressure, but adds no reversible term
        to the marched stagnation pressure.

        Parameters
        ----------
        inlet_state : FluidState
            ``p`` and ``T`` are taken as the inlet stagnation pressure and
            temperature.
        mass_flow : float
            [kg/s]
        transport
            Provides ``rho``, ``mu`` and ``cp`` as ``get_<name>(T, p)`` or
            ``<name>(T, p)``.
        inlet_swirl : float, optional
            Circumferential velocity about the centreline at the inlet mean
            radius :attr:`r_mean` [m/s], positive in either sense as long as
            it is used consistently. It is carried as angular momentum,
            unchanged by a change of radius, and worn down by wall friction.
            Default 0.
        friction_model : {'fluids', 'native'}, optional
            Backend for the friction factor and the curvature multiplier.
        bend_length_scale : {'hydraulic', 'outer'}, optional
            The diameter a :class:`RoundedBend` around a core is given to the
            circular-pipe correlation: the hydraulic diameter ``D - d``, or
            the bore. No annular-bend correlation exists here to choose
            between them; solving with both brackets the bend loss. Without a
            core both are the bore.

        Warns
        -----
        PipeCorrelationWarning
            When a :class:`RoundedBend` turns round a core, since its loss
            comes from a circular-pipe correlation, and when swirl enters a
            :class:`RoundedBend`, since the bend correlations are for flow
            without swirl.
        """
        if bend_length_scale not in _BEND_LENGTH_SCALES:
            raise ValueError(f"bend_length_scale must be one of {_BEND_LENGTH_SCALES}")
        bends = [c if isinstance(c, RoundedBend) else None for c in self.connections]
        if any(b is not None for b in bends):
            if self.has_core:
                warnings.warn(
                    "a rounded bend round a core uses a circular-pipe correlation on one length "
                    "scale of the annulus; no annular-bend correlation is applied",
                    PipeCorrelationWarning, stacklevel=2,
                )
            if inlet_swirl:
                warnings.warn(
                    "swirl enters a rounded bend; the bend correlation is for flow without swirl",
                    PipeCorrelationWarning, stacklevel=2,
                )
        Do, Di = self.diameter, self.core_diameter
        mid_Do, mid_Di = 0.5*(Do[:-1] + Do[1:]), 0.5*(Di[:-1] + Di[1:])
        mid_D = mid_Do - mid_Di
        laminar = [annulus_laminar_factor(ratio) for ratio in mid_Di/mid_Do]
        return _march(
            lengths=self.path.lengths, area=self.area,
            mid_area=math.pi/4.0*(mid_Do**2 - mid_Di**2), hydraulic_diameter=mid_D,
            r_mean=self.r_mean, mid_r_mean=0.25*(mid_Di + mid_Do), laminar_factor=laminar,
            curvature_radius=self.path.curvature_radius, bends=bends,
            bend_diameter=mid_D if bend_length_scale == "hydraulic" else mid_Do,
            area_change=self.area_change_loss_coefficient, roughness=self.roughness,
            inlet_state=inlet_state, mass_flow=mass_flow, inlet_swirl=inlet_swirl,
            transport=transport, friction_model=friction_model,
            stations=self.stations, connections=self.connections,
        )
