"""The duct carrying collected flow from a volute to its pipe flange.

A volute delivers into piping whose diameter is usually fixed by a standard,
and the direction it delivers in varies with how the machine is packaged: some
outlets run straight on from the full-collection plane, some curve round until
they point radially outward, and real hardware sits anywhere between. The exit
cone is the transition that reconciles the passage with that flange.

Direction is given by two angles.

**First angle** turns within the plane the volute sits in, measured from
*tangential*. Tangential is fixed at the volute outlet -- the direction the
flow is already travelling when it leaves the full-collection plane -- because
by the time the duct reaches its flange it has moved round the volute, and
there is no longer room for it to be tangential to the passage there. Radial
is not fixed that way: it always means straight out from the machine axis, so
it is evaluated where the flange actually ends up. Having travelled round to
azimuth ``phi``, that radial direction is only ``pi/2 - phi`` away from
tangential, and it is that angle -- not a right angle -- that a first angle may
not exceed. Beyond it the duct would have turned past radial and be heading
back inward.

**Second angle** lifts the outlet out of the volute's plane, positive toward
``+z``. It applies on top of the first angle, so a radial outlet with a second
angle still leaves in the plane through the machine axis and the flange.

The two together fix the flange's normal, and with the arc length they fix the
whole duct. A purely tangential outlet has nothing to turn through, so it runs
straight.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.optimize import brentq

from .path import _rotate

__all__ = [
    "ExitCone",
    "ResolvedExitCone",
]

#: Directions the first angle may name instead of a number.
NAMED_ANGLES = ("tangential", "radial")

#: Below this turn the duct is treated as straight rather than as a huge arc.
_STRAIGHT_TOLERANCE = 1.0e-9


@dataclass(frozen=True)
class ExitCone:
    """How a volute's outlet duct leaves, and what it grows into.

    Parameters
    ----------
    final_diameter : float
        Diameter of the flange the duct delivers to [m], normally set by a
        piping standard rather than by the volute.
    first_angle : float or {'tangential', 'radial'}, optional
        Turn within the volute's plane, measured from tangential at the volute
        outlet [rad]. ``'tangential'`` is zero and runs straight on;
        ``'radial'`` turns until the flange points straight out from the
        machine axis. A number must lie between the two, since a larger turn
        would point the duct back inward.
    second_angle : float, optional
        Elevation out of the volute's plane [rad], positive upward. Applies on
        top of the first angle.
    arc_length_parameter : float, optional
        Centreline length of the duct, as a multiple of the volute outlet
        diameter. This is the distance the cross-section has to grow in, so it
        must be at least one diameter.

    Notes
    -----
    Exit-cone geometry does not feed back into the loss analysis. Aungier's
    exit cone loss is already a lumped allowance for this whole transition,
    including the turn out of the volute's log-spiral path, so charging a
    separate bend loss on top of it would count the same loss twice. Only
    :attr:`final_diameter` reaches the analysis, as station 7.
    """

    final_diameter: float
    first_angle: float | str = "tangential"
    second_angle: float = 0.0
    arc_length_parameter: float = 2.0

    def __post_init__(self):
        diameter = float(self.final_diameter)
        if not math.isfinite(diameter) or diameter <= 0.0:
            raise ValueError("final_diameter must be finite and positive")
        object.__setattr__(self, "final_diameter", diameter)

        first = self.first_angle
        if isinstance(first, str):
            if first not in NAMED_ANGLES:
                raise ValueError(
                    f"first_angle must be a number or one of {NAMED_ANGLES}")
        else:
            first = float(first)
            if not math.isfinite(first):
                raise ValueError("first_angle must be finite")
            if first < 0.0:
                raise ValueError(
                    "first_angle must be at least zero; a negative turn would "
                    "point the duct back against the collected flow")
            object.__setattr__(self, "first_angle", first)

        second = float(self.second_angle)
        if not math.isfinite(second) or abs(second) > math.pi / 2.0:
            raise ValueError("second_angle must be finite and within +/- pi/2")
        object.__setattr__(self, "second_angle", second)

        parameter = float(self.arc_length_parameter)
        if not math.isfinite(parameter) or parameter < 1.0:
            raise ValueError(
                "arc_length_parameter must be at least 1: the cross-section "
                "needs at least one outlet diameter to grow in")
        object.__setattr__(self, "arc_length_parameter", parameter)


@dataclass(frozen=True)
class ResolvedExitCone:
    """An exit cone placed in machine coordinates.

    Attributes
    ----------
    start, end : np.ndarray
        Centres of the volute outlet and of the flange [m, shape (3,)].
    start_direction, end_direction : np.ndarray
        Unit normals there [shape (3,)]. ``start_direction`` is tangential at
        the volute outlet.
    first_angle, second_angle : float
        The resolved turn within the volute's plane and the elevation [rad].
    radial_angle : float
        The turn that would have put the flange radial [rad]. Equal to
        :attr:`first_angle` when ``'radial'`` was asked for, and the upper
        bound on it otherwise.
    end_azimuth : float
        Azimuth the flange ends up at, measured from the volute outlet [rad].
    arc_length, turn_angle, bend_radius : float
        Centreline length [m], total turn [rad], and centreline bend radius
        [m], which is infinite for a straight duct.
    outlet_diameter, final_diameter : float
        Volute outlet bore diameter and flange diameter [m].
    """

    start: np.ndarray
    end: np.ndarray
    start_direction: np.ndarray
    end_direction: np.ndarray
    first_angle: float
    second_angle: float
    radial_angle: float
    end_azimuth: float
    arc_length: float
    turn_angle: float
    bend_radius: float
    outlet_diameter: float
    final_diameter: float

    @property
    def is_straight(self) -> bool:
        return not math.isfinite(self.bend_radius)

    def frames(self, n_stations=17, *, start_frame=None):
        """Centreline points and a non-twisting cross-section frame.

        Parameters
        ----------
        start_frame : tuple of np.ndarray, optional
            Override for the frame at the volute outlet, as ``(u0, v0)`` unit
            vectors perpendicular to :attr:`start_direction`. Defaults to
            :meth:`_start_frame`, the untilted ``(radial, axial)`` pair. A
            volute with a tilted section passes its own tilted section-plane
            basis here instead, so the duct's first ring is bit-identical to
            the volute's own last one rather than merely tangent to the same
            curve at different points along it.

        Returns
        -------
        points, u, v : np.ndarray
            Centreline [m, shape (N,3)] and two unit vectors spanning each
            station's cross-section plane [shape (N,3)]. The arc is planar, so
            rotating the starting frame about the arc's own axis carries it
            along without twist, whatever that starting frame is.
        """
        n_stations = int(n_stations)
        if n_stations < 2:
            raise ValueError("n_stations must be at least two")
        fraction = np.linspace(0.0, 1.0, n_stations)
        u0, v0 = self._start_frame() if start_frame is None else start_frame
        if self.is_straight:
            points = self.start + np.outer(fraction * self.arc_length,
                                           self.start_direction)
            return points, np.tile(u0, (n_stations, 1)), np.tile(v0, (n_stations, 1))
        axis = self._axis()
        centre = self.start + self.bend_radius * np.cross(axis, self.start_direction)
        spoke = self.start - centre
        points = np.array([centre + _rotate(spoke, axis, angle)
                           for angle in fraction * self.turn_angle])
        u = np.array([_rotate(u0, axis, angle) for angle in fraction * self.turn_angle])
        v = np.array([_rotate(v0, axis, angle) for angle in fraction * self.turn_angle])
        return points, u, v

    def _axis(self):
        axis = np.cross(self.start_direction, self.end_direction)
        return axis / np.linalg.norm(axis)

    def _start_frame(self):
        """The volute's own section frame, which spans the outlet plane.

        A volute section lies in the plane through the machine axis, so its
        normal is the tangential direction the duct starts along. Reusing it
        means the duct's first ring is exactly the volute's last one.
        """
        axial = np.array([0.0, 0.0, 1.0])
        radial = np.cross(self.start_direction, axial)
        radial /= np.linalg.norm(radial)
        return radial, axial


def _direction(first, second, tangential, radial, axial):
    """Unit flange normal for a turn within the plane plus an elevation."""
    horizontal = math.cos(first)*tangential + math.sin(first)*radial
    return math.cos(second)*horizontal + math.sin(second)*axial


def _place(start, start_direction, end_direction, arc_length):
    """Return ``(end, turn, bend_radius)`` for an arc of the given length."""
    turn = math.acos(float(np.clip(np.dot(start_direction, end_direction), -1.0, 1.0)))
    if turn < _STRAIGHT_TOLERANCE:
        return start + arc_length*start_direction, 0.0, math.inf
    axis = np.cross(start_direction, end_direction)
    axis = axis / np.linalg.norm(axis)
    bend_radius = arc_length / turn
    centre = start + bend_radius*np.cross(axis, start_direction)
    return centre + _rotate(start - centre, axis, turn), turn, bend_radius


def resolve_exit_cone(cone, *, start, outlet_diameter, max_section_radius):
    """Place an :class:`ExitCone` on a volute's full-collection plane.

    Parameters
    ----------
    cone : ExitCone
        The requested outlet.
    start : array_like
        Centre of the volute outlet in machine coordinates [m, shape (3,)].
    outlet_diameter : float
        Volute outlet bore diameter [m], which scales the arc length.
    max_section_radius : float
        Largest distance from the outlet centre to the volute section boundary
        [m]. Used to reject a bend tight enough to fold the duct into itself.

    Returns
    -------
    ResolvedExitCone
    """
    start = np.asarray(start, dtype=float)
    if start.shape != (3,) or not np.isfinite(start).all():
        raise ValueError("start must be a finite (3,) point")
    azimuth = math.atan2(start[1], start[0])
    tangential = np.array([-math.sin(azimuth), math.cos(azimuth), 0.0])
    radial = np.array([math.cos(azimuth), math.sin(azimuth), 0.0])
    axial = np.array([0.0, 0.0, 1.0])
    arc_length = cone.arc_length_parameter * outlet_diameter
    second = cone.second_angle

    def flange_azimuth(first):
        """How far round the volute the flange ends up, from the outlet."""
        end, _, _ = _place(start, tangential,
                           _direction(first, second, tangential, radial, axial),
                           arc_length)
        return math.atan2(np.dot(end, tangential), np.dot(end, radial))

    # Radial is implicit: turning further carries the flange further round, and
    # carrying it round brings radial closer to where the duct already points.
    def radial_residual(first):
        return first + flange_azimuth(first) - math.pi/2.0

    radial_angle = brentq(radial_residual, 0.0, math.pi/2.0, xtol=1e-12, rtol=1e-12)

    if cone.first_angle == "tangential":
        first = 0.0
    elif cone.first_angle == "radial":
        first = radial_angle
    else:
        first = float(cone.first_angle)
        if first > radial_angle + 1.0e-9:
            raise ValueError(
                f"first_angle {math.degrees(first):.1f} deg exceeds the radial "
                f"angle {math.degrees(radial_angle):.1f} deg for this duct, so "
                "the outlet would point back inward. Radial is closer than a "
                "right angle because the duct travels round the volute; "
                "shortening arc_length_parameter moves it further out."
            )
        first = min(first, radial_angle)

    end_direction = _direction(first, second, tangential, radial, axial)
    end, turn, bend_radius = _place(start, tangential, end_direction, arc_length)
    duct_radius = max(max_section_radius, cone.final_diameter/2.0)
    if bend_radius <= duct_radius:
        raise ValueError(
            f"a {math.degrees(turn):.1f} deg turn in {arc_length*1e3:.1f} mm "
            f"gives a centreline bend radius of {bend_radius*1e3:.1f} mm, "
            f"inside the duct's own {duct_radius*1e3:.1f} mm radius, so the "
            "duct would fold through itself. Raise arc_length_parameter."
        )
    return ResolvedExitCone(
        start=start, end=end, start_direction=tangential,
        end_direction=end_direction, first_angle=first, second_angle=second,
        radial_angle=radial_angle, end_azimuth=flange_azimuth(first),
        arc_length=arc_length, turn_angle=turn, bend_radius=bend_radius,
        outlet_diameter=outlet_diameter, final_diameter=cone.final_diameter,
    )
