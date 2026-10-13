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

Tangential means the direction the collected flow is travelling, so it reverses
with the volute's rotation sense and a clockwise volute's duct leaves the other
way round. Radial and the elevation do not reverse: they are fixed by the
machine axis rather than by which way the flow runs round it.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.optimize import brentq

from .path import _rotate
from .port import Port, PortSite, ResolvedPort, _STRAIGHT_TOLERANCE

__all__ = [
    "ExitCone",
    "ResolvedExitCone",
]

#: Directions the first angle may name instead of a number.
NAMED_ANGLES = ("tangential", "radial")


@dataclass(frozen=True)
class ExitCone(Port):
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
    Exit-cone geometry reaches the loss analysis three ways.
    :attr:`final_diameter` sets station 7 and, with the arc length, the angle
    the cone opens out at -- which is most of what decides its diffusion loss.
    The arc length also sets the wall area friction is charged over, and the
    turn sets the bend loss.

    That makes :attr:`arc_length_parameter` a real design variable rather than
    a drawing detail: a longer duct opens out more gradually and recovers more
    of the velocity head, at the cost of more wetted wall and, for a given
    turn, a gentler bend. :meth:`~pyskyfire.ducts.VoluteResult.loss_breakdown`
    reports the three apart so the trade can be seen.
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

    def resolve(self, site):
        """Place this cone on a passage. See :class:`ResolvedExitCone`.

        The arc length is a multiple of the bore the duct leaves from, so a
        cone carries its length with it and takes its scale from wherever it
        is fitted.
        """
        return resolve_exit_cone(
            self, start=site.start, outlet_diameter=site.outlet_diameter,
            max_section_radius=site.max_section_radius, sign=site.sign,
        )


@dataclass(frozen=True)
class ResolvedExitCone(ResolvedPort):
    """An exit cone placed in machine coordinates.

    Everything a duct has in common with any other is on
    :class:`~pyskyfire.ducts.port.ResolvedPort`; what a cone adds is the pair
    of angles it was steered by, and where round the volute they carried it.
    ``start_direction`` is tangential at the volute outlet.

    Attributes
    ----------
    first_angle, second_angle : float
        The resolved turn within the volute's plane and the elevation [rad].
    radial_angle : float
        The turn that would have put the flange radial [rad]. Equal to
        :attr:`first_angle` when ``'radial'`` was asked for, and the upper
        bound on it otherwise.
    end_azimuth : float
        Azimuth the flange ends up at, measured from the volute outlet [rad].
    """

    first_angle: float
    second_angle: float
    radial_angle: float
    end_azimuth: float


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


def resolve_exit_cone(cone, *, start, outlet_diameter, max_section_radius,
                      sign=1.0):
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
    sign : float, optional
        Rotation sense of the volute the duct leaves, from
        :func:`~pyskyfire.common.rotation.rotation_sign`. Tangential means the
        direction the collected flow is travelling, so it reverses with the
        sense; radial and the elevation do not, being fixed by the machine
        axis rather than by which way the flow goes round it.

    Returns
    -------
    ResolvedExitCone
    """
    start = np.asarray(start, dtype=float)
    if start.shape != (3,) or not np.isfinite(start).all():
        raise ValueError("start must be a finite (3,) point")
    azimuth = math.atan2(start[1], start[0])
    tangential = sign * np.array(
        [-math.sin(azimuth), math.cos(azimuth), 0.0])
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
