"""How a collecting passage hands its flow to a pipe.

A volute and a chamber manifold both end at a flange, and both need a duct to
get there, but they leave in different directions. A volute's flow is already
travelling round the passage, so its duct leaves tangentially and turns; a
manifold has no such direction, so its duct runs straight out of the dome its
section is built on, aimed wherever the packaging wants it. Those are two ports,
not two classes of machine, and either passage may carry either one -- a
manifold fed tangentially is a real arrangement, and it should not cost a second
implementation of the same duct.

So a port is stated in two parts. :class:`Port` is the request, carrying
whatever the designer chose. :class:`PortSite` is what the passage knows about
where the port attaches, handed over at resolve time. The result is a
:class:`ResolvedPort`, which is the duct placed in machine coordinates and is
all the loss model and the mesher ever read.

The site carries more than any one port uses. A tangential duct wants the
rotation sense, since tangential reverses with it; a radial one wants the dome
it leaves through. Each port asks for the fields it needs and says so plainly
when they are missing, which is what lets a passage offer both without knowing
which it has been given.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
import math

import numpy as np

from .path import _rotate

__all__ = ["Port", "PortSite", "DomeSite", "ResolvedPort", "ResolvedRadialPort",
           "RadialPort"]

#: Below this turn a duct is treated as straight rather than as a huge arc.
_STRAIGHT_TOLERANCE = 1.0e-9

#: Angles this far past a limit are taken as sitting on it.
_ANGLE_TOLERANCE = 1.0e-12


@dataclass(frozen=True)
class DomeSite:
    """The half-ellipse a stub leaves through, in machine coordinates.

    The ellipse is the one the passage's own section is built on, centred where
    its two axes cross, which is the point a stub is aimed from. Only part of it
    is passage boundary: past :attr:`arc_limit` the section's fillets have taken
    over and the ellipse runs outboard of the metal, so a stub aimed out there
    would leave through the side rather than the dome.

    Attributes
    ----------
    normal : np.ndarray
        The way the section faces [shape (3,)], along the ellipse's own
        first axis.
    in_plane : np.ndarray
        The section plane's other axis [shape (3,)], across the normal and
        along the wall the section sits against.
    semi_normal, semi_in_plane : float
        Ellipse semi-axes along those two directions [m].
    arc_limit : float
        Ellipse parameter the dome runs to either side of ``normal`` [rad].
        Under a quarter turn, since the fillets take over before the ellipse
        reaches its widest.
    tilt : float, optional
        Angle of ``normal`` from machine-radial [rad], positive toward ``+z``.
        A stub stated against the plane the passage rings rather than against
        the section is turned through this.
    shrinkage : float, optional
        How fast the dome's depth falls away either side of the port
        [m/rad of azimuth]. Read only when drawing, to sink a stub far enough
        that the passage does not show daylight under it.
    """

    normal: np.ndarray
    in_plane: np.ndarray
    semi_normal: float
    semi_in_plane: float
    arc_limit: float
    tilt: float = 0.0
    shrinkage: float = 0.0

    def __post_init__(self):
        for name in ("normal", "in_plane"):
            vector = np.asarray(getattr(self, name), dtype=float)
            norm = np.linalg.norm(vector)
            if (vector.shape != (3,) or not np.isfinite(vector).all()
                    or norm == 0.0):
                raise ValueError(f"{name} must be a finite nonzero (3,) direction")
            object.__setattr__(self, name, vector/norm)
        for name in ("semi_normal", "semi_in_plane", "arc_limit"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
            object.__setattr__(self, name, value)
        if self.arc_limit >= math.pi/2.0:
            raise ValueError(
                "arc_limit must be under a quarter turn: the fillets take over "
                "before the ellipse reaches its widest point")
        object.__setattr__(self, "tilt", float(self.tilt))
        object.__setattr__(self, "shrinkage", float(self.shrinkage))

    @property
    def angle_limit(self) -> float:
        """Furthest off the normal a stub may be aimed and still clear the fillet [rad]."""
        return math.atan2(self.semi_in_plane*math.sin(self.arc_limit),
                          self.semi_normal*math.cos(self.arc_limit))

    def direction(self, angle):
        """Unit direction of a stub ``angle`` off the normal [shape (3,)]."""
        return (math.cos(angle)*self.normal + math.sin(angle)*self.in_plane)

    def across(self, angle):
        """Unit direction across such a stub, within the section plane."""
        return (-math.sin(angle)*self.normal + math.cos(angle)*self.in_plane)

    def exit_distance(self, angle):
        """How far the dome sits from the ellipse centre along a stub [m]."""
        return 1.0/math.hypot(math.cos(angle)/self.semi_normal,
                              math.sin(angle)/self.semi_in_plane)

    def clearance(self, angle):
        """Half the room the dome leaves across a stub aimed ``angle`` off the normal [m].

        The tangents to the ellipse parallel to the stub bound it where they
        touch the dome itself. Where a tangent touches out past
        :attr:`arc_limit` it bounds nothing, because the metal there is fillet
        rather than dome, and the arc's own end takes over as the bound. So this
        is how far the dome arc reaches across the stub, on whichever side
        reaches less.
        """
        sin_angle, cos_angle = math.sin(angle), math.cos(angle)
        semi_n, semi_p, limit = (self.semi_normal, self.semi_in_plane,
                                 self.arc_limit)

        def offset(parameter):
            return (-semi_n*math.cos(parameter)*sin_angle
                    + semi_p*math.sin(parameter)*cos_angle)

        # Where the ellipse's tangent runs parallel to the stub, which is where
        # the arc reaches furthest across it when it reaches there at all.
        touch = math.atan2(semi_p*cos_angle, -semi_n*sin_angle)
        parameters = [limit, -limit]
        parameters += [touch + turn for turn in (-math.pi, 0.0, math.pi)
                       if -limit < touch + turn < limit]
        offsets = [offset(parameter) for parameter in parameters]
        return min(max(offsets), -min(offsets))

    def generator_exit(self, angle, offsets):
        """Where lines parallel to a stub leave the dome [m].

        ``offsets`` are distances across the stub within the section plane, as
        the walls of a round tube sit at. What comes back is each line's
        crossing measured along the stub from where its own centreline crosses,
        which is how far the tube's wall has to run back to meet the passage.
        """
        semi_n, semi_p = self.semi_normal, self.semi_in_plane
        sin_angle, cos_angle = math.sin(angle), math.cos(angle)
        offsets = np.asarray(offsets, dtype=float)
        distance = self.exit_distance(angle)
        along = distance*cos_angle - offsets*sin_angle
        across = distance*sin_angle + offsets*cos_angle
        quadratic = cos_angle**2/semi_n**2 + sin_angle**2/semi_p**2
        linear = 2.0*(along*cos_angle/semi_n**2 + across*sin_angle/semi_p**2)
        constant = along**2/semi_n**2 + across**2/semi_p**2 - 1.0
        root = np.sqrt(np.maximum(linear**2 - 4.0*quadratic*constant, 0.0))
        return (-linear + root)/(2.0*quadratic)

    def deepest_generator(self, angle, radius):
        """The furthest along a stub of ``radius`` its own dome reaches [m].

        Measured from where the centreline leaves the dome. Normally negative,
        because a curved dome falls away from the centreline's crossing to
        either side; positive for a stub tilted far enough that the dome's crest
        in the stub's own direction sits under the flange rather than on its
        axis.
        """
        sin_angle, cos_angle = math.sin(angle), math.cos(angle)
        semi_n, semi_p = self.semi_normal, self.semi_in_plane
        crest = (sin_angle*cos_angle*(semi_p**2 - semi_n**2)
                 / math.hypot(semi_n*cos_angle, semi_p*sin_angle))
        # The crossing is unimodal in the offset, so the deepest point under the
        # flange is the crest where the flange covers it and the rim otherwise.
        offsets = np.array([min(max(crest, -radius), radius), -radius, radius])
        return float(self.generator_exit(angle, offsets).max())


@dataclass(frozen=True)
class PortSite:
    """What a passage knows about where its port attaches.

    Attributes
    ----------
    start : np.ndarray
        Centre of the bore the port leaves from [m, shape (3,)]. For a section
        built on a dome this is the ellipse centre, which is where a stub is
        aimed from rather than where it starts.
    outlet_diameter : float
        Diameter of that bore [m].
    max_section_radius : float
        Furthest the passage's own section reaches from ``start`` [m], used to
        reject a bend tight enough to fold the duct through itself.
    sign : float, optional
        Rotation sense of the passage, from
        :func:`~pyskyfire.common.rotation.rotation_sign`. Tangential reverses
        with it; radial does not.
    dome : DomeSite, optional
        The half-ellipse a stub would leave through. A volute has no need of
        it, and neither has a section family that offers a stub nowhere to sit.
    """

    start: np.ndarray
    outlet_diameter: float
    max_section_radius: float
    sign: float = 1.0
    dome: DomeSite | None = None

    def require(self, name):
        """Return a field a port needs, or say which port wanted it."""
        value = getattr(self, name)
        if value is None:
            raise ValueError(
                f"this port needs the attachment site's {name!r}, which the "
                "passage it is being fitted to does not supply")
        return value


@dataclass(frozen=True)
class ResolvedPort:
    """A duct placed in machine coordinates.

    Attributes
    ----------
    start, end : np.ndarray
        Centres of the passage outlet and of the flange [m, shape (3,)].
    start_direction, end_direction : np.ndarray
        Unit normals there [shape (3,)].
    arc_length, turn_angle, bend_radius : float
        Centreline length [m], total turn [rad], and centreline bend radius
        [m], which is infinite for a straight duct.
    outlet_diameter, final_diameter : float
        Bore the duct leaves from and the flange it delivers to [m].
    """

    start: np.ndarray
    end: np.ndarray
    start_direction: np.ndarray
    end_direction: np.ndarray
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
            Override for the frame at the outlet, as ``(u0, v0)`` unit vectors
            perpendicular to :attr:`start_direction`. Defaults to
            :meth:`_start_frame`. A passage whose section is tilted passes its
            own basis here instead, so the duct's first ring is bit-identical
            to the passage's last one rather than merely tangent to the same
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

    def sweep(self, n_points=64, n_stations=9):
        """The duct's own wall, as closed boundary rings [m, shape (N, M, 3)].

        The duct's own frame, not the passage's: a stub leaves along the way its
        section faces, so that section's outward direction *is* the duct axis
        and spanning the cross-section with it would flatten the tube into the
        plane it should be normal to.
        """
        n_points = int(n_points)
        if n_points < 3:
            raise ValueError("n_points must be at least three")
        points, u, v = self.frames(n_stations)
        angles = np.linspace(0.0, 2.0*math.pi, n_points, endpoint=False)
        spoke = (np.cos(angles)[None, :, None]*u[:, None, :]
                 + np.sin(angles)[None, :, None]*v[:, None, :])
        return points[:, None, :] + (self.final_diameter/2.0)*spoke

    def _axis(self):
        axis = np.cross(self.start_direction, self.end_direction)
        return axis / np.linalg.norm(axis)

    def _start_frame(self):
        """A frame spanning the plane the duct starts in.

        Radial is taken from where the outlet sits rather than from the
        direction it leaves along, so a clockwise volute -- whose tangential
        runs the other way -- still gets radial pointing outward. Where the
        duct leaves along a direction radial is not perpendicular to, as a
        manifold's does, radial is squared up against it first.
        """
        axial = np.array([0.0, 0.0, 1.0])
        radial = self.start - np.dot(self.start, axial) * axial
        if np.linalg.norm(radial) < 1.0e-12:
            radial = np.array([1.0, 0.0, 0.0])
        direction = self.start_direction
        radial = radial - np.dot(radial, direction) * direction
        norm = np.linalg.norm(radial)
        if norm < 1.0e-12:
            fallback = (np.array([1.0, 0.0, 0.0]) if abs(direction[0]) < 0.9
                        else np.array([0.0, 1.0, 0.0]))
            radial = fallback - np.dot(fallback, direction) * direction
            norm = np.linalg.norm(radial)
        radial = radial / norm
        return radial, np.cross(radial, direction)


class Port(ABC):
    """A request for the duct carrying collected flow to its flange."""

    @abstractmethod
    def resolve(self, site: PortSite) -> ResolvedPort:
        """Place this port on a passage, given what the passage knows."""


@dataclass(frozen=True)
class ResolvedRadialPort(ResolvedPort):
    """A stub placed in machine coordinates, with the dome it leaves through.

    Attributes
    ----------
    angle : float
        Resolved tilt off the section's normal [rad], whichever reference the
        stub was stated against.
    angle_reference : str
        The reference it was stated against.
    centre : np.ndarray
        Ellipse centre the stub was aimed from [m, shape (3,)]. It is
        :attr:`~ResolvedPort.start` that sits on the dome, one exit distance
        further along; the length the designer gave is measured from here.
    dome : DomeSite
        The half-ellipse it was fitted to.
    """

    angle: float
    angle_reference: str
    centre: np.ndarray
    dome: DomeSite

    def _start_frame(self):
        """Across the section plane first, then out of it.

        A stub's cross-section splits naturally into the direction the dome
        curves away in, which lies in the section plane, and the azimuthal one
        the swept passage carries straight on through. Trimming the tube to the
        passage only moves it along the first.
        """
        across = self.dome.across(self.angle)
        return across, np.cross(self.start_direction, across)

    def sweep(self, n_points=64, n_stations=9):
        """The stub's wall, trimmed back onto the dome it leaves [m, (N, M, 3)].

        The first ring follows the dome rather than lying flat against it, which
        is what stops the tube floating a chord above the passage. It is sunk a
        little further than this station's own ellipse asks for, because the
        section shrinks either side of the port and a tube buried slightly in
        metal reads better than one with daylight under it.
        """
        n_points, n_stations = int(n_points), int(n_stations)
        if n_points < 3:
            raise ValueError("n_points must be at least three")
        if n_stations < 3:
            raise ValueError("n_stations must be at least three")
        radius = self.final_diameter/2.0
        across, azimuthal = self._start_frame()
        angles = np.linspace(0.0, 2.0*math.pi, n_points, endpoint=False)
        spoke = radius*(np.cos(angles)[:, None]*across
                        + np.sin(angles)[:, None]*azimuthal)
        sunk = (self.dome.generator_exit(self.angle, radius*np.cos(angles))
                - self._drawing_margin())
        rim = self.start + spoke + sunk[:, None]*self.start_direction
        # The rest of the tube starts ahead of every trimmed vertex, so the ring
        # that closes the trim band cannot invert against it.
        travel = np.linspace(max(0.0, float(sunk.max())), self.arc_length,
                             n_stations - 1)
        tube = (self.start
                + travel[:, None]*self.start_direction)[:, None, :] + spoke
        return np.concatenate((rim[None], tube))

    def _drawing_margin(self):
        """How much further than its own ellipse to sink the first ring [m].

        The trim is taken against the port station's ellipse alone, which is
        cheap but leaves a gap: the passage is swept round a circle and its
        section shrinks away from the port, so both effects pull the real
        surface inboard of that ellipse across the width the flange covers.
        """
        radius = float(np.hypot(self.centre[0], self.centre[1]))
        azimuth = (self.final_diameter/2.0)/radius
        return radius*(1.0 - math.cos(azimuth)) + self.dome.shrinkage*azimuth


@dataclass(frozen=True)
class RadialPort(Port):
    """A straight stub leaving the dome a manifold's section is built on.

    This is the chamber-manifold outlet. There is no swirl to take off
    tangentially, so the duct runs straight from the dome to its flange. It is
    aimed from the ellipse centre -- where the section's two axes cross, which
    is a defined point in a way the curved surface is not -- and steered by a
    single angle in the section's own plane.

    Parameters
    ----------
    length : float
        How far the flange sits from the ellipse centre [m]. Measured from
        there rather than from the surface, so it does not move when the angle
        or the area law does; it must be long enough to put the whole flange
        clear of the dome.
    angle : float, optional
        Tilt in the section plane [rad], positive toward ``+z``, matching how
        :meth:`~pyskyfire.ducts.half_moon.HalfMoonSection.meridional_frame`
        reports the section's own lean. Zero leaves along whichever reference
        was chosen.
    angle_reference : {'normal', 'manifold_plane'}, optional
        What the angle is measured from. ``'normal'`` is the way the section
        faces, so zero points straight out of the dome however the wall happens
        to lean there. ``'manifold_plane'`` is the plane the manifold rings,
        so zero points straight out from the machine axis and the same number
        gives the same hardware direction at every station.
    diameter : float, optional
        Flange diameter outright [m]. May not exceed the room the dome leaves
        across the stub, since a wider one would leave through the fillet.
    diameter_fraction : float, optional
        Flange diameter as a fraction of that room, used when ``diameter`` is
        not given. The default fills the dome nearly to its edge, leaving a
        little metal at the corners.

    Notes
    -----
    The angle both aims the stub and sizes it. How much room the dome leaves
    across a tilted stub is not the ellipse's own semi-axis: it is how far the
    part of the ellipse that is actually passage boundary reaches across the
    direction chosen, which closes down as the stub turns toward the fillets.
    The stub stays centred on the ray it was aimed along, so a tilted one is
    sized by whichever side of it runs out of dome first; letting it shift off
    that ray to take the room on the other side is a further freedom, not taken
    here.
    """

    length: float
    angle: float = 0.0
    angle_reference: str = "normal"
    diameter: float | None = None
    diameter_fraction: float = 0.95

    def __post_init__(self):
        length = float(self.length)
        if not math.isfinite(length) or length <= 0.0:
            raise ValueError("length must be finite and positive")
        object.__setattr__(self, "length", length)
        angle = float(self.angle)
        if not math.isfinite(angle):
            raise ValueError("angle must be finite")
        object.__setattr__(self, "angle", angle)
        if self.angle_reference not in ("normal", "manifold_plane"):
            raise ValueError(
                "angle_reference must be 'normal', measured from the way the "
                "section faces, or 'manifold_plane', measured from the plane "
                "the manifold rings")
        if self.diameter is not None:
            diameter = float(self.diameter)
            if not math.isfinite(diameter) or diameter <= 0.0:
                raise ValueError("diameter must be finite and positive")
            object.__setattr__(self, "diameter", diameter)
        fraction = float(self.diameter_fraction)
        if not math.isfinite(fraction) or not 0.0 < fraction <= 1.0:
            raise ValueError(
                "diameter_fraction must lie in (0, 1]: a stub wider than the "
                "dome it sits on would leave the passage through its side")
        object.__setattr__(self, "diameter_fraction", fraction)

    def resolve(self, site: PortSite) -> ResolvedRadialPort:
        dome = site.require("dome")
        angle = (self.angle if self.angle_reference == "normal"
                 else self.angle - dome.tilt)
        limit = dome.angle_limit
        if abs(angle) > limit + _ANGLE_TOLERANCE:
            raise ValueError(
                f"a stub {math.degrees(angle):.1f} deg off the section normal "
                f"would leave through the fillet rather than the dome, which "
                f"carries one to {math.degrees(limit):.1f} deg either side"
                + ("" if self.angle_reference == "normal" else
                   f" and leans {math.degrees(dome.tilt):.1f} deg itself"))

        clearance = dome.clearance(angle)
        if clearance <= 0.0:
            raise ValueError(
                f"the dome leaves no room across a stub {math.degrees(angle):.1f} "
                "deg off its normal, which is aimed at the very end of it")
        if self.diameter is None:
            diameter = self.diameter_fraction*2.0*clearance
        else:
            diameter = self.diameter
            if diameter > 2.0*clearance:
                raise ValueError(
                    f"a {diameter*1e3:.1f} mm stub does not fit the "
                    f"{2.0*clearance*1e3:.1f} mm the dome leaves across a stub "
                    f"{math.degrees(angle):.1f} deg off its normal; it would "
                    "leave through the fillet rather than the dome")

        centre = np.asarray(site.start, dtype=float)
        if centre.shape != (3,) or not np.isfinite(centre).all():
            raise ValueError("start must be a finite (3,) point")
        direction = dome.direction(angle)
        surface = dome.exit_distance(angle)
        flange = surface + dome.deepest_generator(angle, diameter/2.0)
        if self.length <= flange:
            raise ValueError(
                f"a {self.length*1e3:.1f} mm stub does not clear the dome, "
                f"which its flange reaches out to {flange*1e3:.1f} mm from the "
                "ellipse centre the length is measured from; lengthen it")
        return ResolvedRadialPort(
            start=centre + surface*direction,
            end=centre + self.length*direction,
            start_direction=direction, end_direction=direction,
            arc_length=self.length - surface, turn_angle=0.0,
            bend_radius=math.inf,
            outlet_diameter=float(site.outlet_diameter),
            final_diameter=diameter,
            angle=angle, angle_reference=self.angle_reference,
            centre=centre, dome=dome,
        )
