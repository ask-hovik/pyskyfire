"""A ring that shares flow out to a circumference, or gathers it back in.

Where a volute collects swirling flow and is sized to conserve the angular
momentum it takes in, a manifold has no swirl to conserve. Fluid arrives
through a port and has to travel round to reach the far channels, or comes in
all round and has to travel to the port; either way the circumferential velocity
is a cost rather than a design variable. So a manifold is not a volute with a
different area law, and it is not named the way one is.

Branches, not a direction
-------------------------
A volute can only run one way, because the swirl it collects already has a
sense. A manifold can be built either way, and the choice is a topology rather
than a direction.

``branches=1``
    The flow runs one way round the whole ring, from a dead point at the tongue
    to the port a full revolution later. This needs a ``rotation``, because
    which way it runs is a real choice about the hardware.
``branches=2``
    The flow splits at a dead point and runs both ways, each half gathering its
    own semicircle and meeting at the port. This is mirror-symmetric, so no
    rotation applies and asking for one is refused: there is no sense to name.

Two branches is the usual arrangement for a chamber manifold, because it halves
the distance the fluid travels and so the head it spends getting there.

Multiple ports are the obvious next step and would put a dead point between
each neighbouring pair. The machinery counts branches rather than assuming two,
so that arrives as ``branches=2*n_ports`` with the dead points steerable, not as
a second class.

What is resolved
----------------
Only one branch. Every branch of a given manifold is the same arc sized the
same way, so the ring is that one arc repeated -- mirrored for a split
manifold, rotated for a multi-port one -- and :meth:`Manifold.sweep` is what
puts the whole ring together.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from pyskyfire.common.rotation import parse_rotation, rotation_sign

from .manifold_law import ManifoldLaw
from .port import DomeSite, Port, PortSite
from .volute_section import (VoluteSection, section_frames,
                             solve_section_size)

__all__ = ["Manifold", "ManifoldGeometry"]


@dataclass(frozen=True)
class ManifoldGeometry:
    """Resolved passage geometry along one branch.

    Attributes
    ----------
    theta : np.ndarray
        Angle from the dead point [rad, shape (N,)], ending at the port.
    collected : np.ndarray
        Fraction of the manifold's whole mass flow through each station,
        ending at ``1/branches`` because a branch only ever carries its share.
    size, aspect, area : np.ndarray
        Section size parameter [m], aspect ratio, and passage area [m^2].
    axial, radial, tilt : np.ndarray
        Where each section sits and how it leans, from the section family.
    clipped : np.ndarray
        True where the law asked for less area than the family can build even
        at its smallest size and flattest aspect [bool, shape (N,)].
    """

    theta: np.ndarray
    collected: np.ndarray
    size: np.ndarray
    aspect: np.ndarray
    area: np.ndarray
    axial: np.ndarray
    radial: np.ndarray
    tilt: np.ndarray
    clipped: np.ndarray

    @property
    def dead_clip_angle(self) -> float:
        """Angle up to which the law could not be satisfied [rad].

        Zero when it is satisfied everywhere. A finite value means the passage
        near the dead point is larger than the law wants, because the section
        family cannot shrink any further -- the manifold's version of a
        volute's tongue.
        """
        if not self.clipped.any():
            return 0.0
        return float(self.theta[self.clipped].max())


class Manifold:
    """A distributing or collecting ring, resolved from a section and a law.

    Parameters
    ----------
    section : VoluteSection
        Cross-section family. A
        :class:`~pyskyfire.ducts.half_moon.HalfMoonSection` carries its own
        placement on the chamber wall; any other family needs ``radius``.
    area_law : ManifoldLaw
        How the passage grows as flow is gathered. Stated rather than
        defaulted, because a uniform ring and a tapered one are different
        pieces of hardware.
    branches : int, optional
        How many arcs the circumference is divided into. One runs the whole way
        round; two split at a dead point and meet at the port.
    rotation : str, optional
        Which way a single-branch manifold runs. Required for ``branches=1``
        and refused for anything else, which is mirror-symmetric and has no
        sense to name.
    port : Port, optional
        The duct carrying flow to or from its flange, normally a
        :class:`~pyskyfire.ducts.port.RadialPort`, which leaves the dome the
        section is built on and may be steered within the section plane. A
        manifold fed tangentially may be given an
        :class:`~pyskyfire.ducts.exit_cone.ExitCone` instead.
    radius, tilt : float, optional
        Where the section sits, for a family that does not place itself.
    n_stations : int, optional
        Stations along one branch, including both the dead point and the port.

    Notes
    -----
    Geometry is resolved on construction and carries no operating point, so one
    manifold can be reported at several mass flows.
    """

    def __init__(self, section, area_law, *, branches=2, rotation=None,
                 port=None, radius=None, tilt=0.0, n_stations=181):
        if not isinstance(section, VoluteSection):
            raise TypeError("section must be a VoluteSection")
        if not isinstance(area_law, ManifoldLaw):
            raise TypeError(
                "area_law must be a ManifoldLaw. A volute's AreaLaw answers to "
                "swirl, which a manifold has none of")
        branches = int(branches)
        if branches not in (1, 2):
            raise ValueError(
                "branches must be 1 or 2; more arrive with multiple ports, "
                "which are not implemented yet")
        if branches == 1:
            if rotation is None:
                raise ValueError(
                    "a single-branch manifold runs one way round, so it needs "
                    "a rotation")
            self.rotation = parse_rotation(rotation)
        else:
            if rotation is not None:
                raise ValueError(
                    f"a {branches}-branch manifold is mirror-symmetric, so it "
                    "has no rotation to name; the flow runs both ways from its "
                    "dead point")
            self.rotation = None
        if port is not None and not isinstance(port, Port):
            raise TypeError("port must be a Port, such as a RadialPort")
        n_stations = int(n_stations)
        if n_stations < 3:
            raise ValueError("n_stations must be at least three")

        self.section = section
        self.area_law = area_law
        self.branches = branches
        self._sign = 1.0 if self.rotation is None else rotation_sign(self.rotation)
        self._radius = radius
        self._tilt = float(tilt)
        self.geometry = self._resolve(n_stations)
        self.port = port
        self.port_geometry = None if port is None else self._resolve_port(port)

    # -- geometry ---------------------------------------------------------

    def _place(self, sizes, aspects):
        """Where each section sits, from the family or from the manifold."""
        if hasattr(self.section, "meridional_frame"):
            if self._radius is not None:
                raise ValueError(
                    f"{type(self.section).__name__} places itself on the wall "
                    "it was built against, so radius must not also be given")
            return self.section.meridional_frame(sizes, aspects)
        if self._radius is None:
            raise ValueError(
                f"{type(self.section).__name__} does not place itself, so the "
                "manifold needs a radius to sit it on")
        shape = np.shape(sizes)
        return (np.zeros(shape), np.full(shape, float(self._radius)),
                np.full(shape, self._tilt))

    def _resolve(self, n_stations):
        span = 2.0*math.pi/self.branches
        theta = np.linspace(0.0, span, n_stations)
        collected = theta/(2.0*math.pi)
        sizes = np.empty(n_stations)
        aspects = np.empty(n_stations)
        clipped = np.zeros(n_stations, dtype=bool)
        for index, gathered in enumerate(collected):
            def residual(size, aspect, gathered=float(gathered)):
                area, _, _ = self.section.metrics(size, aspect)
                return self.area_law.residual(gathered, float(area[0]))
            sizes[index], aspects[index], clipped[index] = solve_section_size(
                self.section, residual)
        area, _, _ = self.section.metrics(sizes, aspects)
        axial, radial, tilt = self._place(sizes, aspects)
        return ManifoldGeometry(
            theta=theta, collected=collected, size=sizes, aspect=aspects,
            area=area, axial=axial, radial=radial, tilt=tilt, clipped=clipped)

    def display_geometry(self, n_stations=None):
        """Branch geometry on a chosen station grid, for drawing only.

        A swept manifold facets in azimuth long before it facets around its
        section -- the ring is metres round and the section centimetres -- so
        a smooth picture needs more stations than the area law was resolved on,
        and needing them is a drawing matter rather than a design one.

        The sizes the law resolved are read between stations, and everything
        that follows from a size is then computed at the new angles exactly:
        where the section sits, how it leans, and the area it ends up with. So
        the passage drawn agrees with the resolved one at every design station
        and is interpolated only in between. ``None`` gives the resolved grid
        itself, which is what the design says and nothing more.

        Returns
        -------
        ManifoldGeometry
            Geometry along one branch, on the requested grid.
        """
        geometry = self.geometry
        if n_stations is None:
            return geometry
        n_stations = int(n_stations)
        if n_stations < 3:
            raise ValueError("n_stations must be at least three")
        if n_stations == len(geometry.theta):
            return geometry
        theta = np.linspace(0.0, 2.0*math.pi/self.branches, n_stations)
        sizes = np.interp(theta, geometry.theta, geometry.size)
        aspects = np.interp(theta, geometry.theta, geometry.aspect)
        area, _, _ = self.section.metrics(sizes, aspects)
        axial, radial, tilt = self._place(sizes, aspects)
        # A station between two clipped ones is clipped too; the flag marks a
        # stretch of the ring, not a grid point.
        clipped = np.interp(theta, geometry.theta,
                            geometry.clipped.astype(float)) > 0.0
        return ManifoldGeometry(
            theta=theta, collected=theta/(2.0*math.pi), size=sizes,
            aspect=aspects, area=area, axial=axial, radial=radial, tilt=tilt,
            clipped=clipped)

    def sweep_stations(self, n_stations=None):
        """How many stations :meth:`sweep` returns on a given grid.

        A split manifold sweeps both branches, so this is more than the branch
        that was resolved has: the two meet at the dead point and share it.
        """
        per_branch = len(self.display_geometry(n_stations).theta)
        return self.branches*(per_branch - 1) + 1

    def _ring(self, n_stations=None):
        """One branch expanded to the whole circumference.

        A split manifold is its branch and that branch's mirror image, which in
        azimuth is the same arc run backwards. Nothing about the section itself
        mirrors: the plane the mirror acts in contains the machine axis, so it
        carries each section plane onto another one and leaves the shape in it
        alone.
        """
        geometry = self.display_geometry(n_stations)
        fields = (geometry.theta, geometry.size, geometry.aspect,
                  geometry.axial, geometry.radial, geometry.tilt)
        if self.branches == 1:
            return fields
        # Everything but the angle repeats under the mirror; the angle negates.
        back = slice(None, 0, -1)
        return tuple(np.concatenate((-values[back] if values is geometry.theta
                                     else values[back], values))
                     for values in fields)

    def sweep(self, n_points=64, n_stations=None):
        """Section boundaries swept into machine coordinates.

        The machine axis is ``z`` and the dead point sits at azimuth zero. A
        split manifold is swept over the whole ring, both branches together, so
        what comes back is the closed cavity rather than half of it.

        ``n_points`` samples each boundary piece of the section; ``n_stations``
        asks for a finer grid along the branch than the law was resolved on,
        for drawing -- see :meth:`display_geometry`.

        Returns
        -------
        np.ndarray
            Closed boundary rings [m, shape (N, M, 3)], one per station.
        """
        theta, sizes, aspects, axial, radial_offset, tilt = self._ring(n_stations)
        ring = self.section.polygon(sizes, n_points=n_points, aspect=aspects)
        radial, u_hat, v_hat = section_frames(theta, tilt=tilt, sign=self._sign)
        axis = np.zeros_like(radial)
        axis[:, 2] = 1.0
        origin = radial_offset[:, None]*radial + axial[:, None]*axis
        u, v = ring[..., 0], ring[..., 1]
        return (origin[:, None, :]
                + u[..., None]*u_hat[:, None, :]
                + v[..., None]*v_hat[:, None, :])

    def sweep_slit(self, n_points=16, n_stations=None):
        """The reserved slit swept into machine coordinates.

        A strip inside the chamber-facing face, not a separate face: it is
        where the channels break out, and the rest of that face is metal
        against the wall. Drawn apart so a layout can be checked against where
        the passes actually start and end.

        Returns
        -------
        np.ndarray
            Points [m, shape (N, n_points, 3)], one row per swept station.
        """
        return self._sweep_piece(self._reserved("slit_boundary", n_stations),
                                 n_points, n_stations)

    def sweep_margins(self, n_points=16, n_stations=None):
        """The chamber-facing wall either side of the slit [m].

        Returns the two strips in the order the face is traversed, so with
        :meth:`sweep_slit` between them they tile it without overlapping.
        """
        return tuple(self._sweep_piece(curve, n_points, n_stations)
                     for curve in self._reserved("margin_boundaries", n_stations))

    def _reserved(self, name, n_stations=None):
        """A boundary piece the section reserves, sized for the whole ring."""
        if not hasattr(self.section, name):
            raise ValueError(
                f"{type(self.section).__name__} reserves no slit; only a "
                "section measured against a wall has one")
        _, sizes, aspects, *_ = self._ring(n_stations)
        return getattr(self.section, name)(sizes, aspects)

    def _sweep_piece(self, curve, n_points, n_stations=None):
        """Place one boundary piece of every station into machine coordinates."""
        theta, _, _, axial, radial_offset, tilt = self._ring(n_stations)
        local = curve.evaluate(np.linspace(0.0, 1.0, int(n_points)))
        radial, u_hat, v_hat = section_frames(theta, tilt=tilt, sign=self._sign)
        axis = np.zeros_like(radial)
        axis[:, 2] = 1.0
        origin = radial_offset[:, None]*radial + axial[:, None]*axis
        return (origin[:, None, :]
                + local[..., :1]*u_hat[:, None, :]
                + local[..., 1:2]*v_hat[:, None, :])

    def fan_origins(self, n_stations=None):
        """Star-shaped fan origins for each swept station, in machine coordinates.

        One per station of :meth:`sweep`, so an end face can be closed off as a
        fan about a point that sees its whole boundary.
        """
        theta, sizes, aspects, axial, radial_offset, tilt = self._ring(n_stations)
        local = self.section.fan_origin(sizes, aspects)
        radial, u_hat, v_hat = section_frames(theta, tilt=tilt, sign=self._sign)
        axis = np.zeros_like(radial)
        axis[:, 2] = 1.0
        origin = radial_offset[:, None]*radial + axial[:, None]*axis
        return origin + local[:, :1]*u_hat + local[:, 1:2]*v_hat

    # -- the port ---------------------------------------------------------

    def _port_frame(self):
        """Outward direction and origin of the section carrying the port."""
        geometry = self.geometry
        last = slice(-1, None)
        radial, u_hat, v_hat = section_frames(
            geometry.theta[last], tilt=geometry.tilt[last], sign=self._sign)
        axis = np.array([0.0, 0.0, 1.0])
        origin = (float(geometry.radial[-1])*radial[0]
                  + float(geometry.axial[-1])*axis)
        return origin, u_hat[0], v_hat[0]

    def _dome_shrinkage(self):
        """How fast the dome's depth falls away from the port [m/rad].

        Read only when drawing a stub, which is trimmed against the port
        station's ellipse and has to be sunk past it by whatever the neighbouring
        sections fall short. A dome that grows away from the port leaves nothing
        to make up.
        """
        geometry = self.geometry
        depth, _, _ = self.section.dome(geometry.size[-2:], geometry.aspect[-2:])
        span = float(geometry.theta[-1] - geometry.theta[-2])
        return max(0.0, float(depth[1] - depth[0])/span)

    def _resolve_port(self, port):
        sizes = self.geometry.size[-1:]
        aspects = self.geometry.aspect[-1:]
        origin, u_hat, v_hat = self._port_frame()

        # A family with a dome offers a stub somewhere to sit, aimed from the
        # centre of the ellipse the section is built on -- which is the section
        # origin. One without still takes a tangential duct, which leaves from
        # the bore centre the way a volute's does, so the site simply carries no
        # dome for it to ask after.
        if hasattr(self.section, "dome"):
            depth, half_width, arc_limit = self.section.dome(sizes, aspects)
            local = np.zeros(2)
            dome = DomeSite(
                normal=u_hat, in_plane=v_hat,
                semi_normal=float(depth[0]),
                semi_in_plane=float(half_width[0]),
                arc_limit=float(arc_limit[0]),
                tilt=float(self.geometry.tilt[-1]),
                shrinkage=self._dome_shrinkage(),
            )
        else:
            local = self.section.fan_origin(sizes, aspects)[0]
            dome = None
        start = origin + local[0]*u_hat + local[1]*v_hat

        ring = self.section.polygon(sizes, n_points=96, aspect=aspects)[0]
        return port.resolve(PortSite(
            start=start,
            outlet_diameter=math.sqrt(4.0*float(self.geometry.area[-1])/math.pi),
            max_section_radius=float(
                np.linalg.norm(ring - local, axis=-1).max()),
            sign=self._sign,
            dome=dome,
        ))

    def sweep_port(self, n_points=64, n_stations=9):
        """The port duct swept into machine coordinates.

        A radial port is a plain stub, so this is the tube it cuts, trimmed back
        onto the dome it leaves rather than laid flat against it. What comes back
        is the cavity a printed part would leave for the duct, not a fairing into
        the passage.

        Returns
        -------
        np.ndarray
            Closed boundary rings [m, shape (n_stations, n_points, 3)].
        """
        if self.port_geometry is None:
            raise ValueError("this manifold has no port")
        return self.port_geometry.sweep(n_points, n_stations)

    # -- reporting --------------------------------------------------------

    def collected_velocity(self, mass_flow, density):
        """Through-flow velocity at each station of a branch [m/s].

        Only the share the branch has gathered by that station is passing
        through it, so this is what the area there actually has to carry. A
        uniform ring ramps from nothing to its maximum at the port; a
        constant-velocity ring holds this flat, which is the whole point of it.
        """
        mass_flow = float(mass_flow)
        density = float(density)
        if not np.isfinite(mass_flow) or mass_flow <= 0.0:
            raise ValueError("mass_flow must be finite and positive")
        if not np.isfinite(density) or density <= 0.0:
            raise ValueError("density must be finite and positive")
        return self.geometry.collected*mass_flow/(density*self.geometry.area)

    def __repr__(self):
        sense = "" if self.rotation is None else f", rotation={self.rotation!r}"
        return (f"Manifold({type(self.section).__name__}, "
                f"{self.area_law!r}, branches={self.branches}{sense})")
