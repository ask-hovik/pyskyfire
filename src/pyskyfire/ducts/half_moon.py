"""The half-moon manifold section: a chamber wall on one side, an ellipse on the other.

A manifold wrapped around a thrust chamber is not a free shape. One side of it
is the chamber itself -- whatever the wall happens to do at that station -- and
the other is a dome the designer picks. So the section is an arbitrary inner
curve, a half-ellipse outboard of it, and a fillet at each end blending the two.
Real hardware of this kind reads as a D at the cylindrical part of a chamber and
as a kidney where the wall curves.

Measured along the wall
-----------------------
The slit and the width are arc lengths along the wall, not axial extents. The
slit is space reserved on the wall surface for cooling channels to break out
into, and on a sloped wall the surface runs further than the station span does:
a nozzle at 28 degrees gives only 12.4 mm of wall for every 14 mm of ``x``.
Reserving it axially would quietly hand back a tenth of the room the channels
need. On a cylindrical wall the two measures coincide.

The slit is a reservation, not modelled geometry. Where the channels actually
break through, how they are divided and what the fairing looks like is CAD
work; what pyskyfire draws is the cavity they open into, with enough of the wall
kept clear that the breakout has somewhere to go.

How it is placed
----------------
An anchor station fixes the slit, and it stays there at every azimuth --
channels come out where they come out. ``slit_placement`` says where the slit
sits about that anchor, and ``exposure_placement`` says how the wall left over
once the slit is accounted for is divided either side of it. That division is
held at every azimuth too, so as the area law closes the section down the
exposed wall shrinks onto the slit from both sides in a fixed ratio, and reaches
it exactly at :attr:`minimum_size`.

How it is sized
---------------
The same three regimes a volute works through, for the same reason. ``size`` is
the exposed arc of wall, and shrinking it is the first and cheapest way to take
area out. Once it has closed onto the slit there is no wall left to give up, so
the ellipse's aspect ratio takes over and flattens the dome instead. When that
reaches ``max_aspect_ratio`` the family has run out of knobs, and the area law
gets clipped to the smallest section it can build.

Closing the ends
----------------
The fillets are circular arcs of the radius given -- a real corner radius, the
kind a print or a cutter leaves -- tangent to the wall at each end of the
exposed arc and tangent to the ellipse. Where they touch the ellipse is solved
rather than assumed, because the offset curve of an ellipse is not an ellipse.
The ellipse is carried on the chord between the two ends of the exposed wall and
is centred on it, which is what closes the construction: the chord fixes the
frame, one tangency solve fixes the ellipse, and both fillets follow from it.

Working on the chord means the fillets are tangent to the chord rather than to
the wall's own slope at the point they meet it. The two differ by how far the
wall turns across half a section, which for a manifold on a chamber is a couple
of degrees at most, and the inner boundary still follows the wall exactly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator

from pyskyfire.common.section_geometry import SectionArc, SectionEllipseArc

from .volute_section import VoluteSection

__all__ = ["InnerProfile", "HalfMoonSection"]

#: Bisection steps for the ellipse tangency. The residual is smooth and
#: monotone on its bracket, so this is comfortably past double precision.
_TANGENCY_STEPS = 60


@dataclass(frozen=True)
class InnerProfile:
    """A wall curve in the meridional plane, measured by arc length.

    Attributes
    ----------
    arc : np.ndarray
        Arc length from the anchor [m, shape (M,)], strictly increasing and
        spanning zero. Negative upstream of the anchor, positive downstream.
    points : np.ndarray
        Matching ``(axial, radial)`` positions [m, shape (M, 2)].

    Notes
    -----
    Build one with :meth:`~pyskyfire.regen.outer_wall.OuterWall.profile` to take
    the wall off a resolved thrust chamber, or directly to state a wall
    pyskyfire did not resolve.
    """

    arc: np.ndarray
    points: np.ndarray

    def __post_init__(self):
        arc = np.asarray(self.arc, dtype=float)
        points = np.asarray(self.points, dtype=float)
        if arc.ndim != 1 or len(arc) < 4:
            raise ValueError("arc must be one-dimensional with at least four samples")
        if points.shape != (len(arc), 2):
            raise ValueError("points must have shape (len(arc), 2)")
        if not np.isfinite(arc).all() or not np.isfinite(points).all():
            raise ValueError("the profile must be finite throughout")
        if np.any(np.diff(arc) <= 0.0):
            raise ValueError("arc must be strictly increasing")
        if not arc[0] < 0.0 < arc[-1]:
            raise ValueError(
                "the profile must reach either side of its anchor, so that a "
                "section placed there has wall to grow onto")
        object.__setattr__(self, "arc", arc)
        object.__setattr__(self, "points", points)
        object.__setattr__(self, "_spline", PchipInterpolator(arc, points))

    @property
    def span(self):
        """Arc length available upstream and downstream of the anchor [m]."""
        return float(self.arc[0]), float(self.arc[-1])

    def at(self, s):
        """Meridional position at arc length ``s`` [m, shape (..., 2)]."""
        return self._spline(s)

    def tangent(self, s):
        """Unit tangent at arc length ``s``, pointing downstream."""
        slope = self._spline.derivative()(s)
        return slope/np.linalg.norm(slope, axis=-1, keepdims=True)


@dataclass(frozen=True)
class _WallCurve:
    """The exposed wall, as a boundary piece in each section's own frame.

    Sampled straight off the profile at whatever parameter the quadrature asks
    for, rather than through a polyline, so the enclosed area carries no
    faceting error. Traversed downstream-to-upstream, which is the direction
    that leaves the closed boundary counterclockwise.
    """

    profile: InnerProfile
    lower: np.ndarray
    upper: np.ndarray
    origin: np.ndarray
    axes: np.ndarray

    def _arc(self, fraction):
        fraction = np.asarray(fraction, dtype=float)
        return (self.upper[:, None]
                + fraction[None, :]*(self.lower - self.upper)[:, None])

    def _project(self, points):
        return np.einsum("nij,nkj->nki", self.axes, points)

    def evaluate(self, fraction):
        arc = self._arc(fraction)
        return self._project(self.profile.at(arc) - self.origin[:, None, :])

    def derivative(self, fraction):
        arc = self._arc(fraction)
        sweep = (self.lower - self.upper)[:, None, None]
        return sweep*self._project(self.profile.tangent(arc))


class HalfMoonSection(VoluteSection):
    """A manifold section: chamber wall inboard, half-ellipse outboard.

    Parameters
    ----------
    profile : InnerProfile
        The chamber wall this section sits against, anchored where the slit is.
    slit_height : float
        Arc of wall reserved for the cooling channels to break out into [m].
        The section may never close below it.
    fillet_radius : float
        Corner radius blending each end of the exposed wall into the ellipse
        [m]. Zero leaves a sharp corner.
    slit_placement : float, optional
        Where the slit sits about the anchor, as the fraction of it lying
        upstream. ``0.5`` centres it.
    exposure_placement : float, optional
        How the exposed wall outside the slit divides, as the fraction of it
        lying upstream of the slit. Held at every azimuth, so ``0.75`` keeps
        three quarters of the spare wall upstream all the way round.
    aspect_ratio, max_aspect_ratio : float, optional
        Ellipse width over depth at the design point, and the flattest it may
        be pushed once the exposed wall has closed onto the slit.
    """

    #: The wall is drawn first, and it is the face the channels break out
    #: through rather than a wall of the cavity.
    opening_index = 0

    def __init__(self, profile, slit_height, fillet_radius, *,
                 slit_placement=0.5, exposure_placement=0.5,
                 aspect_ratio=1.0, max_aspect_ratio=None):
        if not isinstance(profile, InnerProfile):
            raise TypeError("profile must be an InnerProfile")
        super().__init__(slit_height, aspect_ratio, max_aspect_ratio)
        fillet_radius = float(fillet_radius)
        if not np.isfinite(fillet_radius) or fillet_radius < 0.0:
            raise ValueError("fillet_radius must be finite and nonnegative")
        for name, value in (("slit_placement", slit_placement),
                            ("exposure_placement", exposure_placement)):
            value = float(value)
            if not np.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must lie between zero and one")
        self.profile = profile
        self.fillet_radius = fillet_radius
        self.slit_placement = float(slit_placement)
        self.exposure_placement = float(exposure_placement)

    @property
    def minimum_size(self) -> float:
        """The slit itself, where the exposed wall has nothing left to give."""
        return self.inlet_width

    @property
    def maximum_size(self) -> float:
        """The most wall there is to spread onto, given where the slit sits.

        The window opens either side of the fixed slit in the fixed ratio, so
        whichever side runs out of wall first is what caps the section.
        """
        low, high = self.slit
        reach_low, reach_high = self.profile.span
        spare = []
        if self.exposure_placement > 0.0:
            spare.append((low - reach_low)/self.exposure_placement)
        if self.exposure_placement < 1.0:
            spare.append((reach_high - high)/(1.0 - self.exposure_placement))
        return self.inlet_width + min(spare)

    # -- placement --------------------------------------------------------

    @property
    def slit(self):
        """Arc-length bounds of the reserved slit [m], about the anchor."""
        width = self.inlet_width
        return -self.slit_placement*width, (1.0 - self.slit_placement)*width

    def window(self, size):
        """Arc-length bounds of the exposed wall at each size [m].

        The slit stays put and the spare wall divides about it in the fixed
        ratio, so the window closes onto the slit as ``size`` reaches
        :attr:`minimum_size`.
        """
        sizes, _ = self._sizes(size)
        if np.any(sizes < self.inlet_width - 1.0e-15):
            raise ValueError(
                f"size must be at least the slit height {self.inlet_width:g} m, "
                "which the section may never close over")
        low, high = self.slit
        spare = sizes - self.inlet_width
        lower = low - self.exposure_placement*spare
        upper = high + (1.0 - self.exposure_placement)*spare
        reach_low, reach_high = self.profile.span
        if np.any(lower < reach_low) or np.any(upper > reach_high):
            raise ValueError(
                "the section reaches past the wall it was given; extend the "
                "InnerProfile or shrink the manifold")
        return lower, upper

    def _frame(self, size):
        """Chord frame of the exposed wall: origin, axes, and chord length."""
        lower, upper = self.window(size)
        ends = self.profile.at(np.stack((lower, upper), axis=-1))
        origin = ends.mean(axis=1)
        along = ends[:, 1] - ends[:, 0]
        chord = np.linalg.norm(along, axis=-1)
        along = along/chord[:, None]
        # Outward is the wall's own outward side: the chamber axis is +x and
        # radius is the second coordinate, so rotating the along-wall direction
        # by a quarter turn one way points away from the axis.
        outward = np.column_stack((-along[:, 1], along[:, 0]))
        axes = np.stack((outward, along), axis=1)
        return lower, upper, origin, axes, chord

    # -- the ellipse the fillets land on ----------------------------------

    def _tangency(self, chord, aspects):
        """Ellipse parameter where each fillet touches it.

        The fillet centres sit one radius inboard of the chord at each of its
        ends, so requiring the ellipse to be tangent to both -- and centred on
        the chord -- leaves a single equation. It is monotone on ``(0, pi/2)``,
        running to plus infinity as the parameter approaches the chord and to
        ``-fillet_radius`` at the outboard extreme, so the bracket is fixed.
        """
        radius = self.fillet_radius
        low = np.full_like(chord, 1.0e-12)
        high = np.full_like(chord, np.pi/2.0 - 1.0e-12)

        def residual(t):
            sin, cos = np.sin(t), np.cos(t)
            scale = np.sqrt(aspects**2*cos**2 + sin**2)
            return (chord*cos/sin/(2.0*aspects)
                    + radius*cos*(1.0/aspects - aspects)/scale
                    - radius)

        for _ in range(_TANGENCY_STEPS):
            middle = 0.5*(low + high)
            positive = residual(middle) > 0.0
            low = np.where(positive, middle, low)
            high = np.where(positive, high, middle)
        return 0.5*(low + high)

    def _layout(self, size, aspect=None):
        """Resolve the chord frame, the ellipse, and both fillet arcs."""
        sizes, aspects = self._sizes(size, aspect)
        lower, upper, origin, axes, chord = self._frame(sizes)
        tau = self._tangency(chord, aspects)
        scale = np.sqrt(aspects**2*np.cos(tau)**2 + np.sin(tau)**2)
        semi_v = chord/(2.0*np.sin(tau)) + self.fillet_radius/scale
        semi_u = semi_v/aspects
        # Where each fillet centre sees its own tangency point.
        phi = np.arctan2(np.sin(tau), aspects*np.cos(tau))
        return (lower, upper, origin, axes, chord, semi_u, semi_v, tau, phi)

    # -- the boundary -----------------------------------------------------

    def boundary(self, size, aspect=None):
        lower, upper, origin, axes, chord, semi_u, semi_v, tau, phi = \
            self._layout(size, aspect)
        half = chord/2.0
        zero = np.zeros_like(chord)
        radii = np.full_like(chord, self.fillet_radius)
        centre = np.full_like(chord, self.fillet_radius)
        pi = np.full_like(chord, np.pi)
        return (
            # The wall itself, downstream end back to upstream end.
            _WallCurve(self.profile, lower, upper, origin, axes),
            # Upstream fillet, off the wall and round onto the ellipse.
            SectionArc(np.column_stack((centre, -half)), radii,
                       pi, 2.0*np.pi - phi),
            # The dome, the long way round past its outboard extreme.
            SectionEllipseArc(np.column_stack((zero, zero)), semi_u, semi_v,
                              -tau, tau),
            # Downstream fillet, back off the ellipse onto the wall.
            SectionArc(np.column_stack((centre, half)), radii, phi, pi),
        )

    def fan_origin(self, size, aspect=None):
        """A point on the centreline that sees the whole boundary.

        Outboard of the wall's own bulge, since a wall curving into the cavity
        would otherwise hide part of itself from a point on the chord.
        """
        lower, upper, origin, axes, chord, semi_u, semi_v, _, _ = \
            self._layout(size, aspect)
        wall = _WallCurve(self.profile, lower, upper, origin, axes)
        crest = wall.evaluate(np.linspace(0.0, 1.0, 33))[..., 0].max(axis=1)
        return np.column_stack((0.5*(np.maximum(crest, 0.0) + semi_u),
                                np.zeros_like(chord)))

    def slit_boundary(self, size, aspect=None):
        """The reserved slit alone, as a boundary piece in the section frame.

        The same wall curve the section is drawn with, read only between the
        slit's own bounds. What comes back is the strip the cooling channels
        break out through, which is otherwise indistinguishable from the rest
        of the chamber-facing face.
        """
        lower, upper, origin, axes, _ = self._frame(self._sizes(size, aspect)[0])
        low, high = self.slit
        return _WallCurve(self.profile, np.full_like(lower, low),
                          np.full_like(upper, high), origin, axes)

    def margin_boundaries(self, size, aspect=None):
        """The chamber-facing wall either side of the slit, as boundary pieces.

        Returned in the order the face is traversed, so the two of them and
        :meth:`slit_boundary` between them tile it exactly. That matters for
        drawing: the slit is where channels break out and the rest of the face
        is metal against the chamber, and the two want telling apart without
        laying one surface on top of the other.
        """
        lower, upper, origin, axes, _ = self._frame(self._sizes(size, aspect)[0])
        low, high = self.slit
        return (_WallCurve(self.profile, np.full_like(lower, high), upper,
                           origin, axes),
                _WallCurve(self.profile, lower, np.full_like(upper, low),
                           origin, axes))

    def dome(self, size, aspect=None):
        """The outboard dome, as ``(depth, half_width, arc_limit)``.

        The depth is how far the section reaches out along the normal from its
        chord and the half-width is the ellipse axis across it, both semi-axes
        of the ellipse the section is built on. The arc limit is the ellipse
        parameter the fillets take over at: only within it is the ellipse
        passage boundary, so it is what bounds where a port may leave through.
        """
        *_, semi_u, semi_v, tau, _ = self._layout(size, aspect)
        return semi_u, semi_v, tau

    def meridional_frame(self, size, aspect=None):
        """Where each section sits on the wall, as ``(axial, radial, tilt)``.

        A manifold's section slides along the wall as the area law resizes it,
        and leans with whatever the wall is doing there, so it carries its own
        placement rather than taking one from the passage. ``tilt`` is the angle
        of the section's outward direction from radial, positive toward ``+x``,
        which is what a sweep needs to lay the section plane down correctly.
        """
        _, _, origin, axes, _ = self._frame(self._sizes(size, aspect)[0])
        outward = axes[:, 0]
        return origin[:, 0], origin[:, 1], np.arctan2(outward[:, 0], outward[:, 1])

    def exposed_wall(self, size, aspect=None):
        """Arc of wall left clear between the two fillets [m].

        Equal to ``size``: the fillets are tangent to the chord at the ends of
        the exposed arc, so they begin where it ends rather than eating into
        it. Reported so that a design can be checked against the room the
        channel breakout actually needs.
        """
        lower, upper = self.window(size)
        return upper - lower
