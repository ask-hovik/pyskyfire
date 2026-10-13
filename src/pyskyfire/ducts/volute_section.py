"""Volute cross-section families.

A volute cross-section is a closed planar region measured in a section plane
that contains the machine axis. Its two coordinates are

``u``
    Distance outward from the volute inlet plane, so ``u = 0`` sits at the
    inlet radius and the passage grows toward larger ``u``.
``v``
    Distance across the inlet slot, zero on the slot centreline. For an
    untilted volute this is the axial direction; a tilted inlet rotates the
    whole plane rigidly, which leaves every quantity computed here unchanged.

Sections are a one-parameter family. The inlet slot width is fixed by the
hardware feeding the volute -- a pump discharge, or the channel-exit band of a
regenerative circuit -- so a section grows only through a single ``size``
parameter, and the area law is what picks ``size`` at each azimuth.

Every section reports the three quantities the design and loss models need:
enclosed area, the radial position of the area centroid, and the boundary
itself. They come from one Green's-theorem integral over the exact arcs and
lines, so none of them carry a polygon-resolution error.

These classes deliberately do not subclass the cooling-channel
:class:`~pyskyfire.regen.cross_section.ChannelSection`. That class carries
thermal exposure fractions, fin conduction, wall thicknesses and sector
placement, none of which a volute has, and it is sized by a different
constraint. Only the planar primitives in
:mod:`pyskyfire.common.section_geometry` are shared.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from scipy.optimize import brentq

from pyskyfire.common.section_geometry import SectionArc, SectionEllipseArc, SectionLine

__all__ = [
    "VoluteSection",
    "EllipticalSlottedSection",
    "section_frames",
    "solve_section_size",
]


# Gauss-Legendre nodes on [0, 1], matching the order used by the shared
# section-geometry metrics so both agree to round-off on common shapes.
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(24)
_NODES = (_NODES + 1.0) / 2.0
_WEIGHTS = _WEIGHTS / 2.0


def _closed_boundary_metrics(pieces):
    """Area, centroid and per-piece length of a closed planar boundary.

    ``pieces`` are :class:`SectionArc` or :class:`SectionLine` objects in
    counterclockwise order, each batched over the same number of sections.
    Green's theorem is integrated along the exact curves rather than a sampled
    polygon.

    Returns
    -------
    area : np.ndarray
        Enclosed area of each section [m^2, shape (N,)].
    centroid : np.ndarray
        Area centroid in ``(u, v)`` [m, shape (N, 2)].
    lengths : np.ndarray
        Arc length of each piece [m, shape (N, len(pieces))].
    """
    area = None
    moment = None
    lengths = []
    for piece in pieces:
        points = piece.evaluate(_NODES)
        derivatives = piece.derivative(_NODES)
        if area is None:
            area = np.zeros(points.shape[0])
            moment = np.zeros((points.shape[0], 2))
        u, v = points[..., 0], points[..., 1]
        du, dv = derivatives[..., 0], derivatives[..., 1]
        area += ((u * dv - v * du) @ _WEIGHTS) / 2.0
        moment[:, 0] += (u * u * dv) @ _WEIGHTS / 2.0
        moment[:, 1] -= (v * v * du) @ _WEIGHTS / 2.0
        lengths.append(np.linalg.norm(derivatives, axis=-1) @ _WEIGHTS)
    if not np.all(np.isfinite(area)) or np.any(area <= 0.0):
        raise ValueError(
            "Volute boundary must enclose a finite positive area; a "
            "non-positive result means the pieces are not a simple "
            "counterclockwise closed curve"
        )
    return area, moment / area[:, None], np.stack(lengths, axis=1)


def solve_section_size(section, residual):
    """The smallest section of a family satisfying an area law.

    Works down three regimes as the law asks for less and less area, and is
    shared by every passage that sizes a section this way -- a volute closing
    on its tongue and a manifold closing on its dead point face the same
    problem, because both laws want zero area somewhere no section can deliver.

    Away from that point the design aspect ratio can hold the law, and only the
    size is solved for. Closer in, the law wants less area than even the
    smallest section at that aspect: the passage can never be narrower than the
    opening feeding it, so the size has bottomed out. Flattening is the way past
    that -- holding the size at its floor and stretching the section across the
    opening shrinks its area without narrowing the passage -- so the aspect
    ratio is solved for instead, up to the family's limit. Only when the
    flattest allowed section is still too big does the law get clipped, and it
    is clipped to that flattest section rather than to the design-aspect one,
    which is as close to the law as the family can get.

    Parameters
    ----------
    section : VoluteSection
        The family being sized.
    residual : callable
        ``residual(size, aspect)``, negative where the trial section is too
        small and increasing with size.

    Returns
    -------
    size, aspect, clipped
    """
    design = section.aspect_ratio
    ceiling = section.max_aspect_ratio
    minimum = section.minimum_size

    if residual(minimum, design) < 0.0:
        ceiling_size = section.maximum_size
        high = min(max(minimum*2.0, minimum + section.inlet_width), ceiling_size)
        for _ in range(60):
            if residual(high, design) > 0.0:
                break
            if high >= ceiling_size:
                raise RuntimeError(
                    f"the area law wants more area than this section can build. "
                    f"At its largest, {ceiling_size:g} m, it still falls short; "
                    "give the section more room to grow onto, or ask the law "
                    "for less area")
            high = min(high*2.0, ceiling_size)
        else:
            raise RuntimeError(
                "could not bracket a section size; the area law does not "
                "appear to be satisfiable by this section")
        size = brentq(residual, minimum, high, args=(design,),
                      xtol=1e-14, rtol=1e-12)
        return size, design, False

    if ceiling > design and residual(minimum, ceiling) < 0.0:
        aspect = brentq(lambda value: residual(minimum, value), design, ceiling,
                        xtol=1e-12, rtol=1e-12)
        return minimum, aspect, False

    return minimum, ceiling, True


def section_frames(theta, *, tilt=0.0, sign=1.0):
    """Bases of the section planes at one or more azimuths about ``+z``.

    A section lies in the plane through the machine axis at its azimuth, so the
    ``(u, v)`` coordinates it is measured in map to a radial and an axial
    direction there, rotated together by ``tilt``.

    ``sign`` is the rotation sense from
    :func:`~pyskyfire.common.rotation.rotation_sign`: ``+1`` winds
    counterclockwise about ``+z`` as azimuth increases, ``-1`` clockwise. It
    mirrors where each section sits without touching the section's own shape,
    since the mirror plane contains the axis and so maps every section plane to
    another one.

    ``tilt`` stays absolute under that mirror: which axial way a discharge
    points is a packaging choice independent of which way the machine spins.

    Returns
    -------
    radial, u_hat, v_hat : np.ndarray
        Unit vectors [shape (N, 3)], where ``radial`` locates the inlet and the
        other two span the section plane.
    """
    theta = np.atleast_1d(np.asarray(theta, dtype=float))
    radial = np.column_stack(
        (np.cos(theta), sign * np.sin(theta), np.zeros_like(theta)))
    axial = np.zeros_like(radial)
    axial[:, 2] = 1.0
    # Tilt may vary station to station: a manifold's section leans with the
    # wall it sits on, where a volute's holds one inclination all the way round.
    tilt = np.broadcast_to(np.atleast_1d(np.asarray(tilt, dtype=float)),
                           theta.shape)[:, None]
    u_hat = np.cos(tilt) * radial + np.sin(tilt) * axial
    v_hat = -np.sin(tilt) * radial + np.cos(tilt) * axial
    return radial, u_hat, v_hat


class VoluteSection(ABC):
    """A volute cross-section family with a fixed inlet slot width.

    Subclasses define a shape parameterised by a positive ``size`` and an
    ``aspect``. All methods accept a scalar or a one-dimensional array of
    sizes, and return matching arrays; ``aspect`` defaults to the family's own
    :attr:`aspect_ratio` where it is not given.

    Attributes
    ----------
    inlet_width : float
        Width of the inlet slot ``b5`` [m]. Fixed by the feeding hardware and
        held constant as the section grows, which is the constraint that makes
        the family one-parameter in ``size``.
    aspect_ratio : float
        The shape the family takes at its design point.
    max_aspect_ratio : float
        The flattest the family may be pushed when the area law asks for less
        area than ``size`` alone can deliver. Equal to :attr:`aspect_ratio` for
        a family with no second degree of freedom, which then simply has a
        floor on the area it can produce.
    """

    #: Which boundary piece is the opening the flow enters or leaves through,
    #: rather than a wall. A volute section is closed by its slot mouth, so the
    #: opening is the last piece; a manifold section leads with the chamber
    #: wall it reserves, so that family overrides this.
    opening_index = -1

    def __init__(self, inlet_width, aspect_ratio=1.0, max_aspect_ratio=None):
        inlet_width = float(inlet_width)
        if not np.isfinite(inlet_width) or inlet_width <= 0.0:
            raise ValueError("inlet_width must be finite and positive")
        aspect_ratio = float(aspect_ratio)
        if not np.isfinite(aspect_ratio) or aspect_ratio <= 0.0:
            raise ValueError("aspect_ratio must be finite and positive")
        max_aspect_ratio = (aspect_ratio if max_aspect_ratio is None
                            else float(max_aspect_ratio))
        if not np.isfinite(max_aspect_ratio) or max_aspect_ratio < aspect_ratio:
            raise ValueError(
                "max_aspect_ratio must be finite and at least aspect_ratio; "
                "flattening only ever makes a section thinner than its design "
                f"shape, so {max_aspect_ratio:g} leaves no room below "
                f"{aspect_ratio:g}")
        self.inlet_width = inlet_width
        self.aspect_ratio = aspect_ratio
        self.max_aspect_ratio = max_aspect_ratio

    @property
    @abstractmethod
    def minimum_size(self) -> float:
        """Smallest ``size`` that still yields a valid section [m]."""

    @property
    def maximum_size(self) -> float:
        """Largest ``size`` the family can build [m].

        Unbounded for a family that can simply keep growing. A section measured
        against something finite -- a manifold's section has only so much
        chamber wall to spread onto -- reports that limit here, so a law asking
        for more area than exists is told so rather than being chased past the
        end of the geometry.
        """
        return float("inf")

    @abstractmethod
    def boundary(self, size, aspect=None):
        """Counterclockwise boundary pieces for each size and aspect."""

    def _sizes(self, size, aspect=None):
        """Broadcast ``size`` and ``aspect`` into a matching pair of arrays."""
        sizes = np.atleast_1d(np.asarray(size, dtype=float))
        if sizes.ndim != 1:
            raise ValueError("size must be a scalar or a one-dimensional array")
        if not np.all(np.isfinite(sizes)) or np.any(sizes <= 0.0):
            raise ValueError("size must be finite and positive")
        if aspect is None:
            aspect = self.aspect_ratio
        aspects = np.broadcast_to(
            np.atleast_1d(np.asarray(aspect, dtype=float)), sizes.shape
        ).astype(float, copy=False)
        if not np.all(np.isfinite(aspects)) or np.any(aspects <= 0.0):
            raise ValueError("aspect must be finite and positive")
        return sizes, aspects

    def metrics(self, size, aspect=None):
        """Return ``(area, centroid, lengths)`` for each size and aspect."""
        return _closed_boundary_metrics(self.boundary(size, aspect))

    def area(self, size, aspect=None):
        """Enclosed passage area [m^2]."""
        return self.metrics(size, aspect)[0]

    def centroid_offset(self, size, aspect=None):
        """Area-centroid distance outboard of the inlet plane [m].

        The volute mean radius ``r_c`` of Aungier's design equation is the
        inlet radius plus this offset.
        """
        return self.metrics(size, aspect)[1][:, 0]

    def perimeter(self, size, aspect=None, *, wetted=True):
        """Boundary length [m].

        With ``wetted`` set, the inlet slot mouth is excluded, since it is an
        opening into the feeding passage rather than a wall.
        """
        lengths = self.metrics(size, aspect)[2]
        return lengths[:, :-1].sum(axis=1) if wetted else lengths.sum(axis=1)

    def fan_origin(self, size, aspect=None):
        """A point each section is star-shaped about, in ``(u, v)`` [m].

        Every boundary point is visible from it, so a triangle fan to the
        boundary is a valid triangulation of the section. Meshing uses this to
        cap the ends of a swept passage. The default is the area centroid,
        which is correct for a convex family; a section with a reentrant
        boundary must override it.

        Returns
        -------
        np.ndarray
            One origin per size [m, shape (N, 2)].
        """
        centroid = self.metrics(size, aspect)[1]
        return centroid

    def radial_profile(self, size, angles, aspect=None, *, samples=720):
        """Distance from :meth:`fan_origin` to the boundary at each angle [m].

        The section is star-shaped about that origin, so the boundary is a
        single-valued ``r(psi)`` about it. That is what lets an exit cone morph
        a section into a circle: blend this radius toward a constant one.

        Parameters
        ----------
        size : float or array_like
            Section size, one per section.
        angles : array_like
            Angles about the fan origin [rad], measured from the ``+u`` axis.
            Taken modulo a full turn.
        samples : int, optional
            Boundary samples per piece used to build the profile.

        Returns
        -------
        np.ndarray
            Radii [m, shape (N, len(angles))].
        """
        angles = np.asarray(angles, dtype=float)
        if angles.ndim != 1:
            raise ValueError("angles must be one-dimensional")
        # Drop the repeated closing point; the periodic interpolation re-closes.
        ring = self.polygon(size, n_points=int(samples), aspect=aspect)[:, :-1]
        origin = self.fan_origin(size, aspect)
        offset = ring - origin[:, None, :]
        psi = np.arctan2(offset[..., 1], offset[..., 0])
        radius = np.linalg.norm(offset, axis=-1)
        profile = np.empty((ring.shape[0], len(angles)))
        for index in range(ring.shape[0]):
            # Star-shapedness makes psi monotonic along a counterclockwise
            # traversal, so unwrapping gives exactly one increasing period.
            turned = np.unwrap(psi[index])
            keep = np.concatenate(([True], np.diff(turned) > 0.0))
            profile[index] = np.interp(angles, turned[keep], radius[index][keep],
                                       period=2.0 * np.pi)
        return profile

    def polygon(self, size, n_points=64, aspect=None):
        """Sampled closed boundary for plotting.

        Returns
        -------
        np.ndarray
            Boundary points in ``(u, v)`` [m, shape (N, M, 2)], one closed ring
            per size.
        """
        n_points = int(n_points)
        if n_points < 2:
            raise ValueError("n_points must be at least two")
        pieces = self.boundary(size, aspect)
        rings = [piece.evaluate(np.linspace(0.0, 1.0, n_points)) for piece in pieces]
        ring = np.concatenate(rings, axis=1)
        return np.concatenate([ring, ring[:, :1]], axis=1)



class EllipticalSlottedSection(VoluteSection):
    """An elliptical bore opening into a parallel-sided inlet slot.

    This is the cast pump-volute shape: an elliptical passage fed through a
    slot of fixed width ``inlet_width``, with the junction between each slot
    wall and the bore rounded by a fillet of user-set radius. At the default
    aspect ratio of one the ellipse is a circle, which is the common case.

    Parameters
    ----------
    inlet_width : float
        Slot width ``b5`` [m].
    fillet_radius : float
        Radius of the two wall-to-bore fillets [m]. Zero gives sharp corners.
    throat_length : float, optional
        Length of straight parallel slot wall between the inlet plane and the
        start of each fillet [m]. Defaults to zero, which starts the fillet at
        the inlet plane.
    aspect_ratio : float, optional
        Bore semi-axis across the slot divided by its semi-axis outward from
        the inlet plane. One is a circle. Above one the bore is stretched
        across the slot and squashed radially; below one it reaches further
        out than it is wide. Defaults to one.
    max_aspect_ratio : float, optional
        How far the aspect ratio may be pushed above its design value where
        the area law asks for less area than ``size`` alone can give -- see
        Notes. Defaults to four.

    Notes
    -----
    ``size`` is the semi-axis across the slot, so it plays exactly the role the
    bore radius does for a circle and the slot still fixes its floor at
    ``inlet_width / 2``. Writing ``c`` for it, ``a = c / aspect`` for the
    semi-axis outward from the inlet plane, ``b`` for the inlet width, ``r_f``
    for the fillet radius and ``L`` for the throat length, each fillet is
    placed by requiring it to be tangent to its slot wall and to the bore.

    The fillets round a *reflex* corner -- the passage opens outward as it
    leaves the slot -- so each fillet centre sits outside the passage at
    ``v = +/- (b/2 + r_f)`` and touches the bore externally. Writing the
    tangency point as the ellipse parameter ``t`` and the ellipse's outward
    unit normal there as ``n(t)``, the fillet centre is ``E(t) + r_f*n(t)``,
    and requiring its ``v`` to equal ``-(b/2 + r_f)`` fixes ``t`` alone::

        c*sin(t) + r_f*sin(t) / (c*N(t)) + b/2 + r_f = 0
        N(t) = sqrt(cos(t)**2/a**2 + sin(t)**2/c**2)

    which is solved by bisection on ``t`` in ``(-pi, -pi/2)``, where the left
    end is positive and the right end is ``b/2 - c <= 0``. The bore centre then
    follows explicitly. For a circle this reduces exactly to the closed form
    ``u_c = L + sqrt((c + r_f)**2 - (b/2 + r_f)**2)``.

    Requiring ``size > b/2`` -- a bore at least as wide as the slot it is fed
    through -- is enough for the construction to be well posed at any aspect
    ratio.

    A shaped volute wants a vanishing area at the tongue, and no section whose
    bore must span the slot can shrink that far. The aspect ratio is the way
    out: holding ``size`` at its floor and flattening the ellipse shrinks the
    area roughly as ``1 / aspect`` without ever narrowing the bore below the
    slot. :class:`~pyskyfire.ducts.Volute` uses that automatically, so the
    passage near the tongue tracks the area law further than a circle could
    before it has to give up and clip. ``max_aspect_ratio`` bounds how thin a
    sliver that is allowed to produce.

    The bore reaches inboard of the inlet plane when its centre sits closer to
    that plane than its own outward semi-axis. Such a section is an
    intermediate rather than a fully external volute, and it is reported
    through :meth:`inboard_overhang` rather than rejected, since it is valid
    geometry that may simply foul neighbouring hardware. Raising
    ``throat_length`` pushes the bore back outboard.
    """

    #: Bisection steps used to place the fillet tangency. Each halves a bracket
    #: of pi/2, so this reaches the double-precision floor on the resulting
    #: point positions. Bisection rather than anything faster-converging on
    #: purpose: as the bore closes on the slot the root slides onto the end of
    #: the bracket and meets it tangentially, where Newton stalls and a
    #: fixed-point iteration on the same residual diverges outright.
    _TANGENCY_STEPS = 52

    def __init__(self, inlet_width, fillet_radius, throat_length=0.0,
                 aspect_ratio=1.0, max_aspect_ratio=4.0):
        super().__init__(inlet_width, aspect_ratio, max_aspect_ratio)
        fillet_radius = float(fillet_radius)
        throat_length = float(throat_length)
        if not np.isfinite(fillet_radius) or fillet_radius < 0.0:
            raise ValueError("fillet_radius must be finite and nonnegative")
        if not np.isfinite(throat_length) or throat_length < 0.0:
            raise ValueError("throat_length must be finite and nonnegative")
        self.fillet_radius = fillet_radius
        self.throat_length = throat_length

    @property
    def minimum_size(self) -> float:
        """Half the inlet width, where the bore just spans the slot."""
        return self.inlet_width / 2.0

    def _tangency(self, semi_u, semi_v, h):
        """Ellipse parameter where each fillet touches the bore.

        Solved numerically rather than in closed form because the offset curve
        of an ellipse is not an ellipse. The bracket is fixed and its ends are
        known to straddle the root: at ``t = -pi`` the residual is ``h > 0``,
        and at ``t = -pi/2`` it is ``b/2 - c <= 0``. That second end is an
        equality exactly when the bore has closed down onto the slot, and there
        the residual meets zero tangentially -- a double root, which is why
        this bisects rather than reaching for faster convergence.
        """
        r_f = self.fillet_radius
        inverse_u, inverse_v = 1.0 / semi_u**2, 1.0 / semi_v**2
        low = np.full_like(semi_v, -np.pi)
        high = np.full_like(semi_v, -np.pi / 2.0)

        def residual(t):
            sin, cos = np.sin(t), np.cos(t)
            normal = np.sqrt(cos**2 * inverse_u + sin**2 * inverse_v)
            return semi_v * sin + r_f * sin / (semi_v * normal) + h

        for _ in range(self._TANGENCY_STEPS):
            middle = 0.5 * (low + high)
            positive = residual(middle) > 0.0
            low = np.where(positive, middle, low)
            high = np.where(positive, high, middle)
        return 0.5 * (low + high)

    def _layout(self, size, aspect=None):
        """Resolve bore semi-axes and centre, and both tangency angles."""
        semi_v, aspects = self._sizes(size, aspect)
        if np.any(semi_v < self.minimum_size):
            raise ValueError(
                f"size must be at least inlet_width/2 = {self.minimum_size:g} m "
                "so that the bore is wider than the slot feeding it"
            )
        semi_u = semi_v / aspects
        r_f = self.fillet_radius
        h = self.inlet_width / 2.0 + r_f
        t = self._tangency(semi_u, semi_v, h)
        sin, cos = np.sin(t), np.cos(t)
        normal = np.sqrt(cos**2 / semi_u**2 + sin**2 / semi_v**2)
        normal_u = cos / (semi_u * normal)
        normal_v = sin / (semi_v * normal)
        u_f = np.full_like(semi_v, self.throat_length)
        # The fillet centre is the tangency point pushed out along the normal,
        # and its u is pinned to the throat, so the bore centre follows.
        u_c = u_f - semi_u * cos - r_f * normal_u
        # Angle from a fillet centre to where it touches the bore: the
        # ellipse's inward normal there.
        psi = np.arctan2(-normal_v, -normal_u)
        return semi_u, semi_v, r_f, h, u_f, u_c, t, psi

    def fan_origin(self, size, aspect=None):
        """The bore centre, which sees the whole boundary.

        The area centroid does not, in general. A long throat drags it inboard
        of where the fillets begin, and from inside the slot a line to a fillet
        -- which lies beyond the slot wall it is tangent to -- leaves the
        passage through that wall. The bore centre is always outboard of the
        fillet tangency, so the whole boundary stays in view.
        """
        _, _, _, _, _, u_c, _, _ = self._layout(size, aspect)
        return np.column_stack((u_c, np.zeros_like(u_c)))

    def inboard_overhang(self, size, aspect=None):
        """How far the bore reaches inboard of the inlet plane [m].

        Zero for a fully external volute, positive for an intermediate one.
        """
        semi_u, _, _, _, _, u_c, _, _ = self._layout(size, aspect)
        return np.maximum(semi_u - u_c, 0.0)

    def boundary(self, size, aspect=None):
        semi_u, semi_v, r_f, h, u_f, u_c, t, psi = self._layout(size, aspect)
        half = self.inlet_width / 2.0
        zero = np.zeros_like(semi_v)
        fillet_radii = np.full_like(semi_v, r_f)
        lower_fillet = np.column_stack((u_f, np.full_like(semi_v, -h)))
        upper_fillet = np.column_stack((u_f, np.full_like(semi_v, h)))
        return (
            # Lower slot wall, outward from the mouth.
            SectionLine(np.column_stack((zero, np.full_like(semi_v, -half))),
                        np.column_stack((u_f, np.full_like(semi_v, -half)))),
            # Lower fillet, from the wall around to the bore.
            SectionArc(lower_fillet, fillet_radii,
                       np.full_like(semi_v, np.pi / 2.0), psi),
            # The bore, the long way around past the outboard extreme. The
            # upper tangency mirrors the lower one, at parameter -t.
            SectionEllipseArc(np.column_stack((u_c, zero)), semi_u, semi_v, t, -t),
            # Upper fillet, from the bore back down to the wall.
            SectionArc(upper_fillet, fillet_radii,
                       -psi, np.full_like(semi_v, -np.pi / 2.0)),
            # Upper slot wall, inward to the mouth.
            SectionLine(np.column_stack((u_f, np.full_like(semi_v, half))),
                        np.column_stack((zero, np.full_like(semi_v, half)))),
            # The slot mouth itself, closing the section on the inlet plane.
            SectionLine(np.column_stack((zero, np.full_like(semi_v, half))),
                        np.column_stack((zero, np.full_like(semi_v, -half)))),
        )
