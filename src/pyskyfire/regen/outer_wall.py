"""The cold-side envelope of a cooling circuit, read like a contour.

:class:`~pyskyfire.regen.contour.Contour` describes where the hot gas sees the
chamber. This describes the other face of the same hardware: the outermost
surface the cooling channels and their walls reach. It is what a manifold
wrapped around the chamber has to sit against, and it is stated the same way --
a radius as a function of axial station, plus the direction the wall faces
there -- so the two can be read side by side.

A channelled wall is not a smooth cylinder, so "the outer radius" is two
numbers rather than one.

``envelope``
    The outermost point the hardware reaches. Rounded tubes bulge, so this is
    the crest of a bump. Anything that has to *clear* the chamber -- a manifold
    cavity, a jacket, a clamp -- is sized against this one, which is why it is
    the default.
``equivalent``
    The radius of the plain cylinder enclosing the same cross-sectional area as
    the real scalloped surface, so it sits between the crests and the troughs.
    This is the one that belongs in area, volume and mass bookkeeping, where
    counting out to the crests would claim material that is not there.

The two are equal wherever the outer surface really is a cylinder, which is the
case for a squared channel: its roof is an arc at constant radius and the rib
tops between channels sit at that same radius. They separate only for sections
whose exterior is scalloped.

Both are surfaces of revolution even when the channels are helical. A helical
channel presents a narrower section in the plane normal to its own path, and
that narrowing is already carried in the resolved section; sweeping it about
the axis still traces one radius per station. So the outer wall stays straight
in exactly the sense the contour does.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator

__all__ = ["OuterWall", "OUTER_WALL_KINDS"]

#: The two ways the outer radius of a channelled wall can be meant.
OUTER_WALL_KINDS = ("envelope", "equivalent")

#: Samples along the exterior roof curve used to find its outermost point.
_ROOF_SAMPLES = 513

# Gauss-Legendre nodes on [0, 1], matching the order the shared section metrics
# use. The roof is an arc or a line, so the swept area is exact rather than
# chorded -- a polyline would quietly under-measure every convex roof and stop
# a squared channel's two radii from landing on each other exactly.
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(24)
_NODES = (_NODES + 1.0)/2.0
_WEIGHTS = _WEIGHTS/2.0


def _sector(section, prof):
    """Measure one channel's sector: ``(envelope, swept_area, angular_span)``.

    The exterior roof spans exactly one channel pitch -- a rounded tube's outer
    lobe is tangent to both sector rays, and a squared channel's offset roof
    carries half a rib on each side -- so one roof curve closes one sector
    against the axis. The sector's two bounding rays are radial and sweep no
    area, which leaves the roof alone carrying the whole integral.

    The angular span is returned alongside the area because a circuit does not
    always tile the circle on its own: where two circuits interleave, each owns
    only part of the circumference, and the area has to be shared out by angle
    rather than multiplied up by one circuit's channel count.
    """
    try:
        roof = section.surfaces(prof, exterior=True)["roof"]
    except NotImplementedError:
        raise ValueError(
            f"{type(section).__name__} does not define an exterior boundary, so "
            "it has no outer wall to measure. Give the circuit a section that "
            "closes on its cold side, or supply the outer wall explicitly."
        ) from None

    base = float(prof.cross_section_base_point_path[0, 1])

    # Boundaries already carry the transverse frame -- a squared section builds
    # them there and only maps back through its projection to *measure* the
    # channel in its own tilted plane. Radii about the engine axis want the
    # transverse frame, so the curves are read as they come. The check that
    # pins this down is that a squared roof is an arc at constant radius, so
    # its two radii have to land on each other at any wall slope, and they only
    # do in this frame.
    def placed(fractions):
        points = roof.evaluate(fractions)[0]
        return (np.column_stack((base + points[:, 0], points[:, 1])),
                roof.derivative(fractions)[0])

    dense, _ = placed(np.linspace(0.0, 1.0, _ROOF_SAMPLES))
    envelope = float(np.hypot(dense[:, 0], dense[:, 1]).max())
    azimuth = np.arctan2(dense[:, 1], dense[:, 0])
    span = float(azimuth.max() - azimuth.min())

    # Green's theorem about the axis, over the roof alone.
    quadrature, slopes = placed(_NODES)
    swept = 0.5*abs((quadrature[:, 0]*slopes[:, 1]
                     - quadrature[:, 1]*slopes[:, 0]) @ _WEIGHTS)
    if span <= 0.0:
        raise ValueError("the exterior roof sweeps no angle about the axis")
    return envelope, float(swept), span


class OuterWall:
    """Outer radius of a channelled wall against axial station.

    Parameters
    ----------
    xs : array_like
        Axial stations [m], strictly increasing.
    envelope, equivalent : array_like
        The two radii at each station [m].
    name : str, optional
        Label for error messages.

    Notes
    -----
    Build one with :meth:`from_circuit` or :meth:`from_chamber` rather than
    directly, unless the wall is being stated by hand for a part pyskyfire did
    not resolve.
    """

    def __init__(self, xs, envelope, equivalent, *, name=None):
        xs = np.asarray(xs, dtype=float)
        radii = {kind: np.asarray(values, dtype=float)
                 for kind, values in (("envelope", envelope),
                                      ("equivalent", equivalent))}
        if xs.ndim != 1 or len(xs) < 2:
            raise ValueError("xs must be one-dimensional with at least two stations")
        if np.any(np.diff(xs) <= 0.0):
            raise ValueError("xs must be strictly increasing")
        for kind, values in radii.items():
            if values.shape != xs.shape:
                raise ValueError(f"{kind} must match xs in shape")
            if not np.isfinite(values).all() or np.any(values <= 0.0):
                raise ValueError(f"{kind} radii must be finite and positive")
        if np.any(radii["envelope"] < radii["equivalent"] - 1.0e-12):
            raise ValueError(
                "the envelope radius cannot fall inside the equivalent radius")

        self.name = name
        self.xs = xs
        self._radii = radii
        self._splines = {kind: PchipInterpolator(xs, values)
                         for kind, values in radii.items()}

    # -- construction -----------------------------------------------------

    @staticmethod
    def _measure(circuit, x):
        """``(envelope, equivalent, pitch)`` for one circuit at one ``x``.

        The equivalent radius is normalised by the angle the circuit's own
        channel sweeps, so it holds whether or not the circuit tiles the circle
        by itself. ``pitch`` is the total angle the circuit occupies, which is
        what weights it against any circuit interleaved with it.

        An outer liner is added to both radii rather than folded into the swept
        area: it is a layer laid over the finished surface, so it raises the
        crest and the mean by its own thickness.
        """
        prof = circuit.section_profiles_at(float(x))
        envelope, swept, span = _sector(circuit.cross_section, prof)
        liner = float(circuit.outer_wall_thickness(float(x)))
        return (envelope + liner,
                np.sqrt(2.0*swept/span) + liner,
                circuit.n_channels*float(prof.theta[0]))

    @classmethod
    def from_circuit(cls, circuit, *, n_stations=None):
        """Measure the outer wall of one resolved cooling circuit.

        The equivalent radius is normalised by the angle the circuit's own
        channel sweeps, so it answers "what cylinder encloses the same area per
        radian as this roof" without assuming the circuit tiles the circle by
        itself.
        """
        if circuit.x_domain is None:
            raise ValueError(
                "the circuit has no resolved geometry yet; build its "
                "ThrustChamber first")
        xs = np.asarray(circuit.x_domain, dtype=float)
        if n_stations is not None:
            xs = np.linspace(xs[0], xs[-1], int(n_stations))
        measured = np.array([cls._measure(circuit, x)[:2] for x in xs])
        return cls(xs, measured[:, 0], measured[:, 1], name=circuit.name)

    @classmethod
    def from_chamber(cls, chamber, *, n_stations=None):
        """Measure the outer wall across every circuit of a thrust chamber.

        Where circuits overlap in ``x`` they are interleaved around the
        circumference rather than stacked on top of each other, so the two
        radii combine differently. The envelope is whichever circuit reaches
        furthest out, since that is the crest a surrounding part has to clear.
        The equivalent radius is built from the area they jointly enclose, each
        circuit contributing over the angle its own channels occupy -- taking
        the larger of the two would claim the whole circumference looks like
        the deeper circuit, and taking one circuit's channel count would assume
        it had the circle to itself.
        """
        circuits = list(chamber.cooling_circuits)
        if not circuits:
            raise ValueError("this chamber has no cooling circuits to measure")
        spans = [(float(c.x_domain[0]), float(c.x_domain[-1])) for c in circuits]
        xs = np.unique(np.concatenate(
            [np.asarray(c.x_domain, dtype=float) for c in circuits]
            if n_stations is None else
            [np.linspace(lo, hi, int(n_stations)) for lo, hi in spans]))

        envelope = np.full(len(xs), -np.inf)
        moment = np.zeros(len(xs))
        angle = np.zeros(len(xs))
        for circuit, (lo, hi) in zip(circuits, spans, strict=True):
            for index, x in enumerate(xs):
                if not lo <= x <= hi:
                    continue
                crest, equivalent, pitch = cls._measure(circuit, x)
                envelope[index] = max(envelope[index], crest)
                # Areas add, and a radius contributes area as its square.
                moment[index] += pitch*equivalent**2
                angle[index] += pitch
        if not np.isfinite(envelope).all():
            raise ValueError(
                "the cooling circuits leave a gap in x, so the outer wall is "
                "undefined there; measure the circuits separately")
        return cls(xs, envelope, np.sqrt(moment/angle),
                   name=getattr(chamber, "name", None))

    # -- evaluation -------------------------------------------------------

    def _spline(self, kind):
        if kind not in OUTER_WALL_KINDS:
            raise ValueError(f"kind must be one of {OUTER_WALL_KINDS}; got {kind!r}")
        return self._splines[kind]

    @property
    def span(self):
        """First and last axial station the wall is defined over [m]."""
        return float(self.xs[0]), float(self.xs[-1])

    def r(self, x, kind="envelope"):
        """Outer radius at ``x`` [m]."""
        return self._spline(kind)(x)

    def dr_dx(self, x, kind="envelope"):
        """Slope of the outer wall at ``x``."""
        return self._spline(kind).derivative()(x)

    def normal(self, x, kind="envelope"):
        """Outward unit normal in the meridional plane, as ``(n_x, n_r)``.

        Radially outward where the wall runs parallel to the axis, and tipped
        toward the narrowing end where it does not.
        """
        slope = np.asarray(self.dr_dx(x, kind), dtype=float)
        scale = np.sqrt(1.0 + slope**2)
        return np.stack((-slope/scale, np.ones_like(slope)/scale), axis=-1)

    def profile(self, x, kind="envelope", *, reach=None, samples=401):
        """The wall around station ``x``, ready for a manifold section.

        Returns a :class:`~pyskyfire.ducts.half_moon.InnerProfile`: the wall
        measured by arc length from ``x`` rather than by station, because that
        is the room a cooling-channel breakout actually has on the surface. On
        a sloped wall the surface runs further than the station span does.

        Parameters
        ----------
        x : float
            Anchor station [m], strictly inside the wall's span.
        kind : {'envelope', 'equivalent'}, optional
            Which radius the wall follows. The default clears the crests,
            which is what a cavity wrapped around the chamber has to do.
        reach : float, optional
            Arc length to carry either side of the anchor [m]. Defaults to
            whatever the wall offers, which is all a section can ever use.
        samples : int, optional
            Stations sampled across that reach.
        """
        from pyskyfire.ducts import InnerProfile

        low, high = self.span
        if not low < x < high:
            raise ValueError(
                f"x must lie strictly inside the wall's span {low:g}..{high:g} m, "
                "so that a section anchored there has wall on both sides")
        stations = np.linspace(low, high, int(samples))
        points = np.column_stack((stations, self.r(stations, kind)))
        steps = np.linalg.norm(np.diff(points, axis=0), axis=1)
        arc = np.concatenate(([0.0], np.cumsum(steps)))
        arc = arc - np.interp(x, stations, arc)
        if reach is not None:
            reach = float(reach)
            if not np.isfinite(reach) or reach <= 0.0:
                raise ValueError("reach must be finite and positive")
            keep = np.abs(arc) <= reach
            # Keep the samples straddling the anchor even at a tiny reach.
            keep[max(np.argmax(arc >= 0.0) - 1, 0)] = True
            keep[min(np.argmax(arc >= 0.0) + 1, len(arc) - 1)] = True
            arc, points = arc[keep], points[keep]
        return InnerProfile(arc, points)

    def normal_angle(self, x, kind="envelope"):
        """Signed angle of the outward normal from radial [rad].

        Positive toward ``+x``. Unlike
        :meth:`~pyskyfire.regen.contour.Contour.normal_angle`, which reports an
        unsigned angle, this keeps the sign: a manifold has to know which way
        along the axis the wall tips, not merely by how much.
        """
        return np.arctan(-np.asarray(self.dr_dx(x, kind), dtype=float))
