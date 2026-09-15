"""Geometry-only transverse sections and hot-wall surfaces.

No solver dimensions, paths or section analytics are modified here. Oblique
sections represent the local prismatic approximation, not exact CAD sweep cuts.
"""

from dataclasses import dataclass, replace
import warnings

import numpy as np

from .channel_construction import build_channel_construction
from .channel_helpers import build_axis_normal_section_planes
from pyskyfire.common.section_geometry import SectionLine


@dataclass(frozen=True)
class PlanarLine:
    """Directed line, sharing the sampling/export contract of PlanarArc."""

    start: np.ndarray
    end: np.ndarray

    def point(self, fraction):
        return self.start + np.asarray(fraction)[..., None]*(self.end-self.start)

    def tangent(self, fraction):
        vector = self.end-self.start
        return vector/np.linalg.norm(vector)

    def trimmed(self, left, right):
        return PlanarLine(self.point(left), self.point(right))

    def transformed(self, matrix, origin):
        return PlanarLine(origin+matrix@self.start, origin+matrix@self.end)

    def to_data(self):
        return {"kind": "line", "start_uv_mm": (self.start*1000).tolist(),
                "end_uv_mm": (self.end*1000).tolist()}


@dataclass(frozen=True)
class PlanarArc:
    """Affine circular arc: c + a*cos(t) + b*sin(t), in metres."""

    center: np.ndarray
    a: np.ndarray
    b: np.ndarray
    start: float
    end: float

    def point(self, fraction):
        t = self.start + np.asarray(fraction) * (self.end - self.start)
        return (
            self.center + np.cos(t)[..., None] * self.a + np.sin(t)[..., None] * self.b
        )

    def tangent(self, fraction):
        t = self.start + fraction * (self.end - self.start)
        v = (-np.sin(t) * self.a + np.cos(t) * self.b) * (self.end - self.start)
        return v / np.linalg.norm(v)

    def trimmed(self, left, right):
        delta = self.end - self.start
        return replace(
            self, start=self.start + left * delta, end=self.end - right * delta
        )

    def transformed(self, matrix, origin):
        return replace(
            self,
            center=origin + matrix @ self.center,
            a=matrix @ self.a,
            b=matrix @ self.b,
        )

    def to_data(self):
        """Principal-axis ellipse parameters, retaining directed endpoints."""
        basis = np.column_stack((self.a, self.b))
        u, radii, _ = np.linalg.svd(basis)
        major = u[:, 0]
        minor = np.array([-major[1], major[0]])
        axes = np.column_stack((major, minor))
        q = (axes.T @ basis) / radii[:, None]
        start_vector = q @ np.array([np.cos(self.start), np.sin(self.start)])
        start = np.arctan2(start_vector[1], start_vector[0])
        end = start + np.sign(np.linalg.det(q)) * (self.end - self.start)
        circle = np.isclose(radii[0], radii[1], rtol=1e-12, atol=0)
        return {
            "kind": "circle" if circle else "ellipse",
            "center_uv_mm": (self.center * 1000).tolist(),
            "axisA_uv_mm": (self.a * 1000).tolist(),
            "axisB_uv_mm": (self.b * 1000).tolist(),
            "startAngle_deg": float(np.degrees(self.start)),
            "endAngle_deg": float(np.degrees(self.end)),
            "radius_mm": float(radii[0] * 1000),
            "majorRadius_mm": float(radii[0] * 1000),
            "minorRadius_mm": float(radii[1] * 1000),
            "majorAxis": major.tolist(),
            "startParameter": float(min(start, end) / (2 * np.pi) % 1),
            "endParameter": float(max(start, end) / (2 * np.pi) % 1),
        }


@dataclass(frozen=True)
class TangentBezier:
    control_points: np.ndarray

    @classmethod
    def between(cls, first, second):
        p, q = first.point(1), second.point(0)
        handle = np.linalg.norm(q - p) / 3
        if handle <= 1e-14:
            raise ValueError(
                "Blend has coincident endpoints; increase hot-wall trim fraction"
            )
        return cls(
            np.array(
                [p, p + handle * first.tangent(1), q - handle * second.tangent(0), q]
            )
        )

    def point(self, fraction):
        t = np.asarray(fraction)[..., None]
        p = self.control_points
        return (
            (1 - t) ** 3 * p[0]
            + 3 * (1 - t) ** 2 * t * p[1]
            + 3 * (1 - t) * t**2 * p[2]
            + t**3 * p[3]
        )

    def tangent(self, fraction):
        p = self.control_points
        t = float(fraction)
        v = (
            3 * (1 - t) ** 2 * (p[1] - p[0])
            + 6 * (1 - t) * t * (p[2] - p[1])
            + 3 * t**2 * (p[3] - p[2])
        )
        return v / np.linalg.norm(v)

    def to_data(self):
        return {
            "kind": "bezier",
            "controlPoints_uv_mm": (self.control_points * 1000).tolist(),
        }


@dataclass(frozen=True)
class TransverseConstruction:
    source: object
    planes: object
    inner_arcs: tuple
    outer_arcs: tuple

    def polygon(self, index, n_arc_points=8):
        uv = np.vstack(
            [
                arc.point(np.linspace(0, 1, n_arc_points))
                for arc in (self.inner_arcs[index], self.outer_arcs[index])
            ]
        )
        return (
            self.planes.origins[index]
            + uv[:, :1] * self.planes.x_axes[index]
            + uv[:, 1:] * self.planes.y_axes[index]
        )


def transverse_matrix(construction, planes, index):
    if construction.section_type == "squared":
        return np.eye(2)
    flow = construction.flow_planes
    normal = flow.normals[index]
    # Intersect the locally extruded section with X=constant. This includes
    # the rotated lobe bisector as well as the slope/helix shear term.
    basis = np.column_stack((flow.x_axes[index], flow.y_axes[index]))
    projected = basis - normal[:, None] * (basis[0] / normal[0])[None, :]
    target = np.stack((planes.x_axes[index], planes.y_axes[index]))
    return target @ projected


def transverse_construction(construction):
    planes = build_axis_normal_section_planes(
        construction.profiles.cross_section_base_point_path
    )
    families = []
    for source in (construction.inner_arc, construction.outer_arc):
        arcs = []
        for i in range(len(planes)):
            matrix = transverse_matrix(construction, planes, i)
            if isinstance(source, SectionLine):
                arcs.append(PlanarLine(matrix@source.starts[i], matrix@source.ends[i]))
                continue
            radius = source.radii[i]
            arcs.append(
                PlanarArc(
                    matrix @ source.centers[i],
                    matrix[:, 0] * radius,
                    matrix[:, 1] * radius,
                    source.start_angles[i],
                    source.end_angles[i],
                )
            )
        families.append(tuple(arcs))
    return TransverseConstruction(construction, planes, *families)


def _span_union(spans):
    union = []
    for a, b in sorted(spans):
        if union and a <= union[-1][1]:
            union[-1][1] = max(union[-1][1], b)
        else:
            union.append([a, b])
    return union


def hot_wall_regions(spans, squared, contour_span, *, transition_fraction=0.02):
    """Partition covered spans by construction type and active channel seeds.

    With rounded channels, changes into/out of shared spans leave a centered
    open transition zone. Its width is a fraction of the full contour length,
    matching the former blend-zone convention. Empty spans are never connected.
    Squared-only spans can share a revolve regardless of channel count or winding.
    """
    lo, hi = map(float, contour_span)
    length = hi - lo
    if not np.isfinite(length) or length <= 0:
        raise ValueError("Full contour must have positive finite axial length")
    fraction = float(transition_fraction)
    if not np.isfinite(fraction) or not 0 < fraction <= 1:
        raise ValueError("Hot-wall transition fraction must be in (0, 1]")
    half_width = fraction * length / 2
    if (
        not spans
        or len(spans) != len(squared)
        or any(not np.isfinite(a + b) or b <= a or a < lo or b > hi for a, b in spans)
    ):
        raise ValueError(
            "Circuit spans must be nonempty and lie within the full contour"
        )
    union = _span_union(spans)
    uncovered = max(0.0, 1 - sum(b - a for a, b in union) / length)
    if uncovered > 0.1 + 1e-12:
        warnings.warn(
            f"{uncovered:.1%} of the full contour length is not covered by cooling circuits",
            UserWarning,
            stacklevel=2,
        )
    ends = sorted({v for span in spans for v in span})
    cells = []
    for a, b in zip(ends, ends[1:]):
        active = tuple(
            i for i, (start, stop) in enumerate(spans) if start <= a and b <= stop
        )
        cells.append((a, b, active))
    sites = []
    if not all(squared):
        for before, after in zip(cells, cells[1:]):
            if before[2] and after[2] and max(len(before[2]), len(after[2])) > 1:
                sites.append(before[1])
    exclusions = _span_union(
        [(max(lo, x - half_width), min(hi, x + half_width)) for x in sites]
    )
    regions = []
    for a, b, active in cells:
        if not active:
            continue
        pieces = [(a, b)]
        for left, right in exclusions:
            remaining = []
            for start, stop in pieces:
                if right <= start or left >= stop:
                    remaining.append((start, stop))
                else:
                    if start < left:
                        remaining.append((start, left))
                    if right < stop:
                        remaining.append((right, stop))
            pieces = remaining
        mode = "revolve" if all(squared[i] for i in active) else "loft"
        for start, stop in pieces:
            if (
                regions
                and mode == regions[-1]["mode"] == "revolve"
                and regions[-1]["end"] == start
            ):
                regions[-1]["end"] = stop
            else:
                regions.append(dict(start=start, end=stop, active=active, mode=mode))
    return regions, exclusions, uncovered


def _interval_bezier(contour, a, b):
    """The cubic Bezier through four evenly spaced samples of one interval."""
    x = np.linspace(a, b, 4)
    r = np.asarray(contour.r(x))
    inner = np.linalg.solve(
        [[12, 6], [6, 12]],
        [27 * r[1] - 8 * r[0] - r[3], 27 * r[2] - r[0] - 8 * r[3]],
    )
    return TangentBezier(np.column_stack(([r[0], *inner, r[3]], x)))


def _secant_tangent(contour, x, step):
    """Unit (radius, axial) tangent of the meridian, from a short secant.

    ``step`` is signed: a negative step reads the interval behind ``x``. The
    result always points in the direction of increasing axial position, so a
    contour corner shows up as two differing tangents at one knot.
    """
    r_here, r_there = np.asarray(contour.r(np.array([x, x + step])), dtype=float)
    tangent = np.array([r_there - r_here, step]) * np.sign(step)
    return tangent / np.linalg.norm(tangent)


def contour_corner_knots(contour, grid, *, tolerance=1e-3):
    """Interior grid positions where the meridian's tangent jumps.

    The throat carries independent slopes on either side, and interpolated
    contours corner at every input point. A spline may only be fitted through
    a smooth run, so these are the positions the run has to be cut at.
    """
    corners = []
    for previous, x, following in zip(grid, grid[1:], grid[2:]):
        step = 1e-4 * min(x - previous, following - x)
        behind = _secant_tangent(contour, x, -step)
        ahead = _secant_tangent(contour, x, step)
        if np.dot(behind, ahead) < np.cos(tolerance):
            corners.append(x)
    return corners


def hot_wall_revolve(contour, start, end):
    """Meridional splines in (radius, axial) mm, one per smooth contour run.

    Each spline interpolates the contour's own knots, so CAD receives one
    curve where the contour is smooth instead of a chain of short segments
    that reads as faceted. Runs are cut at tangent discontinuities, which
    keeps the throat corner sharp, and a trimmed span contributes its exact
    endpoints without depending on coolant sampling.
    """
    knots = np.asarray(contour.xs)
    grid = np.unique(np.r_[start, knots[(knots > start) & (knots < end)], end])
    cuts = [start, *contour_corner_knots(contour, grid), end]
    curves = []
    for a, b in zip(cuts, cuts[1:]):
        x = grid[(grid >= a) & (grid <= b)]
        if len(x) < 2:
            continue
        # A run of a single interval has too few points to fit a spline
        # through, so it keeps the cubic Bezier that reproduces it exactly.
        if len(x) == 2:
            curves.append(_interval_bezier(contour, x[0], x[1]).to_data())
            continue
        points = np.column_stack((np.asarray(contour.r(x), dtype=float), x)) * 1000
        curves.append({"kind": "spline", "points_uv_mm": points.tolist()})
    return {"span_mm": [float(start * 1000), float(end * 1000)], "curves": curves}


def _rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s], [s, c]])


def _hot_arc(construction, transverse, index, rotation):
    prof = construction.profiles
    if construction.section_type == "squared":
        radius = construction.reference_path[index, 1]
        theta = construction.reference_path[index, 2] + rotation
        pitch = construction.pitch[index]
        return PlanarArc(
            np.zeros(2),
            np.array([radius, 0.0]),
            np.array([0.0, radius]),
            theta - pitch / 2,
            theta + pitch / 2,
        )
    d, half = prof.sector_radius[index], prof.theta[index] / 2
    liner = np.broadcast_to(prof.t_hot_layer, (len(transverse.planes),))[index]
    radius = d * np.sin(half) + liner
    root = np.sqrt(max(0.0, radius**2 - (d * np.sin(half)) ** 2))
    intersection = d * np.cos(half) - root
    upper = np.arctan2(intersection * np.sin(half), intersection * np.cos(half) - d)
    arc = PlanarArc(
        np.zeros(2),
        np.array([radius, 0.0]),
        np.array([0.0, radius]),
        2 * np.pi - upper,
        upper,
    )
    planes = transverse.planes
    world = np.column_stack((planes.x_axes[index, 1:], planes.y_axes[index, 1:]))
    matrix = world @ transverse_matrix(construction, planes, index)
    return arc.transformed(
        _rotation(rotation) @ matrix, _rotation(rotation) @ planes.origins[index, 1:]
    )


def join_hot_wall(pieces, trim_fraction):
    """Join angularly ordered (id, type, arc, reference angle) pieces."""
    if not np.isfinite(trim_fraction) or not 0 < trim_fraction < 0.5:
        raise ValueError("Hot-wall trim fraction must be in (0, 0.5)")
    pieces = sorted(pieces, key=lambda p: p[3] % (2 * np.pi))
    if not pieces:
        raise ValueError("A hot-wall profile requires active cooling channels")
    # Keep the profile seam on a stable channel as a helix crosses theta=0.
    # This supplies deterministic sketch ordering without adding loft guides.
    seam = min(range(len(pieces)), key=lambda i: pieces[i][0])
    pieces = pieces[seam:] + pieces[:seam]
    arcs = [p[2] for p in pieces]
    compatible = []
    for i, piece in enumerate(pieces):
        j = (i + 1) % len(pieces)
        same = piece[1] == pieces[j][1] == "squared" and np.isclose(
            np.linalg.norm(arcs[i].a), np.linalg.norm(arcs[j].a), rtol=1e-10, atol=1e-12
        )
        compatible.append(same)
        if same:
            # Complete the wall over ribs and unallocated angular space.
            end = arcs[i].end
            start = arcs[j].start + 2 * np.pi * round(
                (end - arcs[j].start) / (2 * np.pi)
            )
            seam = (end + start) / 2
            arcs[i] = replace(arcs[i], end=seam)
            arcs[j] = replace(arcs[j], start=arcs[j].start + seam - start)
    arcs = [
        arc.trimmed(
            0 if compatible[i - 1] else trim_fraction,
            0 if compatible[i] else trim_fraction,
        )
        for i, arc in enumerate(arcs)
    ]
    result = []
    for i, arc in enumerate(arcs):
        result.append((pieces[i][0], arc))
        if not compatible[i]:
            result.append(
                (
                    pieces[i][0] + "_blend",
                    TangentBezier.between(arc, arcs[(i + 1) % len(arcs)]),
                )
            )
    for (_, first), (_, second) in zip(result, result[1:] + result[:1]):
        if not np.allclose(first.point(1), second.point(0), atol=1e-10, rtol=0):
            raise ValueError("Hot-wall profile is not closed")
    # A simple chamber contour must wind once without folding in azimuth.
    points = np.concatenate(
        [curve.point(np.linspace(0, 1, 17)[:-1]) for _, curve in result]
    )
    points = np.vstack([points, points[0]])
    angles = np.unwrap(np.arctan2(points[:, 1], points[:, 0]))
    if np.any(np.diff(angles) < -1e-9) or not np.isclose(
        angles[-1] - angles[0], 2 * np.pi
    ):
        raise ValueError(
            "Hot-wall profile folds or overlaps; adjust trim fraction or channel layout"
        )
    return result


def _station_index(stations, x):
    """Index of the construction sample at ``x``.

    Every hot-wall station is one of the samples its channel was constructed
    on, so this is an exact lookup. Silently taking a neighbouring station's
    geometry would mis-place the arc, so an inexact match is an error.
    """
    index = int(np.searchsorted(stations, x))
    if index and (
        index == len(stations)
        or abs(stations[index - 1] - x) <= abs(stations[index] - x)
    ):
        index -= 1
    if abs(stations[index] - x) > 1e-9 * max(1.0, abs(x)):
        raise ValueError(f"Hot-wall station x={x:.9g} m is not a construction sample")
    return index


def build_hot_wall(
    chamber,
    resolved,
    *,
    trim_fraction=0.1,
    transition_fraction=0.02,
):
    """Export contour revolves and independent, circumferentially blended lofts."""
    from .cross_section import CrossSectionSquared

    channels = resolved.channels
    # Multiple seeds of one circuit must not create fictitious interleave sites.
    circuit_indices = sorted({c.circuit_index for c in channels})
    spans = []
    squared = []
    for index in circuit_indices:
        seeds = [c for c in channels if c.circuit_index == index]
        spans.append(
            (
                min(c.reference_path[0, 0] for c in seeds),
                max(c.reference_path[-1, 0] for c in seeds),
            )
        )
        squared.append(
            isinstance(
                chamber.cooling_circuits[index].cross_section, CrossSectionSquared
            )
        )
    xs = np.asarray(chamber.contour.xs)
    regions, exclusions, uncovered = hot_wall_regions(
        spans, squared, (xs[0], xs[-1]), transition_fraction=transition_fraction
    )
    sections, runs, revolves = [], [], []
    for region in regions:
        a, b = region["start"], region["end"]
        if region["mode"] == "revolve":
            revolves.append(hot_wall_revolve(chamber.contour, a, b))
            continue
        active_circuits = {circuit_indices[i] for i in region["active"]}
        active = [c for c in channels if c.circuit_index in active_circuits]
        samples = np.concatenate([c.reference_path[:, 0] for c in active])
        stations = np.unique(np.r_[a, samples[(samples > a) & (samples < b)], b])
        cached = []
        for channel in active:
            path = channel.reference_path
            grid = np.unique(np.r_[path[:, 0], stations])
            sample = np.column_stack(
                (grid, chamber.contour.r(grid), np.interp(grid, path[:, 0], path[:, 2]))
            )
            circuit = chamber.cooling_circuits[channel.circuit_index]
            construction = build_channel_construction(circuit, sample)
            cached.append((construction, transverse_construction(construction)))
        indices = []
        for x in stations:
            pieces = []
            for channel, (construction, transverse) in zip(active, cached):
                i = _station_index(construction.reference_path[:, 0], x)
                for instance in range(channel.pattern_count):
                    rotation = instance * channel.pattern_angle_rad
                    arc = _hot_arc(construction, transverse, i, rotation)
                    pieces.append(
                        (
                            f"{channel.name}_{instance}",
                            construction.section_type,
                            arc,
                            construction.reference_path[i, 2] + rotation,
                        )
                    )
            try:
                curves = join_hot_wall(pieces, trim_fraction)
            except ValueError as exc:
                raise ValueError(f"Hot wall at x={x:.9g} m: {exc}") from exc
            indices.append(len(sections))
            sections.append(
                {
                    "name": f"hot_wall_{len(sections)}",
                    "plane": {
                        "origin_mm": [0.0, 0.0, float(x * 1000)],
                        "normal": [0.0, 0.0, 1.0],
                        "xAxis": [1.0, 0.0, 0.0],
                    },
                    "curves": [
                        dict(curve.to_data(), name=name) for name, curve in curves
                    ],
                }
            )
        names = [curve["name"] for curve in sections[indices[0]]["curves"]]
        if any(
            [curve["name"] for curve in sections[i]["curves"]] != names for i in indices
        ):
            raise ValueError("Hot-wall curve order changes within a loft span")
        runs.append({"sectionIndices": indices, "curveNames": names})
    wall = {
        "trimFraction": trim_fraction,
        "interleaveGapFraction": float(transition_fraction),
        "excludedInterleaveIntervals_mm": (
            np.asarray(exclusions).reshape(-1, 2) * 1000
        ).tolist(),
        "uncoveredLengthFraction": uncovered,
        "fullContourAllowed": bool(all(squared)),
        "revolves": revolves,
        "sections": sections,
        "runs": runs,
    }
    if wall["fullContourAllowed"]:
        wall["fullContour"] = hot_wall_revolve(chamber.contour, xs[0], xs[-1])
    return wall
