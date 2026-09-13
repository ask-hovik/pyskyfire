"""Geometry-only transverse sections and continuous hot-wall profiles.

No solver dimensions, paths or section analytics are modified here. Oblique
sections represent the local prismatic approximation, not exact CAD sweep cuts.
"""

from dataclasses import dataclass, replace
import warnings

import numpy as np

from .channel_construction import build_channel_construction
from .channel_helpers import build_axis_normal_section_planes


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

    def segment(self, first, last):
        """The part of this arc between two fractions of its parameter span."""
        return self.trimmed(first, 1 - last)

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


def _de_casteljau(control_points, fraction):
    """Split a cubic Bezier, returning the control points of both halves."""
    p = np.asarray(control_points, dtype=float)
    t = float(fraction)
    a = p[:-1] + t * np.diff(p, axis=0)
    b = a[:-1] + t * np.diff(a, axis=0)
    c = b[0] + t * (b[1] - b[0])
    return (
        np.array([p[0], a[0], b[0], c]),
        np.array([c, b[1], a[2], p[3]]),
    )


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

    def segment(self, first, last):
        """The part of this curve between two parameters, by de Casteljau.

        Subdivision is exact, so a split profile still describes the same
        curve as the profile it was split from.
        """
        p = self.control_points
        if last <= first:
            raise ValueError("A Bezier segment must span an increasing parameter")
        # Split at ``last`` first, then reparameterize ``first`` onto the piece
        # that is left, so both cuts come from one pair of de Casteljau passes.
        head = _de_casteljau(p, last)[0]
        return type(self)(_de_casteljau(head, first / last)[1])

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


def hot_wall_stations(
    spans, samples, contour_span, *, transition_fraction=0.02, transition_length_m=None
):
    """Physical hand-offs, including gaps, independent of sample spacing.

    Each uncovered internal interval is one hand-off. Other internal circuit
    endpoints are individual hand-offs. Intervals are centered and never clipped.
    """
    lo, hi = map(float, contour_span)
    length = hi - lo
    if not np.isfinite(length) or length <= 0:
        raise ValueError("Full contour must have positive finite axial length")
    if transition_length_m is None:
        if not np.isfinite(transition_fraction) or not 0 < transition_fraction <= 1:
            raise ValueError("Hot-wall transition fraction must be in (0, 1]")
        transition = transition_fraction * length
    else:
        transition = float(transition_length_m)
        if not np.isfinite(transition) or transition <= 0:
            raise ValueError("Hot-wall transition length must be positive and finite")
    spans = sorted((float(a), float(b)) for a, b in spans)
    if not spans or any(
        not np.isfinite(a + b) or b <= a or a < lo or b > hi for a, b in spans
    ):
        raise ValueError(
            "Circuit spans must be nonempty and lie within the full contour"
        )
    union = []
    for a, b in spans:
        if union and a <= union[-1][1]:
            union[-1][1] = max(union[-1][1], b)
        else:
            union.append([a, b])
    uncovered = max(0.0, 1 - sum(b - a for a, b in union) / length)
    if uncovered > 0.1 + 1e-12:
        warnings.warn(
            f"{uncovered:.1%} of the full contour length is not covered by cooling circuits",
            UserWarning,
            stacklevel=2,
        )
    start, end = union[0][0], union[-1][1]
    gaps = [(a[1], b[0]) for a, b in zip(union[:-1], union[1:])]
    gap_ends = {v for gap in gaps for v in gap}
    events = [
        (v, v)
        for v in sorted({v for span in spans for v in span})
        if start < v < end and v not in gap_ends
    ]
    events += gaps
    intervals = []
    for a, b in sorted(events):
        mid = (a + b) / 2
        width = max(transition, b - a)
        left, right = mid - width / 2, mid + width / 2
        # Preserve exact gap endpoints when the gap sets the width.
        if b - a >= transition:
            left, right = a, b
        if left < start or right > end:
            raise ValueError(
                "Hot-wall transition extends beyond circuit-covered extent; reduce transition length"
            )
        if intervals and left < intervals[-1][1] - 1e-12 * length:
            raise ValueError(
                "Overlapping hot-wall transition intervals; reduce transition length"
            )
        intervals.append((left, right))
    stations = np.unique(
        np.concatenate(
            [np.asarray(samples), [start, end], np.asarray(intervals).ravel()]
        )
    )
    stations = stations[(stations >= start) & (stations <= end)]
    for a, b in intervals:
        stations = stations[~((stations > a) & (stations < b))]
    return stations, intervals, uncovered


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


def _azimuths(points):
    return np.arctan2(np.asarray(points)[..., 1], np.asarray(points)[..., 0])


def _fraction_at_azimuth(curve, target, samples=129):
    """Fraction at which a curve reaches ``target`` azimuth about the axis.

    ``join_hot_wall`` checks that a profile never folds in azimuth, so each of
    its curves turns monotonically and this root is unique.
    """
    grid = np.linspace(0.0, 1.0, samples)
    angles = np.unwrap(_azimuths(curve.point(grid)))
    if angles[-1] <= angles[0]:
        raise ValueError("Hot-wall curve does not advance in azimuth")
    target = angles[0] + (target - angles[0]) % (2 * np.pi)
    if target > angles[-1]:
        raise ValueError("Azimuth lies outside the curve it was located on")
    index = max(1, int(np.searchsorted(angles, target)))
    low, high = grid[index - 1], grid[index]
    base = angles[index - 1]
    for _ in range(60):
        middle = (low + high) / 2
        raw = float(_azimuths(curve.point(middle)))
        value = base + (raw - base + np.pi) % (2 * np.pi) - np.pi
        if value < target:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def _ordered_profile(curves):
    """Rotate a closed profile so its start azimuths increase from the seam."""
    starts = np.array([float(_azimuths(c.point(0))) % (2 * np.pi) for _, c in curves])
    seam = int(np.argmin(starts))
    return list(curves[seam:]) + list(curves[:seam]), np.roll(starts, -seam)


def _azimuth_extents(curves):
    return np.array(
        [
            np.ptp(np.unwrap(_azimuths(c.point(np.linspace(0, 1, 33)))))
            for _, c in curves
        ]
    )


def _cluster_azimuths(values, tolerance):
    """Group nearly coincident azimuths, closing the group across the seam."""
    order = np.argsort(values)
    groups = []
    for index in order:
        if groups and values[index] - values[groups[-1][-1]] <= tolerance:
            groups[-1].append(index)
        else:
            groups.append([index])
    if len(groups) > 1 and (
        values[groups[0][0]] + 2 * np.pi - values[groups[-1][-1]] <= tolerance
    ):
        groups[0] = groups.pop() + groups[0]
    return groups


def hot_wall_bridge(first, second):
    """Split two profiles of unlike topology onto one set of azimuth cuts.

    A circuit hand-off replaces one channel pattern with another, so the two
    profiles bounding it share no curves and nothing pairs them by index. Both
    do wind once about the axis, though, so azimuth pairs them: cut each at
    every curve start of either profile and the two carry the same number of
    curves, spanning the same angular intervals, in the same order. Splitting
    is exact, and a profile keeps its own curve starts wherever the two nearly
    coincide, so neither side gains a sliver the other's cut would have made.
    """
    profiles = [_ordered_profile(first), _ordered_profile(second)]
    extents = np.concatenate([_azimuth_extents(curves) for curves, _ in profiles])
    if np.any(extents <= 0):
        raise ValueError("Hot-wall bridge profiles must advance in azimuth")
    tolerance = 0.25 * float(extents.min())

    marks = np.concatenate([starts for _, starts in profiles])
    offsets = np.cumsum([0] + [len(starts) for _, starts in profiles])
    sides = np.concatenate(
        [np.full(len(starts), s) for s, (_, starts) in enumerate(profiles)]
    )
    # Groups are already in ascending azimuth, and no group's azimuths reach
    # into another's, so this one order runs both profiles the same way round.
    # Ordering by either profile's own cuts instead would swap a pair wherever
    # a borrowed cut fell just the other side of a curve start.
    groups = _cluster_azimuths(marks, tolerance)
    located = [[], []]
    for group in groups:
        centre = float(np.mean(np.unwrap(marks[group])) % (2 * np.pi))
        for side in (0, 1):
            own = [index for index in group if sides[index] == side]
            if len(own) > 1:
                raise ValueError(
                    "Hot-wall bridge profile has two curve starts at one azimuth"
                )
            # A profile's own curve starts are placed, not searched for. Going
            # back through an azimuth would round them onto the curve before.
            located[side].append(
                (own[0] - offsets[side], 0.0)
                if own
                else _locate_azimuth(*profiles[side], centre)
            )
    names = [f"bridge_span_{k}" for k in range(len(groups))]
    split = []
    for side in (0, 1):
        curves = profiles[side][0]
        pieces = []
        for k, (index, fraction) in enumerate(located[side]):
            following, next_fraction = located[side][(k + 1) % len(groups)]
            end = next_fraction if following == index else 1.0
            pieces.append((names[k], curves[index][1].segment(fraction, end)))
        _check_closed(pieces)
        split.append(pieces)
    return names, split[0], split[1]


def _check_closed(pieces):
    for (_, first), (_, second) in zip(pieces, pieces[1:] + pieces[:1]):
        if not np.allclose(first.point(1), second.point(0), atol=1e-10, rtol=0):
            raise ValueError("Split hot-wall bridge profile is not closed")


def _locate_azimuth(curves, starts, azimuth):
    """Which curve of an ordered profile an azimuth falls on, and where."""
    value = starts[0] + (azimuth - starts[0]) % (2 * np.pi)
    index = int(np.searchsorted(starts, value, side="right")) - 1
    if np.isclose(value, starts[index], rtol=0, atol=1e-12):
        return index, 0.0
    return index, _fraction_at_azimuth(curves[index][1], value)


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


def topology_runs(sections):
    """Group consecutive sections sharing one ordered list of curve names.

    ``join_hot_wall`` pins the seam to a stable channel, so inside a run curve
    ``j`` belongs to the same channel at every station. CAD can then pair
    curves by index instead of matching vertices by proximity, which aliases
    onto a neighbouring channel once a helix advances past half a pitch. A
    circuit hand-off changes the name list and starts a new run.
    """
    runs = []
    for index, section in enumerate(sections):
        names = [curve["name"] for curve in section["curves"]]
        if runs and runs[-1]["curveNames"] == names:
            runs[-1]["sectionIndices"].append(index)
        else:
            runs.append({"sectionIndices": [index], "curveNames": names})
    return runs


def build_hot_wall(
    chamber,
    resolved,
    *,
    trim_fraction=0.1,
    transition_fraction=0.02,
    transition_length_m=None,
):
    channels = resolved.channels
    spans = [(c.reference_path[0, 0], c.reference_path[-1, 0]) for c in channels]
    xs = np.asarray(chamber.contour.xs)
    stations, intervals, uncovered = hot_wall_stations(
        spans,
        np.concatenate([c.reference_path[:, 0] for c in channels]),
        (xs[0], xs[-1]),
        transition_fraction=transition_fraction,
        transition_length_m=transition_length_m,
    )
    cached = []
    for channel, (a, b) in zip(channels, spans):
        path = channel.reference_path
        grid = np.unique(np.r_[path[:, 0], stations[(stations >= a) & (stations <= b)]])
        sample = np.column_stack(
            (
                grid,
                np.interp(grid, path[:, 0], path[:, 1]),
                np.interp(grid, path[:, 0], path[:, 2]),
            )
        )
        circuit = chamber.cooling_circuits[channel.circuit_index]
        if hasattr(chamber.contour, "r"):
            sample[:, 1] = chamber.contour.r(grid)
        construction = build_channel_construction(circuit, sample)
        cached.append((construction, transverse_construction(construction)))
    sections = []
    profiles = []
    for x in stations:
        pieces = []
        for channel, (a, b), (construction, transverse) in zip(channels, spans, cached):
            if not a <= x <= b:
                continue
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
        if not pieces:
            continue  # Internal uncovered spans are bridged by the single loft.
        try:
            curves = join_hot_wall(pieces, trim_fraction)
        except ValueError as exc:
            raise ValueError(f"Hot wall at x={x:.9g} m: {exc}") from exc
        sections.append(
            {
                "name": f"hot_wall_{len(sections)}",
                "plane": {
                    "origin_mm": [0.0, 0.0, float(x * 1000)],
                    "normal": [0.0, 0.0, 1.0],
                    "xAxis": [1.0, 0.0, 0.0],
                },
                "curves": [dict(curve.to_data(), name=name) for name, curve in curves],
            }
        )
        profiles.append(curves)
    runs = topology_runs(sections)
    bridges = []
    for before, after in zip(runs, runs[1:]):
        i, j = before["sectionIndices"][-1], after["sectionIndices"][0]
        try:
            names, first, second = hot_wall_bridge(profiles[i], profiles[j])
        except ValueError as exc:
            raise ValueError(
                f"Hot-wall hand-off between x={stations[i]:.9g} m and "
                f"x={stations[j]:.9g} m: {exc}"
            ) from exc
        bridges.append(
            {
                "beforeSection": i,
                "afterSection": j,
                "curveNames": names,
                "sections": [
                    {
                        "plane": sections[index]["plane"],
                        "curves": [
                            dict(curve.to_data(), name=name) for name, curve in side
                        ],
                    }
                    for index, side in ((i, first), (j, second))
                ],
            }
        )
    return {
        "trimFraction": trim_fraction,
        "transitionIntervals_mm": (
            np.asarray(intervals).reshape(-1, 2) * 1000
        ).tolist(),
        "uncoveredLengthFraction": uncovered,
        "sections": sections,
        "runs": runs,
        "bridges": bridges,
    }
