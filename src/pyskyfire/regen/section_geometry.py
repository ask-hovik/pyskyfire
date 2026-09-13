"""Planar arc geometry used by section analytics and physical construction."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class SectionArc:
    centers: np.ndarray
    radii: np.ndarray
    start_angles: np.ndarray
    end_angles: np.ndarray

    def points(self, index, count):
        angles = np.linspace(self.start_angles[index], self.end_angles[index], count)
        return self.centers[index] + self.radii[index] * np.column_stack(
            (np.cos(angles), np.sin(angles))
        )


_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(24)
_NODES = (_NODES + 1.0) / 2.0
_WEIGHTS = _WEIGHTS / 2.0


def arc_boundary_metrics(inner, outer, transform=None):
    """Area, centroid, and edge lengths of arc/line/arc/line boundaries.

    Green's theorem integrates the actual arcs, without a polygon-resolution
    dependency. An optional linear transform maps the section into the
    coolant-normal plane before measuring it.
    """
    n = len(inner.radii)
    area = np.zeros(n)
    moment = np.zeros((n, 2))
    lengths = []

    def integrate(points, derivatives):
        nonlocal area, moment
        if transform is not None:
            points = np.einsum("nij,nkj->nki", transform, points)
            derivatives = np.einsum("nij,nkj->nki", transform, derivatives)
        x, y = points[..., 0], points[..., 1]
        dx, dy = derivatives[..., 0], derivatives[..., 1]
        area += ((x * dy - y * dx) @ _WEIGHTS) / 2
        moment[:, 0] += (x*x * dy) @ _WEIGHTS / 2
        moment[:, 1] -= (y*y * dx) @ _WEIGHTS / 2
        lengths.append(np.linalg.norm(derivatives, axis=-1) @ _WEIGHTS)

    def arc_points(arc, fraction):
        angle = arc.start_angles[:, None] + (
            arc.end_angles - arc.start_angles
        )[:, None] * fraction
        directions = np.stack((np.cos(angle), np.sin(angle)), axis=-1)
        return arc.centers[:, None, :] + arc.radii[:, None, None] * directions

    for arc, next_arc in ((inner, outer), (outer, inner)):
        angle = arc.start_angles[:, None] + (
            arc.end_angles - arc.start_angles
        )[:, None] * _NODES
        derivatives = arc.radii[:, None, None] * (
            arc.end_angles - arc.start_angles
        )[:, None, None] * np.stack((-np.sin(angle), np.cos(angle)), axis=-1)
        integrate(arc_points(arc, _NODES), derivatives)
        start = arc_points(arc, np.array([1.0]))[:, 0]
        end = arc_points(next_arc, np.array([0.0]))[:, 0]
        delta = end - start
        integrate(start[:, None] + _NODES[None, :, None] * delta[:, None],
                  np.broadcast_to(delta[:, None], (n, len(_NODES), 2)))
    if np.any(~np.isfinite(area)) or np.any(np.abs(area) < 1e-30):
        raise ValueError("Channel boundary must enclose a finite, positive area")
    return np.abs(area), moment / area[:, None], np.stack(lengths, axis=1)


def rounded_gas_sector_area(center_distance, half_angle, liner=0.0):
    """Gas area inside one tube sector in its unprojected section plane.

    The main tube's exterior is tangent to both sector rays. A liner extends
    its gas-facing arc inward; the neighboring arcs meet on the sector rays.
    """
    d, a, t = np.broadcast_arrays(center_distance, half_angle, liner)
    radius = d * np.sin(a) + t
    root = np.sqrt(np.maximum(radius**2 - (d*np.sin(a))**2, 0.0))
    intersection = d * np.cos(a) - root
    angle = np.arctan2(intersection * np.sin(a), intersection * np.cos(a) - d)
    return d * radius * np.sin(angle) - radius**2 * (np.pi - angle)


def rounded_center_distance(target_area, half_angle, liner=0.0):
    """Solve the lobe-center distance for an allocated gas-sector area."""
    area, angle, thickness = np.broadcast_arrays(target_area, half_angle, liner)
    factor = np.sin(angle)*np.cos(angle) - (np.pi/2-angle)*np.sin(angle)**2
    if np.any(area <= 0) or np.any(factor <= 0) or np.any(thickness < 0):
        raise ValueError("Rounded placement requires positive area/pitch and nonnegative liners")
    result = np.sqrt(area/factor)
    for i in np.ndindex(result.shape):
        if thickness[i] == 0:
            continue
        lower = thickness[i] / (1 - np.sin(angle[i]))
        upper = lower + 2 * result[i]
        def residual(d):
            return float(rounded_gas_sector_area(d, angle[i], thickness[i]) - area[i])
        while residual(upper) < 0:
            upper *= 2
        result[i] = brentq(residual, lower, upper, xtol=1e-14)
    return result


def rounded_boundary(prof, base_height):
    radius = prof.sector_radius
    if radius is None:
        radius = prof.cross_section_base_point_path[:, 1]
    q_inner = radius * np.sin(prof.theta/2) - prof.t_wall
    q_outer = q_inner + base_height * np.sin(prof.theta/2)
    zero = np.zeros_like(radius)
    tangent = (np.pi + prof.theta) / 2
    return (
        SectionArc(np.column_stack((zero, zero)), q_inner, tangent, 2*np.pi-tangent),
        SectionArc(np.column_stack((base_height, zero)), q_outer, -tangent, tangent),
    )


def squared_boundary(prof, rib_width):
    radius = prof.cross_section_base_point_path[:, 1]
    projection = prof.projection
    if projection is None:
        projection = np.tile(np.eye(2), (len(radius), 1, 1))
    pitch = prof.azimuthal_pitch
    if pitch is None:
        pitch = prof.theta * projection[:, 1, 1]
    if np.any((pitch <= 0) | (pitch >= np.pi)):
        raise ValueError("Transverse construction pitch must satisfy 0 < pitch < pi")
    outer_radius = radius + prof.h * projection[:, 0, 0]
    # Rib faces run along the local radial direction, so a transverse section
    # shows ribs that point at the engine axis rather than leaning over.
    #
    # Two such faces are helicoids. Their separation within a transverse plane
    # is constant, but the material between them is that separation times the
    # cosine of the helix angle seen locally, and a fixed angular rate is a
    # steeper helix further out: the tightest material is at the channel roof,
    # not the floor. Open the azimuthal offset until that tightest gap is
    # exactly ``rib_width``.
    rate = prof.winding_rate
    twist = 0.0 if rate is None else (outer_radius * rate)**2
    offset = rib_width * np.sqrt(1.0 + twist) / 2

    def intersection(r, side):
        sine = -side * offset / r
        if np.any(np.abs(sine) >= 1):
            raise ValueError("Projected rib face does not intersect the channel radius")
        return side * pitch / 2 + np.arcsin(sine)

    inner_minus, inner_plus = intersection(radius, -1), intersection(radius, 1)
    outer_minus, outer_plus = intersection(outer_radius, -1), intersection(outer_radius, 1)
    if np.any(inner_plus <= inner_minus) or np.any(outer_plus <= outer_minus):
        raise ValueError("Projected rib leaves no open channel at the floor or roof")
    centers = np.column_stack((-radius, np.zeros_like(radius)))
    return (
        SectionArc(centers, radius, inner_minus, inner_plus),
        SectionArc(centers, outer_radius, outer_plus, outer_minus),
    )
