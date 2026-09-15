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

    def evaluate(self, fraction):
        angle = self.start_angles[:, None] + (
            self.end_angles - self.start_angles
        )[:, None] * np.asarray(fraction)
        return self.centers[:, None, :] + self.radii[:, None, None] * np.stack(
            (np.cos(angle), np.sin(angle)), axis=-1
        )

    def derivative(self, fraction):
        angle = self.start_angles[:, None] + (
            self.end_angles - self.start_angles
        )[:, None] * np.asarray(fraction)
        return self.radii[:, None, None] * (
            self.end_angles - self.start_angles
        )[:, None, None] * np.stack((-np.sin(angle), np.cos(angle)), axis=-1)


@dataclass(frozen=True)
class SectionEllipseArc:
    """Elliptical boundary, with the same sampling API as a circular arc.

    ``semi_x`` and ``semi_y`` are the semi-axes along the section plane's own
    two axes, so equal semi-axes reproduce a :class:`SectionArc` exactly. The
    angles are the ellipse's parameter, not the polar angle to the point: a
    point sits at ``center + (semi_x*cos(t), semi_y*sin(t))``.
    """

    centers: np.ndarray
    semi_x: np.ndarray
    semi_y: np.ndarray
    start_angles: np.ndarray
    end_angles: np.ndarray

    def _angle(self, fraction):
        return self.start_angles[:, None] + (
            self.end_angles - self.start_angles
        )[:, None] * np.asarray(fraction)

    def points(self, index, count):
        angles = np.linspace(self.start_angles[index], self.end_angles[index], count)
        return self.centers[index] + np.column_stack(
            (self.semi_x[index] * np.cos(angles), self.semi_y[index] * np.sin(angles))
        )

    def evaluate(self, fraction):
        angle = self._angle(fraction)
        return self.centers[:, None, :] + np.stack(
            (self.semi_x[:, None] * np.cos(angle),
             self.semi_y[:, None] * np.sin(angle)), axis=-1
        )

    def derivative(self, fraction):
        angle = self._angle(fraction)
        sweep = (self.end_angles - self.start_angles)[:, None]
        return sweep[..., None] * np.stack(
            (-self.semi_x[:, None] * np.sin(angle),
             self.semi_y[:, None] * np.cos(angle)), axis=-1
        )


@dataclass(frozen=True)
class SectionLine:
    """Directed straight boundary, with the same sampling API as an arc."""

    starts: np.ndarray
    ends: np.ndarray

    def evaluate(self, fraction):
        return self.starts[:, None, :] + np.asarray(fraction)[None, :, None] * (
            self.ends - self.starts
        )[:, None, :]

    def derivative(self, fraction):
        return np.broadcast_to((self.ends - self.starts)[:, None, :],
                               (len(self.starts), len(fraction), 2))

    def points(self, index, count):
        return self.starts[index] + np.linspace(0, 1, count)[:, None] * (
            self.ends[index] - self.starts[index]
        )


def boundary_surfaces(floor, roof):
    """Four directed surfaces, in floor/right/roof/left integration order."""
    return (floor, SectionLine(floor.evaluate([1])[:, 0], roof.evaluate([0])[:, 0]),
            roof, SectionLine(roof.evaluate([1])[:, 0], floor.evaluate([0])[:, 0]))


def surface_lengths(surface, fraction=1.0, transform=None):
    """Length of a centred fraction of physical surface length [m]."""
    fraction = np.broadcast_to(np.asarray(fraction, float),
                               (len(surface.evaluate([0])),))
    # Fractions refer to physical length even when projection makes an arc
    # elliptical. Rounded angular exposure is converted in its circular plane.
    deriv = surface.derivative(_NODES)
    if transform is not None:
        deriv = np.einsum("nij,nkj->nki", transform, deriv)
    return (np.linalg.norm(deriv, axis=-1) @ _WEIGHTS) * fraction


_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(24)
_NODES = (_NODES + 1.0) / 2.0
_WEIGHTS = _WEIGHTS / 2.0


def arc_boundary_metrics(inner, outer, transform=None):
    """Area, centroid, and edge lengths of arc/line/arc/line boundaries.

    Green's theorem integrates the actual arcs, without a polygon-resolution
    dependency. An optional linear transform maps the section into the
    coolant-normal plane before measuring it.
    """
    n = len(inner.evaluate([0]))
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

    for surface in boundary_surfaces(inner, outer):
        integrate(surface.evaluate(_NODES), surface.derivative(_NODES))
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
    return rounded_width_boundary(2*q_inner*np.cos(prof.theta/2),
                                  2*q_outer*np.cos(prof.theta/2),
                                  base_height, height_input="base",
                                  side_angle=prof.theta/2)


def rounded_width_boundary(floor_width, roof_width, height, *, height_input="full",
                           side_angle=None):
    """Two circular lobes joined tangentially, specified by endpoint chords.

    With full height, choose the solution with the smaller sidewall angle
    (larger separation of lobe centres). Sector construction supplies its
    known angle explicitly, including the zero-length-side circular limit.
    """
    lower, upper, height = np.broadcast_arrays(floor_width, roof_width, height)
    if np.any(~np.isfinite(lower + upper + height)) or np.any(lower <= 0) or np.any(upper <= 0):
        raise ValueError("Rounded widths must be finite and positive; height must be finite")
    delta = (upper - lower) / 2
    if side_angle is not None:
        angle = np.broadcast_to(side_angle, height.shape)
        distance = height
    elif height_input == "base":
        if np.any(height < 0) or np.any((height == 0) & (delta != 0)):
            raise ValueError("Unequal rounded widths require positive base height")
        angle = np.arctan2(delta, height)
        distance = height
    else:
        angle = np.zeros_like(height)
        for index in np.ndindex(height.shape):
            if delta[index] == 0:
                if height[index] < lower[index]:
                    raise ValueError("Rounded height is below the geometric minimum")
                continue
            mean = (lower[index] + upper[index]) / 2
            difference = abs(delta[index])
            def required(a):
                return (mean + difference / np.sin(a)) / np.cos(a)
            ratio = difference/mean
            sine = brentq(lambda z: z**3+2*ratio*z*z-ratio, 0, 1/np.sqrt(2), xtol=1e-15)
            minimum_angle = np.arcsin(sine)
            minimum_height = required(minimum_angle)
            if height[index] < minimum_height * (1 - 1e-12):
                raise ValueError("Rounded height is below the geometric minimum for these widths")
            a = minimum_angle if height[index] <= minimum_height else brentq(
                lambda a: required(a) - height[index], 1e-15, minimum_angle, xtol=1e-15)
            angle[index] = np.copysign(a, delta[index])
        distance = height - (lower + upper) / (2*np.cos(angle))
    q_inner = lower / (2*np.cos(angle))
    q_outer = upper / (2*np.cos(angle))
    zero = np.zeros_like(height)
    tangent = np.pi/2 + angle
    return (
        SectionArc(np.column_stack((zero, zero)), q_inner, tangent, 2*np.pi-tangent),
        SectionArc(np.column_stack((distance, zero)), q_outer, -tangent, tangent),
    )


def squared_boundary(prof, sidewall_half_thickness):
    radius = prof.cross_section_base_point_path[:, 1]
    if np.any(~np.isfinite(prof.h)) or np.any(prof.h <= 0):
        raise ValueError("Squared channel height must be finite and positive")
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
    # exactly twice ``sidewall_half_thickness``.
    rate = prof.winding_rate
    twist = 0.0 if rate is None else (outer_radius * rate)**2
    offset = sidewall_half_thickness * np.sqrt(1.0 + twist)

    def intersection(r, side):
        sine = -side * offset / r
        if np.any(np.abs(sine) >= 1):
            raise ValueError("Projected rib face does not intersect the channel radius")
        return side * pitch / 2 + np.arcsin(sine)

    inner_minus, inner_plus = intersection(radius, -1), intersection(radius, 1)
    outer_minus, outer_plus = intersection(outer_radius, -1), intersection(outer_radius, 1)
    if np.any(inner_plus <= inner_minus) or np.any(outer_plus <= outer_minus):
        raise ValueError("Projected rib leaves no open channel at the floor or roof")
    # Convert occupied-sector geometry into internal endpoint widths, then
    # use the same construction as an independently dimensioned exchanger.
    lower = 2*radius*np.sin(inner_plus)
    upper = 2*outer_radius*np.sin(outer_plus)
    return squared_width_boundary(lower, upper, outer_radius-radius, radius, "arc", "arc")


def squared_width_boundary(floor_width, roof_width, height, radius, floor, roof):
    """Internal chord widths, with floor/roof midpoints separated by height.

    The floor arc has radius ``radius`` and the roof arc ``radius + height``.
    Infinite radius represents a line exactly, without evaluating inf * 0.
    """
    lower, upper, height, radius = np.broadcast_arrays(floor_width, roof_width, height, radius)
    if np.any(~np.isfinite(lower + upper + height)) or np.any(lower <= 0) or np.any(upper <= 0) or np.any(height <= 0):
        raise ValueError("Squared widths and height must be finite and positive")
    if np.any(np.isnan(radius)) or np.any(radius <= 0):
        raise ValueError("Radius must be positive or infinity")
    curves = []
    for width, level, r, kind, sign in ((lower, np.zeros_like(height), radius, floor, 1),
                                      (upper, height, radius+height, roof, -1)):
        if kind == "line" or np.all(np.isinf(r)):
            start = np.column_stack((level, -sign*width/2))
            end = np.column_stack((level, sign*width/2))
            curves.append(SectionLine(start, end))
        else:
            if np.any(~np.isfinite(r)):
                raise ValueError("Use separate profiles for finite and infinite arc radii")
            if np.any(width >= 2*r):
                raise ValueError("Arc chord width must be smaller than its diameter")
            half = np.arcsin(width/(2*r))
            curves.append(SectionArc(np.column_stack((level-r, np.zeros_like(level))),
                                     r, -sign*half, sign*half))
    # Floor must remain below the roof at both sidewall endpoints.
    if np.any(curves[1].evaluate([0])[:, 0, 0] <= curves[0].evaluate([1])[:, 0, 0]):
        raise ValueError("Floor and roof intersect; increase height or reduce curvature")
    return tuple(curves)


def offset_squared_boundary(floor, roof, floor_thickness, roof_thickness, side_thickness):
    """Offset each surface normally and intersect it with the offset sides.

    Thicknesses are normal distances in this section plane. Corners are
    intersections (mitres), so unequal floor and roof thicknesses are allowed.
    """
    _, right, _, left = boundary_surfaces(floor, roof)
    sides = []
    for side in (right, left):
        vector = side.ends-side.starts
        norm = np.linalg.norm(vector, axis=1)
        normal = np.column_stack((-vector[:, 1], vector[:, 0]))/norm[:, None]
        sides.append(SectionLine(side.starts + side_thickness[:, None]*normal,
                                 side.ends + side_thickness[:, None]*normal))

    def intersection(curve, side, near):
        start, vector = side.starts, side.ends-side.starts
        if isinstance(curve, SectionLine):
            other = curve.ends-curve.starts
            cross = lambda a, b: a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]
            divisor = cross(vector, other)
            if np.any(np.abs(divisor) < 1e-30):
                raise ValueError("Offset surfaces do not form a finite corner")
            t = cross(curve.starts-start, other)/divisor
        else:
            delta = start-curve.centers
            a = np.sum(vector**2, axis=1)
            b = 2*np.sum(delta*vector, axis=1)
            c = np.sum(delta**2, axis=1)-curve.radii**2
            discriminant = b*b-4*a*c
            if np.any(discriminant < 0):
                raise ValueError("Offset sidewall does not intersect the end arc")
            root = np.sqrt(discriminant)
            t0, t1 = (-b-root)/(2*a), (-b+root)/(2*a)
            p0, p1 = start+t0[:, None]*vector, start+t1[:, None]*vector
            t = np.where(np.linalg.norm(p0-near, axis=1) < np.linalg.norm(p1-near, axis=1), t0, t1)
        return start+t[:, None]*vector

    result = []
    for curve, thickness, sign, first, second in (
        (floor, floor_thickness, -1, sides[1], sides[0]),
        (roof, roof_thickness, 1, sides[0], sides[1]),
    ):
        if isinstance(curve, SectionArc):
            radii = curve.radii+sign*thickness
            if np.any(radii <= 0):
                raise ValueError("Floor thickness must be smaller than the floor radius")
            offset = SectionArc(curve.centers, radii, curve.start_angles, curve.end_angles)
        else:
            shift = np.column_stack((sign*thickness, np.zeros_like(thickness)))
            offset = SectionLine(curve.starts+shift, curve.ends+shift)
        start = intersection(offset, first, curve.evaluate([0])[:, 0])
        end = intersection(offset, second, curve.evaluate([1])[:, 0])
        if isinstance(offset, SectionLine):
            result.append(SectionLine(start, end))
        else:
            a = np.arctan2((start-offset.centers)[:, 1], (start-offset.centers)[:, 0])
            b = np.arctan2((end-offset.centers)[:, 1], (end-offset.centers)[:, 0])
            result.append(SectionArc(offset.centers, offset.radii, a, b))
    return tuple(result)
