"""Shared geometric helpers for regenerative cooling-channel paths.

The regenerative model stores paths in cylindrical coordinates as
``[axial_position, radius, angle]``.  This module turns those lightweight
arrays into Cartesian section planes only when a consumer needs them.  The
path a consumer passes in decides what the planes mean: CAD sections and
visualization frame themselves on the cross-section base point, not on the
coolant centroid.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


__all__ = [
    "SectionPlanes",
    "build_axis_normal_section_planes",
    "build_section_planes",
    "cylindrical_centerline_to_cartesian",
    "estimate_centerline_tangents",
    "local_section_axes",
    "rotate_cylindrical_centerline",
]

_VECTOR_EPSILON = 1.0e-12


@dataclass(frozen=True)
class SectionPlanes:
    """Orthonormal planes located at the samples of a channel centerline.

    ``origins`` and the three axes use the Cartesian coordinate convention of
    PySkyFire: X is the engine axis, while Y and Z span a radial plane.
    ``normals`` follow either the centerline tangent or the engine axis.
    ``x_axes`` point outward and form the zero-angle direction of a section.
    """

    origins: np.ndarray
    normals: np.ndarray
    x_axes: np.ndarray

    def __post_init__(self) -> None:
        arrays = tuple(
            np.asarray(array, dtype=float)
            for array in (self.origins, self.normals, self.x_axes)
        )
        shapes = [np.asarray(array).shape for array in arrays]
        if not shapes or any(shape != shapes[0] for shape in shapes):
            raise ValueError("Section-plane arrays must all have the same shape")
        if len(shapes[0]) != 2 or shapes[0][1] != 3:
            raise ValueError("Section-plane arrays must have shape (N, 3)")
        if shapes[0][0] == 0:
            raise ValueError("At least one section plane is required")
        if any(not np.all(np.isfinite(array)) for array in arrays):
            raise ValueError("Section-plane arrays must contain only finite values")
        object.__setattr__(self, "origins", arrays[0])
        object.__setattr__(self, "normals", arrays[1])
        object.__setattr__(self, "x_axes", arrays[2])

    def __len__(self) -> int:
        return len(self.origins)

    @property
    def y_axes(self) -> np.ndarray:
        """Return each plane's second in-plane axis."""

        return np.cross(self.normals, self.x_axes)


def _validate_centerline(centerline: np.ndarray) -> np.ndarray:
    centerline = np.asarray(centerline, dtype=float)
    if centerline.ndim != 2 or centerline.shape[1] != 3:
        raise ValueError("A channel centerline must have shape (N, 3)")
    if len(centerline) == 0:
        raise ValueError("A channel centerline must contain at least one point")
    if not np.all(np.isfinite(centerline)):
        raise ValueError("A channel centerline must contain only finite values")
    return centerline


def cylindrical_centerline_to_cartesian(centerline: np.ndarray) -> np.ndarray:
    """Convert a cylindrical ``[x, r, theta]`` centerline to ``[x, y, z]``."""

    centerline = _validate_centerline(centerline)
    axial = centerline[:, 0]
    radius = centerline[:, 1]
    angle = centerline[:, 2]
    return np.column_stack(
        (
            axial,
            radius * np.cos(angle),
            radius * np.sin(angle),
        )
    )


def rotate_cylindrical_centerline(
    centerline: np.ndarray,
    angle_rad: float,
) -> np.ndarray:
    """Return a copy of a cylindrical centerline rotated about the engine axis."""

    centerline = _validate_centerline(centerline)
    if not np.isfinite(angle_rad):
        raise ValueError("Rotation angle must be finite")
    rotated = centerline.copy()
    rotated[:, 2] += angle_rad
    return rotated


def estimate_centerline_tangents(points: np.ndarray) -> np.ndarray:
    """Estimate Cartesian polyline tangents with centered differences.
    """

    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("Cartesian centerline points must have shape (N, 3)")
    if len(points) == 0:
        raise ValueError("At least one Cartesian point is required")
    if len(points) == 1:
        return np.array([[1.0, 0.0, 0.0]])

    tangents = np.gradient(points, axis=0, edge_order=min(2, len(points)-1))

    magnitudes = np.linalg.norm(tangents, axis=1)
    if np.any(magnitudes <= _VECTOR_EPSILON):
        raise ValueError("Consecutive centerline points must not coincide")
    return tangents / magnitudes[:, None]


def local_section_axes(
    point: np.ndarray, tangent: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Construct radial and circumferential axes perpendicular to ``tangent``."""

    point = np.asarray(point, dtype=float)
    tangent = np.asarray(tangent, dtype=float)
    if point.shape != (3,) or tangent.shape != (3,):
        raise ValueError("A section point and tangent must each contain three values")

    tangent_norm = np.linalg.norm(tangent)
    if tangent_norm <= _VECTOR_EPSILON:
        raise ValueError("A section-plane tangent must be nonzero")
    tangent = tangent / tangent_norm

    radial = np.array([0.0, point[1], point[2]], dtype=float)
    radial_norm = np.linalg.norm(radial)
    if radial_norm <= _VECTOR_EPSILON:
        radial = np.array([0.0, 1.0, 0.0])
    else:
        radial /= radial_norm

    radial_component = np.dot(radial, tangent)
    radial = tangent[0] * radial - radial_component * np.array([1.0, 0.0, 0.0])
    magnitude = np.linalg.norm(radial)
    if magnitude <= _VECTOR_EPSILON:
        raise ValueError("Section tangent must have a nonzero meridional component")
    radial /= magnitude

    circumferential = np.cross(tangent, radial)
    circumferential /= np.linalg.norm(circumferential)
    return radial, circumferential


def build_section_planes(centerline: np.ndarray) -> SectionPlanes:
    """Build section planes normal to a cylindrical path."""

    origins = cylindrical_centerline_to_cartesian(centerline)
    path = _validate_centerline(centerline)
    x, radius, angle = path.T
    if len(path) == 1:
        dr = da = np.zeros(1)
    else:
        if np.any(np.diff(x) <= 0):
            raise ValueError("Path axial stations must increase strictly")
        dr = np.gradient(radius, x, edge_order=min(2, len(x)-1))
        da = np.gradient(np.unwrap(angle), x, edge_order=min(2, len(x)-1))
    normals = np.column_stack((np.ones(len(x)),
        dr*np.cos(angle)-radius*da*np.sin(angle),
        dr*np.sin(angle)+radius*da*np.cos(angle)))
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    x_axes = np.empty_like(origins)
    for index, (point, tangent) in enumerate(zip(origins, normals, strict=True)):
        x_axes[index], _ = local_section_axes(point, tangent)
    return SectionPlanes(origins=origins, normals=normals, x_axes=x_axes)


def build_axis_normal_section_planes(centerline: np.ndarray) -> SectionPlanes:
    """Build transverse planes with origins on the supplied path.

    The normal is always engine +X, independent of coolant direction. The
    first in-plane axis points radially away from the engine axis, so the
    engine axis has local coordinates ``(-r, 0)`` in every plane.
    """
    centerline = _validate_centerline(centerline)
    angle = centerline[:, 2]
    return SectionPlanes(
        origins=cylindrical_centerline_to_cartesian(centerline),
        normals=np.tile([1.0, 0.0, 0.0], (len(centerline), 1)),
        x_axes=np.column_stack((np.zeros_like(angle), np.cos(angle), np.sin(angle))),
    )
