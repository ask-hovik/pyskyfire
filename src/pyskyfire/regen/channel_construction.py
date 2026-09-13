"""Resolve nominal paths into physical section anchors, arcs and centroids.

The normal section is projected locally onto an engine-axis-normal plane to
allocate gas area. This is exact for prismatic sections; a curved CAD loft is
an approximation to the sequence of these local sections.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .channel_helpers import SectionPlanes, build_axis_normal_section_planes
from .cross_section import CrossSectionRounded, CrossSectionRoundedInternal, SectionProfiles
from .section_geometry import SectionArc, rounded_center_distance, rounded_gas_sector_area, squared_boundary

__all__ = ["ConstructionArc", "ChannelConstruction", "build_channel_construction"]
ConstructionArc = SectionArc


@dataclass(frozen=True)
class ChannelConstruction:
    """One circuit's resolved section geometry along a nominal winding.

    ``planes`` are what the section boundary is drawn in -- transverse for a
    squared section, flow-normal for a rounded one -- while ``flow_planes`` are
    always flow-normal and are what the centroid is measured in. ``profiles``
    carries the resolved dimensions, including the base-point path the section
    is anchored on, and ``gas_area`` the combustion-gas area each channel's
    sector keeps.
    """

    section_type: str
    reference_path: np.ndarray
    planes: SectionPlanes
    flow_planes: SectionPlanes
    profiles: SectionProfiles
    inner_arc: SectionArc
    outer_arc: SectionArc
    circular: np.ndarray
    alpha: np.ndarray
    gamma: np.ndarray
    pitch: np.ndarray
    gas_area: np.ndarray
    rib_width: np.ndarray | None = None

    def polygon(self, index, n_arc_points=8):
        if n_arc_points < 2:
            raise ValueError("n_arc_points must be at least 2")
        uv = np.vstack((self.inner_arc.points(index, n_arc_points),
                        self.outer_arc.points(index, n_arc_points)))
        return (self.planes.origins[index]
                + uv[:, :1]*self.planes.x_axes[index]
                + uv[:, 1:]*self.planes.y_axes[index])

    def centroid_path(self, cross_section):
        """Return the section centroid in world cylindrical coordinates."""
        if hasattr(cross_section, "centroid"):
            uv = cross_section.centroid(self.profiles)
        else:
            uv = np.column_stack((self.profiles.h/2, np.zeros(len(self.profiles.h))))
        xyz = (self.flow_planes.origins + uv[:, :1]*self.flow_planes.x_axes
               + uv[:, 1:]*self.flow_planes.y_axes)
        base = self.profiles.cross_section_base_point_path
        angle = base[:, 2] + np.arctan2(
            xyz[:, 2]*np.cos(base[:, 2])-xyz[:, 1]*np.sin(base[:, 2]),
            xyz[:, 1]*np.cos(base[:, 2])+xyz[:, 2]*np.sin(base[:, 2]),
        )
        return np.column_stack((xyz[:, 0], np.hypot(xyz[:, 1], xyz[:, 2]), angle))


def _values(value, size):
    return np.broadcast_to(np.asarray(value, dtype=float), (size,)).copy()


def _derivative(values, x):
    if len(x) == 1:
        return np.zeros_like(x)
    if np.any(np.diff(x) <= 0):
        raise ValueError("Reference path stations must increase along engine X")
    return np.gradient(values, x, edge_order=min(2, len(x)-1))


def build_channel_construction(circuit, reference_path=None):
    """Build sections from the ideal-contour winding, including rotated seeds.

    Every angular sector receives its share of the target gas area. Rounded
    lobes are tangent to that sector in the locally projected geometry; liners
    extend the gas-facing arc inward. No translation follows lobe sizing.
    """
    if reference_path is None:
        reference_path = circuit.reference_path
    reference_path = np.asarray(reference_path, dtype=float)
    if reference_path.ndim != 2 or reference_path.shape[1] != 3 or not len(reference_path):
        raise ValueError("reference_path must have shape (N, 3)")
    x, radius, angle = reference_path.T
    if np.any(~np.isfinite(reference_path)) or np.any(radius <= 0):
        raise ValueError("Reference path must be finite with positive radii")
    n = len(x)
    domain = circuit.x_domain
    nominal = getattr(circuit, "reference_path", reference_path)
    slope = _derivative(radius, x)
    contour = getattr(circuit, "contour", None)
    if contour is not None and hasattr(contour, "dr_dx"):
        slope = np.array([float(contour.dr_dx(float(value))) for value in x])
    alpha = np.arctan(slope)
    stored_gamma = getattr(circuit, "helix_angle_vals", None)
    if stored_gamma is None:
        gamma = np.arctan(radius*_derivative(angle, x)/np.sqrt(1+slope**2))
        gamma_nominal = gamma
    else:
        gamma_nominal = np.interp(x, domain, stored_gamma)
        if np.any(~np.isfinite(gamma_nominal)) or np.any(np.abs(gamma_nominal) >= np.pi/2-1e-8):
            raise ValueError("Construction helix angle must be finite and below 90 degrees")
        delta = angle - np.interp(x, nominal[:, 0], nominal[:, 2])
        gamma = np.arctan(np.tan(gamma_nominal) + radius*_derivative(delta, x)/np.sqrt(1+slope**2))
    ca, sa, cg, sg = np.cos(alpha), np.sin(alpha), np.cos(gamma), np.sin(gamma)
    if np.any(cg < 1e-8):
        raise ValueError("Construction helix angle is too close to 90 degrees")
    # Normal-plane (meridional-normal, across-flow) to transverse (radial, azimuthal).
    projection = np.zeros((n, 2, 2))
    projection[:, 0, 0] = 1/ca
    projection[:, 1, 0] = np.tan(alpha)*np.tan(gamma)
    projection[:, 1, 1] = 1/cg
    inverse = np.linalg.inv(projection)
    azimuthal = getattr(circuit, "channel_angular_width", None)
    if azimuthal is None:
        pitch = np.interp(x, domain, circuit.channel_local_sector)/np.cos(gamma_nominal)
    else:
        pitch = np.interp(x, domain, azimuthal)
    if np.any((pitch <= 0) | (pitch >= np.pi)):
        raise ValueError("Transverse construction pitch must satisfy 0 < pitch < pi")
    height = np.interp(x, domain, circuit.channel_heights)
    main = _values(circuit.channel_wall_thickness(x), n)
    liner = _values(circuit.hot_side_layer_thickness(x), n) if hasattr(circuit, "hot_side_layer_thickness") else np.zeros(n)
    if np.any(~np.isfinite(height)) or np.any(height < 0) or np.any(~np.isfinite(main+liner)) or np.any(main < 0) or np.any(liner < 0):
        raise ValueError("Section dimensions must be finite and nonnegative")
    e_r = np.column_stack((np.zeros(n), np.cos(angle), np.sin(angle)))
    e_phi = np.column_stack((np.zeros(n), -np.sin(angle), np.cos(angle)))
    axial = np.tile([1., 0., 0.], (n, 1))
    normals = ca[:, None]*cg[:, None]*axial + sa[:, None]*cg[:, None]*e_r + sg[:, None]*e_phi
    u_axis = -sa[:, None]*axial + ca[:, None]*e_r
    v_axis = np.cross(normals, u_axis)
    rounded = isinstance(circuit.cross_section, (CrossSectionRounded, CrossSectionRoundedInternal))
    gas_area = pitch*radius**2/2
    sector_radius = None
    if rounded:
        b = pitch/2
        rays = np.stack((np.column_stack((np.cos(b), -np.sin(b))),
                         np.column_stack((np.cos(b), np.sin(b)))), axis=1)
        local = np.einsum("nij,nkj->nki", inverse, rays)
        local /= np.linalg.norm(local, axis=2)[:, :, None]
        bisector = local.sum(axis=1)
        bisector /= np.linalg.norm(bisector, axis=1)[:, None]
        a = np.arccos(np.clip(np.einsum("ni,ni->n", local[:, 0], local[:, 1]), -1, 1))/2
        sector_radius = rounded_center_distance(gas_area/np.linalg.det(projection), a, liner)
        center = np.einsum("nij,nj->ni", projection, sector_radius[:, None]*bisector)
        x_axes = bisector[:, :1]*u_axis + bisector[:, 1:]*v_axis
        theta = 2*a
        gas_area = rounded_gas_sector_area(sector_radius, a, liner)*np.linalg.det(projection)
    else:
        center = np.column_stack((radius+(main+liner)/ca, np.zeros(n)))
        x_axes = u_axis
        theta = pitch*cg
    origins = x[:, None]*axial + center[:, :1]*e_r + center[:, 1:]*e_phi
    base_path = np.column_stack((x, np.linalg.norm(center, axis=1), angle+np.arctan2(center[:, 1], center[:, 0])))
    flow_planes = SectionPlanes(origins, normals, x_axes)
    context_builder = getattr(circuit, "_design_input_context", None)
    context = context_builder(x) if context_builder else {"contour": contour}
    context.update(gamma=gamma_nominal.copy(), t_wall_total=main+liner, t_hot_layer=liner)
    profiles = SectionProfiles(
        h=height, theta=theta, t_wall=main,
        cross_section_base_point_path=base_path, sector_radius=sector_radius,
        projection=None if rounded else projection, azimuthal_pitch=pitch,
        # The helix angle is measured on the reference path, so that is the
        # radius its tangent belongs to; every other radius shares the rate.
        winding_rate=np.tan(gamma)/(ca*radius),
        **context,
    )
    cs = circuit.cross_section
    if rounded:
        inner, outer = cs.boundary(profiles)
        circular = cs.base_height(profiles) == 0
        planes = flow_planes
        rib_width = None
    else:
        if np.any(height <= 0):
            raise ValueError("Squared channel heights must be positive")
        rib_width = cs.rib_width(profiles) if hasattr(cs, "rib_width") else np.zeros(n)
        inner, outer = squared_boundary(profiles, rib_width)
        circular = np.zeros(n, dtype=bool)
        planes = build_axis_normal_section_planes(base_path)
    return ChannelConstruction(
        "rounded" if rounded else "squared", reference_path.copy(), planes,
        flow_planes, profiles, inner, outer, circular, alpha, gamma, pitch, gas_area, rib_width,
    )
