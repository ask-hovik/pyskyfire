"""Engineering objects adapted into local-coordinate scene hierarchies."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from pyskyfire.regen.channel_construction import build_channel_construction
from pyskyfire.regen.cad_geometry import transverse_construction
from pyskyfire.regen.channel_layout import ChannelLayout, resolve_analytical_layout, resolve_numerical_layout
from .scene import Anchor, FieldData, MeshData, SceneNode, Transform

__all__ = ["ThrustChamberComponent", "ImpellerComponent", "PipeComponent"]

_COLORS = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

# Attribute, native axial grid, base SI label. No unit conversion is performed.
_REGEN_FIELDS = {
    "p_static": ("x_coolant", "Pa", "Static pressure"),
    "p_stagnation": ("x_coolant", "Pa", "Stagnation pressure"),
    "T_static": ("x_coolant", "K", "Static coolant temperature"),
    "T_stagnation": ("x_coolant", "K", "Stagnation coolant temperature"),
    "velocity": ("x_coolant", "m/s", "Coolant velocity"),
    "coolant_enthalpy": ("x_coolant", "J/kg", "Coolant enthalpy"),
    "coolant_quality": ("x_coolant", "1", "Vapour mass fraction"),
    "coolant_T_sat": ("x_coolant", "K", "Saturation temperature"),
    "qpp_hot": ("x_heat_flux", "W/m²", "Hot-side heat flux"),
    "dQ_dA": ("x_heat_flux", "W/m²", "Heat flux"),
    "T_hot_wall": ("x_wall", "K", "Hot wall temperature"),
    "T_coolant_wall": ("x_wall", "K", "Coolant wall temperature"),
}


def _surface(xyz, *, closed=False):
    """Triangulate an ordered surface grid, retaining the supplied stations."""
    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim != 3 or xyz.shape[2] != 3 or min(xyz.shape[:2]) < 2:
        raise ValueError("Surface grid must have shape (M,N,3), with M,N >= 2")
    m, n, _ = xyz.shape
    faces = []
    for i in range(m - 1):
        for j in range(n if closed else n - 1):
            k = (j + 1) % n
            a, b, c, d = i*n+j, i*n+k, (i+1)*n+k, (i+1)*n+j
            faces.extend(((a, b, c), (a, c, d)))
    return MeshData(xyz.reshape(-1, 3), np.asarray(faces, dtype=np.int64))


def _end_anchors(centers, direction=1):
    if direction >= 0:
        return {"inlet": Anchor.from_direction(centers[0], centers[0]-centers[1]),
                "outlet": Anchor.from_direction(centers[-1], centers[-1]-centers[-2])}
    return {"inlet": Anchor.from_direction(centers[-1], centers[-1]-centers[-2]),
            "outlet": Anchor.from_direction(centers[0], centers[0]-centers[1])}


def _results_by_circuit(data, chamber):
    if data is None:
        return {}
    if hasattr(data, "circuit_index"):
        values = [data]
    else:
        values = data.values() if isinstance(data, Mapping) else data
    result_map = {}
    for result in values:
        index = result.circuit_index
        if not isinstance(index, (int, np.integer)) or not 0 <= index < len(chamber.cooling_circuits):
            raise ValueError(f"Invalid result circuit_index: {index}")
        if index in result_map:
            raise ValueError(f"Duplicate results for circuit {index}")
        if result.circuit_name != chamber.cooling_circuits[index].name:
            raise ValueError(f"Result name does not match circuit {index}")
        result_map[index] = result
    return result_map


def _regen_fields(result, x, perimeter_count):
    fields = {}
    for name, (grid_name, units, label) in _REGEN_FIELDS.items():
        if name in ("T_hot_wall", "T_coolant_wall"):
            wall = getattr(result, "T", None)
            if wall is not None and (np.asarray(wall).ndim != 2 or np.asarray(wall).shape[1] < 3):
                raise ValueError("Wall temperatures must contain coolant, coolant-wall, and hot-wall columns")
            values = None if wall is None else np.asarray(wall)[:, -1 if name == "T_hot_wall" else 1]
        else:
            values = getattr(result, name, None)
        if values is None:
            continue
        grid = getattr(result, grid_name, None)
        grid = np.asarray(result.x if grid is None else grid, dtype=float)
        values = np.asarray(values, dtype=float)
        if grid.ndim != 1 or values.shape != grid.shape or len(grid) < 2 or not np.isfinite(grid).all():
            raise ValueError(f"{name} must match its native {grid_name} grid")
        if np.isinf(values).any():
            raise ValueError(f"{name} contains infinite values")
        if np.all(np.diff(grid) < 0):
            grid, values = grid[::-1], values[::-1]
        elif not np.all(np.diff(grid) > 0):
            raise ValueError(f"{grid_name} must be strictly monotone")
        # Fields outside the solved domain remain unavailable. Do not extrapolate
        # or bridge NaN gaps (e.g. vapour quality outside the boiling dome).
        sampled = np.interp(x, grid, values, left=np.nan, right=np.nan)
        fields[name] = FieldData(np.repeat(sampled, perimeter_count), units, label)
    return fields


@dataclass
class ThrustChamberComponent:
    """Cooling-channel surfaces with optional RegenResult fields.

    Each layout seed is constructed once. Instances share its immutable mesh
    and fields and carry separate rotations and identities. No hot wall is
    constructed. ``circuit_N.inlet/outlet`` anchors refer to the first seed's
    first channel, at its end-section polygon centre, not a manifold port.
    Every channel also exposes its own inlet/outlet anchors.
    """

    thrust_chamber: object
    cooling_data: object = None
    stride: int = 2
    n_arc_pts: int = 8
    layout: str | ChannelLayout = "numerical"
    circuit_colors: object = None

    def build(self):
        if not isinstance(self.stride, int) or self.stride < 1:
            raise ValueError("stride must be an integer of at least 1")
        if not isinstance(self.n_arc_pts, int) or self.n_arc_pts < 2:
            raise ValueError("n_arc_pts must be an integer of at least 2")
        chamber = self.thrust_chamber
        if self.circuit_colors is not None and len(self.circuit_colors) < len(chamber.cooling_circuits):
            raise ValueError("circuit_colors must provide a color for every cooling circuit")
        results = _results_by_circuit(self.cooling_data, chamber)
        if isinstance(self.layout, ChannelLayout):
            layout = self.layout
        elif self.layout == "analytical":
            layout = resolve_analytical_layout(chamber)
        elif self.layout == "numerical":
            layout = resolve_numerical_layout(chamber)
        else:
            raise ValueError("layout must be 'numerical', 'analytical', or ChannelLayout")
        palette = _COLORS if self.circuit_colors is None else self.circuit_colors
        root = SceneNode("chamber", metadata={"kind": "thrust_chamber"})
        groups = {}
        for seed in layout.channels:
            ci = seed.circuit_index
            circuit = chamber.cooling_circuits[ci]
            if ci not in groups:
                groups[ci] = SceneNode(f"circuit_{ci}", metadata={"kind": "cooling_circuit", "label": circuit.name, "circuit_index": ci})
                root.children.append(groups[ci])
            construction = transverse_construction(build_channel_construction(circuit, seed.reference_path))
            count = len(construction.planes)
            if count < 2:
                raise ValueError("Cooling-channel visualization requires at least two stations")
            indices = list(range(0, count, self.stride))
            if indices[-1] != count - 1:
                indices.append(count - 1)
            sections = np.asarray([construction.polygon(i, self.n_arc_pts) for i in indices])
            mesh = _surface(sections, closed=True)
            fields = _regen_fields(results[ci], sections[:, 0, 0], sections.shape[1]) if ci in results else {}
            # Full endpoint sampling avoids making anchor direction depend on stride.
            endpoint_indices = [0, 1, count-2, count-1]
            centers = np.array([construction.polygon(i, self.n_arc_pts).mean(axis=0) for i in endpoint_indices])
            anchors = _end_anchors(centers, getattr(circuit, "direction", 1))
            if not groups[ci].anchors:
                groups[ci].anchors.update(anchors)
                root.anchors.update({f"circuit_{ci}.{key}": value for key, value in anchors.items()})
            for instance in range(seed.pattern_count):
                groups[ci].children.append(SceneNode(
                    f"seed_{seed.seed_index}_channel_{instance}", mesh=mesh, fields=dict(fields),
                    anchors=dict(anchors), color=palette[ci % len(palette)],
                    transform=Transform(rotation=(instance*seed.pattern_angle_rad, 0, 0)),
                    metadata={"kind": "cooling_channel", "circuit_index": ci,
                              "seed_index": seed.seed_index, "instance_index": instance},
                ))
        return root


@dataclass
class ImpellerComponent:
    """Display prepared impeller.geometry arrays, with the existing local x shaft axis.

    Optional ``results`` maps part names (hub, shroud, blade_N) to dictionaries
    of mesh-aligned FieldData. ``inlet`` and ``shaft`` are axial centre datums;
    there is no discrete outlet port on a radially discharging impeller.
    """

    impeller: object
    results: object = None
    n_theta: int = 160

    def build(self):
        from pyskyfire.pump.impeller_new import cylindrical_to_cartesian
        if not isinstance(self.n_theta, int) or self.n_theta < 3:
            raise ValueError("n_theta must be an integer of at least 3")
        g = self.impeller.geometry
        root = SceneNode("impeller", metadata={"kind": "impeller"})
        for name, curve in (("hub", g.hub_curve), ("shroud", g.shroud_curve)):
            theta = np.linspace(0, 2*np.pi, self.n_theta)
            x, angle = np.meshgrid(np.asarray(curve)[:, 0], theta, indexing="ij")
            radius = np.asarray(curve)[:, 1, None]
            xyz = np.stack((x, radius*np.cos(angle), radius*np.sin(angle)), axis=-1)
            root.children.append(SceneNode(name, mesh=_surface(xyz), opacity=0.22))
        for i, blade in enumerate(g.blade_surfaces):
            root.children.append(SceneNode(f"blade_{i}", mesh=_surface(cylindrical_to_cartesian(blade)), color="tomato", opacity=0.75))
        if self.results is not None:
            known = {node.name: node for node in root.children}
            for name, fields in self.results.items():
                if name not in known:
                    raise ValueError(f"Unknown impeller result part: {name}")
                known[name].fields.update(fields)
        origin = (float(g.hub_curve[0, 0]), 0, 0)
        root.anchors = {"shaft": Anchor.from_direction(origin, (1, 0, 0)),
                        "inlet": Anchor.from_direction(origin, (-1, 0, 0))}
        return root


@dataclass
class PipeComponent:
    """Display an explicitly defined pipe; never route or modify a solver model.

    ``pipe.centerline`` is an (N,3) Cartesian polyline in metres. ``pipe.radius``
    is a positive scalar or one radius per station, in metres. The local frame
    is the frame of the supplied coordinates (no implicit recentering).
    ``results`` maps field names to point-associated FieldData with one value
    per centerline station. Their values are repeated around each ring.
    """

    pipe: object
    results: object = None
    n_theta: int = 16

    def build(self):
        p = np.asarray(self.pipe.centerline, dtype=float)
        if p.ndim != 2 or p.shape[1] != 3 or len(p) < 2 or not np.isfinite(p).all():
            raise ValueError("Pipe centerline must be finite with shape (N,3), N >= 2")
        if not isinstance(self.n_theta, int) or self.n_theta < 3:
            raise ValueError("n_theta must be an integer of at least 3")
        radius = np.broadcast_to(np.asarray(self.pipe.radius, dtype=float), (len(p),))
        if not np.isfinite(radius).all() or (radius <= 0).any():
            raise ValueError("Pipe radii must be finite and positive")
        segments = np.diff(p, axis=0)
        lengths = np.linalg.norm(segments, axis=1)
        if (lengths <= 0).any():
            raise ValueError("Pipe centerline must not contain consecutive duplicate points")
        direction = segments / lengths[:, None]
        tangent = np.vstack((direction[0], direction[:-1]+direction[1:], direction[-1]))
        norms = np.linalg.norm(tangent, axis=1)
        if (norms < 1e-10).any():
            raise ValueError("Pipe centerline must not reverse at a station")
        tangent /= norms[:, None]
        theta = np.linspace(0, 2*np.pi, self.n_theta, endpoint=False)
        rings = []
        u = None
        for center, axis, r in zip(p, tangent, radius, strict=True):
            if u is None:
                u = np.cross(axis, np.eye(3)[np.argmin(np.abs(axis))])
            else:
                u = u - np.dot(u, axis)*axis
                if np.linalg.norm(u) < 1e-10:
                    u = np.cross(axis, np.eye(3)[np.argmin(np.abs(axis))])
            u /= np.linalg.norm(u)
            v = np.cross(axis, u)
            rings.append(center+r*(np.cos(theta)[:, None]*u+np.sin(theta)[:, None]*v))
        fields = {}
        if self.results is not None:
            for name, f in self.results.items():
                if f.association != "point" or len(f.values) != len(p):
                    raise ValueError("Pipe fields must have one point value per centerline station")
                fields[name] = FieldData(np.repeat(f.values, self.n_theta), f.units, f.label)
        return SceneNode("pipe", mesh=_surface(rings, closed=True), fields=fields,
                         anchors=_end_anchors(p), metadata={"kind": "pipe"})
