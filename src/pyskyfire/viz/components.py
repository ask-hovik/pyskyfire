"""Engineering objects adapted into local-coordinate scene hierarchies."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
import re

import CoolProp.CoolProp as CP
import numpy as np

from pyskyfire.ducts import AnnularDuct, DisplayPath, DuctResult
from pyskyfire.regen.channel_construction import build_channel_construction
from pyskyfire.regen.cad_geometry import transverse_construction
from pyskyfire.regen.channel_layout import ChannelLayout, resolve_analytical_layout, resolve_numerical_layout
from pyskyfire.thermo.eos import fluid_name
from .scene import DEFAULT_COLOR, Anchor, FieldData, MeshData, SceneNode, Transform

__all__ = ["ThrustChamberComponent", "ImpellerComponent", "DiffuserComponent", "TurbineComponent",
           "PipeComponent", "duct_fields", "AnnularDuctComponent",
           "VoluteComponent", "volute_fields", "ManifoldComponent"]

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


def _fluid_label(name):
    """One name per fluid however CoolProp was asked for it, so groups match.

    Drops the backend prefix and a pure fluid's unit fraction, as
    ``thermo.Fluid`` writes them, then takes CoolProp's own name for an alias
    such as ``H2``. A mixture is spelled by mass, as propellants are mixed,
    rather than by the mole fractions CoolProp takes, so
    ``Ethanol[0.610]&Water[0.390]`` reads ``Ethanol 80 wt% & Water 20 wt%``.
    Fluids CoolProp does not know keep their spelling.
    """
    label = name.split("::", 1)[-1]
    if "&" in label:
        return _mass_label(label)
    label = re.sub(r"\[\s*1(\.0*)?\s*\]$", "", label)
    try:
        label = CP.get_fluid_param_string(label, "name")
    except ValueError:
        pass
    return label


def _mass_label(mixture):
    """A CoolProp ``A[x]&B[y]`` mole spelling as mass percentages."""
    parts = [re.fullmatch(r"\s*([^\[\]]+?)\s*\[\s*([^\[\]]+?)\s*\]\s*", part) for part in mixture.split("&")]
    if not all(parts):
        return mixture
    try:
        names = [CP.get_fluid_param_string(part[1], "name") for part in parts]
        masses = [float(part[2]) * CP.PropsSI("M", name) for part, name in zip(parts, names)]
    except ValueError:
        return mixture
    total = sum(masses)
    return " & ".join(f"{name} {100 * mass / total:.3g} wt%" for name, mass in zip(names, masses))


def _carrying(metadata, fluid):
    """Metadata noting the fluid a result was solved for, when it names one."""
    if fluid is not None:
        metadata["fluid"] = _fluid_label(fluid)
    return metadata


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

    A circuit with a result records its coolant's CoolProp name as
    ``metadata["fluid"]``. ``circuit_colors`` paints each circuit apart;
    without it every channel takes the default colour.
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
        root = SceneNode("chamber", metadata={"kind": "thrust_chamber"})
        groups = {}
        for seed in layout.channels:
            ci = seed.circuit_index
            circuit = chamber.cooling_circuits[ci]
            if ci not in groups:
                fluid = fluid_name(getattr(circuit, "coolant_transport", None)) if ci in results else None
                groups[ci] = SceneNode(f"circuit_{ci}", metadata=_carrying(
                    {"kind": "cooling_circuit", "label": circuit.name, "circuit_index": ci}, fluid))
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
                    anchors=dict(anchors),
                    color=DEFAULT_COLOR if self.circuit_colors is None else self.circuit_colors[ci],
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
        if not isinstance(self.n_theta, int) or self.n_theta < 3:
            raise ValueError("n_theta must be an integer of at least 3")
        g = self.impeller.geometry
        if g is None:
            raise ValueError("Call impeller.build_geometry() before adding it to a scene")
        root = SceneNode("impeller", metadata={"kind": "impeller"})
        for name, curve in (("hub", g.hub_curve), ("shroud", g.shroud_curve)):
            theta = np.linspace(0, 2*np.pi, self.n_theta)
            x, angle = np.meshgrid(np.asarray(curve)[:, 0], theta, indexing="ij")
            radius = np.asarray(curve)[:, 1, None]
            xyz = np.stack((x, radius*np.cos(angle), radius*np.sin(angle)), axis=-1)
            root.children.append(SceneNode(name, mesh=_surface(xyz), opacity=0.22))
        for family, blade, phase in (("blade", g.main_blade, 0.),
                                     ("splitter", g.splitter, g.splitter_pitch_fraction)):
            if blade is None:
                continue
            mesh = _impeller_blade_mesh(blade)
            # The blade is already built in the wheel's own sense, so stepping
            # round to the next passage has to step that way too. It only shows
            # on the splitter, which sits part of a pitch from the main blade;
            # the main blades land on the same ring whichever way they step.
            for i in range(g.blade_count):
                root.children.append(SceneNode(
                    f"{family}_{i}", mesh=mesh,
                    transform=Transform(rotation=(
                        g.rotation_sign*(i+phase)*2*np.pi/g.blade_count, 0, 0)),
                    metadata={"kind": family, "instance_index": i},
                ))
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
class DiffuserComponent:
    """Display a built vaned diffuser: its two walls and every vane.

    ``diffuser.build_geometry()`` must have been called. The diffuser is drawn
    in the impeller's own frame, walls symmetric about the impeller outlet
    mid-plane at local x = 0, so it takes the impeller's rotation and a
    translation to that mid-plane. The ``outlet`` anchor is not a port: the
    diffuser discharges radially all round into its volute.
    """

    diffuser: object
    n_theta: int = 160

    def build(self):
        if not isinstance(self.n_theta, int) or self.n_theta < 3:
            raise ValueError("n_theta must be an integer of at least 3")
        g = self.diffuser.geometry
        if g is None:
            raise ValueError("Call diffuser.build_geometry() before adding it to a scene")
        root = SceneNode("diffuser", metadata={"kind": "diffuser"})
        for name, curve in (("hub", g.hub), ("shroud", g.shroud)):
            theta = np.linspace(0, 2*np.pi, self.n_theta)
            # The vane ends lie in the wall planes. Drawing the walls a hair
            # outside them keeps the renderer from flickering between two
            # coincident surfaces; the design geometry is untouched.
            x, angle = np.meshgrid(1.01*np.asarray(curve)[:, 0], theta, indexing="ij")
            radius = np.asarray(curve)[:, 1, None]
            xyz = np.stack((x, radius*np.cos(angle), radius*np.sin(angle)), axis=-1)
            root.children.append(SceneNode(name, mesh=_surface(xyz), opacity=0.22))
        mesh = _impeller_blade_mesh(g.vane)
        for i in range(g.vane_count):
            root.children.append(SceneNode(
                f"vane_{i}", mesh=mesh,
                transform=Transform(rotation=(g.rotation_sign*i*2*np.pi/g.vane_count, 0, 0)),
                metadata={"kind": "diffuser_vane", "instance_index": i},
            ))
        root.anchors = {"shaft": Anchor.from_direction((0, 0, 0), (1, 0, 0))}
        return root


@dataclass
class TurbineComponent:
    """Display a built turbine: its passage walls and every row's blades.

    ``geometry`` is a :class:`~pyskyfire.turbine.TurbineGeometry`. Walls are
    the passage revolved about the local +x shaft; each row's one built blade
    is repeated round the wheel in the turbine's rotation sense. ``colors``
    paints (stator, rotor) blades apart; without it every part takes the
    default colour. Optional ``results`` maps part names
    (hub, shroud, ``<row>_<i>``) to dictionaries of mesh-aligned FieldData.
    ``inlet`` and ``outlet`` anchors sit on the axis at the passage ends,
    pointing out of the turbine.
    """

    geometry: object
    results: object = None
    n_theta: int = 160
    colors: tuple | None = None

    def build(self):
        if not isinstance(self.n_theta, int) or self.n_theta < 3:
            raise ValueError("n_theta must be an integer of at least 3")
        g = self.geometry
        root = SceneNode("turbine", metadata={"kind": "turbine"})
        for name in ("hub", "shroud"):
            root.children.append(SceneNode(name, mesh=_surface(g.passage.surface(name, self.n_theta)),
                                           opacity=0.22))
        stator_color, rotor_color = (DEFAULT_COLOR,)*2 if self.colors is None else self.colors
        for rg in g.rows:
            mesh = _impeller_blade_mesh(rg.blade)
            count = rg.row.blade_count
            for i in range(count):
                root.children.append(SceneNode(
                    f"{rg.row.name}_{i}", mesh=mesh,
                    color=rotor_color if rg.row.is_rotor else stator_color,
                    transform=Transform(rotation=(g.rotation_sign*i*2*np.pi/count, 0, 0)),
                    metadata={"kind": rg.row.kind, "row": rg.row.name, "instance_index": i},
                ))
        if self.results is not None:
            known = {node.name: node for node in root.children}
            for name, fields in self.results.items():
                if name not in known:
                    raise ValueError(f"Unknown turbine result part: {name}")
                known[name].fields.update(fields)
        inlet = (float(g.passage.hub[0, 0]), 0, 0)
        outlet = (float(g.passage.hub[-1, 0]), 0, 0)
        root.anchors = {"inlet": Anchor.from_direction(inlet, (-1, 0, 0)),
                        "outlet": Anchor.from_direction(outlet, (1, 0, 0)),
                        "shaft": Anchor.from_direction(inlet, (1, 0, 0))}
        return root


def _impeller_blade_mesh(blade):
    """Join closed-edge sides with shared vertices and cap the hub and tip.

    Caps bridge corresponding pressure/suction stations, avoiding a fan across
    the concave planform of a wrapped blade. These are blade solids for viewing;
    endwalls remain flow-path surfaces, not a manufacturing model of the wheel.
    """
    from pyskyfire.common.blade_surface import cylindrical_to_cartesian
    ss = cylindrical_to_cartesian(blade.suction)
    ps = cylindrical_to_cartesian(blade.pressure)
    ring = np.concatenate((ss, ps[:, -2:0:-1]), axis=1)
    surface = _surface(ring, closed=True)
    faces = surface.triangles.tolist()
    n = ss.shape[1]
    nr = ring.shape[1]

    def pressure_index(j):
        return j if j in (0, n-1) else nr-j

    for span, reverse in ((0, False), (len(ring)-1, True)):
        offset = span*nr
        for j in range(n-1):
            a, b, c, d = j, j+1, pressure_index(j+1), pressure_index(j)
            for tri in ((b, a, d), (b, d, c)):
                if len(set(tri)) == 3:
                    faces.append([offset+k for k in (tri[::-1] if reverse else tri)])
    faces = np.asarray(faces, dtype=np.int64)
    p = surface.points
    volume6 = np.einsum("ij,ij->i", p[faces[:,0]], np.cross(p[faces[:,1]],p[faces[:,2]])).sum()
    if volume6 < 0:
        faces = faces[:, ::-1]
    return MeshData(p, faces)


# DuctResult attribute, field name, base SI units, label. Station arrays carry
# one value per PathStation; connection arrays one per connection. Names match
# the thrust chamber's where the quantity is the same, so one ``color_by`` can
# span a duct and a chamber. No unit conversion is performed.
_DUCT_STATION_FIELDS = {
    "p_stagnation": ("p_stagnation", "Pa", "Stagnation pressure"),
    "p_static": ("p_static", "Pa", "Static pressure"),
    "T_stagnation": ("T_stagnation", "K", "Stagnation temperature"),
    "T_static": ("T_static", "K", "Static temperature"),
    "density": ("density", "kg/m³", "Density"),
    "velocity": ("velocity", "m/s", "Velocity"),
    "meridional_velocity": ("meridional_velocity", "m/s", "Through-flow velocity"),
    "swirl_velocity": ("swirl_velocity", "m/s", "Swirl velocity"),
}
_DUCT_INTERVAL_FIELDS = {
    "pressure_loss": ("pressure_loss", "Pa", "Connection pressure loss"),
    "reynolds": ("reynolds", "1", "Reynolds number"),
    "darcy_factor": ("darcy_factor", "1", "Darcy friction factor"),
    "loss_coefficient": ("loss_coefficient", "1", "Loss coefficient"),
}


def duct_fields(result):
    """Adapt a :class:`~pyskyfire.ducts.DuctResult` into scene fields.

    Station quantities become point fields and interval (connection)
    quantities become cell fields, both in the base SI units the solver
    already works in. NaN is retained where an interval does not define a
    quantity, such as the Darcy factor of a lumped bend. :class:`PipeComponent`
    and :class:`AnnularDuctComponent` map either association onto their
    meshes.
    """
    if not isinstance(result, DuctResult):
        raise TypeError("duct_fields expects a ducts.DuctResult")
    fields = {name: FieldData(getattr(result, attribute), units, label)
              for attribute, (name, units, label) in _DUCT_STATION_FIELDS.items()}
    fields.update({name: FieldData(getattr(result, attribute), units, label, "cell")
                   for attribute, (name, units, label) in _DUCT_INTERVAL_FIELDS.items()})
    return fields


def _display_path(pipe, segments):
    """Refined drawing geometry, or the supplied polyline for a plain object."""
    resolve = getattr(pipe, "display_path", None)
    if callable(resolve):
        return resolve(segments=segments)
    points = np.asarray(pipe.centerline, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 2 or not np.isfinite(points).all():
        raise ValueError("Pipe centerline must be finite with shape (N,3), N >= 2")
    radius = np.broadcast_to(np.asarray(pipe.radius, dtype=float), (len(points),))
    return DisplayPath(points, radius, np.arange(len(points), dtype=float),
                       np.arange(len(points) - 1))


def _mapped_pipe_field(f, path, n_theta):
    """Carry one source field onto the refined mesh, keeping its association."""
    if f.association == "point":
        if len(f.values) != path.station_count:
            raise ValueError("Pipe point fields must have one value per centerline station")
        values = np.interp(path.station, np.arange(path.station_count), f.values)
        return FieldData(np.repeat(values, n_theta), f.units, f.label)
    if len(f.values) != path.connection_count:
        raise ValueError("Pipe cell fields must have one value per pipe connection")
    # _surface emits two triangles per ring pair, so a segment owns 2*n_theta cells.
    return FieldData(np.repeat(f.values[path.connection], 2*n_theta), f.units, f.label, "cell")


def _ring_frames(p):
    """Ring centres and two unit vectors spanning each ring's plane, transported along ``p``."""
    legs = np.diff(p, axis=0)
    lengths = np.linalg.norm(legs, axis=1)
    if (lengths <= 0).any():
        raise ValueError("Pipe centerline must not contain consecutive duplicate points")
    direction = legs / lengths[:, None]
    tangent = np.vstack((direction[0], direction[:-1]+direction[1:], direction[-1]))
    norms = np.linalg.norm(tangent, axis=1)
    if (norms < 1e-10).any():
        raise ValueError("Pipe centerline must not reverse at a station")
    tangent /= norms[:, None]
    frames = []
    u = None
    for axis in tangent:
        if u is None:
            u = np.cross(axis, np.eye(3)[np.argmin(np.abs(axis))])
        else:
            u = u - np.dot(u, axis)*axis
            if np.linalg.norm(u) < 1e-10:
                u = np.cross(axis, np.eye(3)[np.argmin(np.abs(axis))])
        u = u/np.linalg.norm(u)
        frames.append((u, np.cross(axis, u)))
    return frames


def _rings(p, frames, radius, n_theta):
    theta = np.linspace(0, 2*np.pi, n_theta, endpoint=False)
    return [center + r*(np.cos(theta)[:, None]*u + np.sin(theta)[:, None]*v)
            for center, (u, v), r in zip(p, frames, radius, strict=True)]


def _spans(present):
    """Point ranges ``(first, last)`` covering every segment with a core at either end."""
    segment = present[:-1] | present[1:]
    spans, start = [], None
    for index, on in enumerate(segment):
        if on and start is None:
            start = index
        if not on and start is not None:
            spans.append((start, index))
            start = None
    if start is not None:
        spans.append((start, len(segment)))
    return spans


@dataclass
class PipeComponent:
    """Display an explicitly defined pipe; never route or modify a solver model.

    The pipe supplies its own geometry. A :class:`~pyskyfire.ducts.Pipe` offers
    ``display_path(segments=...)``, which opens its rounded and continuous bend
    connections into arcs for drawing without touching the solved model;
    ``segments`` sets the drawn sides per quarter turn. Any other object is read through ``centerline``, an (N,3)
    Cartesian polyline in metres, and ``radius``, a positive scalar or one
    radius per station, and is drawn exactly as supplied. The local frame is
    the frame of those coordinates (no implicit recentering).

    A pipe with a core is drawn as a ``"core"`` surface, one per stretch the
    core runs along, inside a translucent ``"bore"`` of ``opacity``, both
    children of the returned node, which carries the bore's end anchors.

    ``results`` is a :class:`~pyskyfire.ducts.DuctResult`, or a mapping of
    field names to FieldData: point-associated with one value per station, or
    cell-associated with one value per connection. Station values are
    interpolated along refined bends and repeated around each ring; connection
    values fill the rings spanning that connection. Every surface carries the
    same fields.
    """

    pipe: object
    results: object = None
    n_theta: int = 16
    segments: int = 16
    opacity: float = 0.35

    def build(self):
        if not isinstance(self.n_theta, int) or self.n_theta < 3:
            raise ValueError("n_theta must be an integer of at least 3")
        path = _display_path(self.pipe, self.segments)
        p = path.points
        if len(p) < 2:
            raise ValueError("Pipe centerline must be finite with shape (N,3), N >= 2")
        radius = path.radius
        if not np.isfinite(radius).all() or (radius <= 0).any():
            raise ValueError("Pipe radii must be finite and positive")
        frames = _ring_frames(p)
        n = self.n_theta
        solved = isinstance(self.results, DuctResult)
        fluid = self.results.fluid if solved else None
        results = duct_fields(self.results) if solved else self.results
        fields = {name: _mapped_pipe_field(f, path, n) for name, f in (results or {}).items()}
        bore = _surface(_rings(p, frames, radius, n), closed=True)
        anchors = _end_anchors(p)
        core = path.core_radius
        if not np.any(core > 0.0):
            return SceneNode("pipe", mesh=bore, fields=fields, anchors=anchors,
                             metadata=_carrying({"kind": "pipe"}, fluid))
        root = SceneNode("pipe", anchors=anchors, metadata=_carrying({"kind": "pipe"}, fluid))
        rings = _rings(p, frames, core, n)
        for number, (first, last) in enumerate(_spans(core > 0.0)):
            # A span's points are a block of n rows per ring, and its
            # segments a block of 2n cells per ring pair.
            span = {name: FieldData(f.values[first*n:(last + 1)*n] if f.association == "point"
                                    else f.values[2*n*first:2*n*last],
                                    f.units, f.label, f.association)
                    for name, f in fields.items()}
            root.children.append(SceneNode("core" if number == 0 else f"core_{number + 1}",
                                           mesh=_surface(rings[first:last + 1], closed=True),
                                           fields=span, metadata={"kind": "pipe_wall"}))
        root.children.append(SceneNode("bore", mesh=bore, fields=dict(fields), opacity=self.opacity,
                                       metadata={"kind": "pipe_wall"}))
        return root


@dataclass
class AnnularDuctComponent:
    """Display an :class:`~pyskyfire.ducts.AnnularDuct`.

    The duct is drawn as its hub and shroud revolved about the local +x
    shaft, with ``inlet`` and ``outlet`` anchors on the axis at its ends as a
    turbine's are. The shroud is translucent so the hub shows. A pipe with a
    core is drawn by :class:`PipeComponent`.

    ``results`` is a :class:`~pyskyfire.ducts.DuctResult`, or a mapping of
    field names to FieldData with one value per station (point) or per
    interval (cell). Both walls carry the same fields.
    """

    duct: object
    results: object = None
    n_theta: int = 64
    opacity: float = 0.35

    def build(self):
        if not isinstance(self.n_theta, int) or self.n_theta < 3:
            raise ValueError("n_theta must be an integer of at least 3")
        if not isinstance(self.duct, AnnularDuct):
            raise TypeError("duct must be a ducts.AnnularDuct")
        solved = isinstance(self.results, DuctResult)
        fluid = self.results.fluid if solved else None
        results = duct_fields(self.results) if solved else self.results
        stations = len(self.duct.s)
        fields = {}
        for name, f in (results or {}).items():
            if f.association == "point":
                if len(f.values) != stations:
                    raise ValueError("Annular duct point fields must have one value per station")
                fields[name] = FieldData(np.repeat(f.values, self.n_theta), f.units, f.label)
            else:
                if len(f.values) != stations - 1:
                    raise ValueError("Annular duct cell fields must have one value per interval")
                # _surface emits two triangles per open ring segment.
                fields[name] = FieldData(np.repeat(f.values, 2*(self.n_theta - 1)), f.units, f.label, "cell")
        root = SceneNode("annular_duct", metadata=_carrying({"kind": "annular_duct"}, fluid))
        for name in ("hub", "shroud"):
            root.children.append(SceneNode(name, mesh=_surface(self.duct.surface(name, self.n_theta)),
                                           fields=dict(fields),
                                           opacity=self.opacity if name == "shroud" else 1.0))
        inlet = (float(self.duct.mean_line[0, 0]), 0, 0)
        outlet = (float(self.duct.mean_line[-1, 0]), 0, 0)
        root.anchors = {"inlet": Anchor.from_direction(inlet, (-1, 0, 0)),
                        "outlet": Anchor.from_direction(outlet, (1, 0, 0))}
        return root


# Per-station quantities a volute passage can be coloured by. Bore radius and
# centroid radius are deliberately absent: they describe the section rather
# than the flow through it, and the design tables already carry them.
_VOLUTE_GEOMETRY_FIELDS = {
    "area": ("area", "m²", "Passage area"),
}


def volute_fields(result, exit_stations=0):
    """Adapt a :class:`~pyskyfire.ducts.VoluteResult` into scene fields.

    Every quantity carries one value per azimuthal station, so all are
    point-associated and :class:`VoluteComponent` repeats each around its
    station's section ring. Units are the base SI units the solver works in.

    ``collected`` marks the fraction of the mass flow gathered so far, and
    ``clipped`` flags the stations near the tongue where the area law asked for
    a section smaller than the family admits.

    ``exit_stations`` appends that many further values covering the exit-cone
    duct, so a component sweeping the duct can colour it too. Pressure and
    velocity carry their real variation along it; the rest hold their
    collection-plane value, since the duct has no azimuth left to vary with.

    The two pressures are worth reading with their model in mind. Aungier
    charges one lumped coefficient per mechanism across the whole passage, so
    how the loss is spread round the volute is an assumption made for display
    -- see :meth:`~pyskyfire.ducts.VoluteResult.stagnation_pressure`. The
    endpoints are what the analysis actually fixes.
    """
    from pyskyfire.ducts import VoluteResult

    if not isinstance(result, VoluteResult):
        raise TypeError("volute_fields expects a ducts.VoluteResult")
    geometry = result.geometry
    stagnation, static = result.stagnation_pressure(), result.static_pressure()
    velocity = result.collected_velocity()
    values = {
        "velocity": (velocity, "m/s", "Collected through-flow velocity"),
        "p_stagnation": (stagnation, "Pa", "Stagnation pressure"),
        "p_static": (static, "Pa", "Static pressure"),
        "collected": (geometry.theta / (2*np.pi), "1", "Mass flow collected"),
        "clipped": (geometry.clipped.astype(float), "1", "Area law clipped at the tongue"),
    }
    for attribute, (name, units, label) in _VOLUTE_GEOMETRY_FIELDS.items():
        values[name] = (np.asarray(getattr(geometry, attribute), dtype=float), units, label)

    exit_stations = int(exit_stations)
    if exit_stations > 0:
        # The duct's own stations continue past the collection plane; drop the
        # first, which is the collection plane the passage already covers.
        duct_stagnation, duct_static = result.exit_cone_pressures(exit_stations + 1)
        duct_velocity = np.linspace(0.0, 1.0, exit_stations + 1)
        eased = duct_velocity**2 * (3.0 - 2.0*duct_velocity)
        duct_velocity = (1.0 - eased)*result.C_6 + eased*result.C_7
        along = {
            "velocity": duct_velocity[1:],
            "p_stagnation": duct_stagnation[1:],
            "p_static": duct_static[1:],
        }
        values = {
            name: (np.concatenate((data, along.get(name, np.full(exit_stations, data[-1])))),
                   units, label)
            for name, (data, units, label) in values.items()
        }
    return {name: FieldData(data, units, label)
            for name, (data, units, label) in values.items()}


def _fan(ring, origin, *, outward):
    """Triangulate one end section as a fan about a point it is star-shaped about."""
    count = len(ring)
    points = np.vstack((ring, origin))
    index = np.arange(count)
    spokes = (index, (index + 1) % count)
    if outward:
        spokes = spokes[::-1]
    return MeshData(points, np.column_stack((np.full(count, count), *spokes)))


def _merge(meshes):
    """Concatenate meshes, offsetting each one's triangle indices."""
    points, triangles, offset = [], [], 0
    for mesh in meshes:
        points.append(mesh.points)
        triangles.append(mesh.triangles + offset)
        offset += len(mesh.points)
    return MeshData(np.vstack(points), np.vstack(triangles))


@dataclass
class VoluteComponent:
    """Display a resolved volute or collector passage.

    The volute supplies its own geometry through ``sweep(n_points=...)``, which
    places each station's cross-section in machine coordinates, already
    inclined by the volute's tilt. The local frame is that frame: the machine
    axis is +z, the tongue sits at azimuth zero, and no recentering is applied.

    The passage is swept from the tongue through one revolution, and on through
    its :class:`~pyskyfire.ducts.ExitCone` duct to the flange when it has one.
    The face at the tongue is a wall -- the cutwater -- so the step visible
    there is the passage meeting itself.

    Three child nodes come back, so the openings can be told apart from the
    walls at a glance:

    ``passage``
        Every wall, including the whole exit-cone duct and the cutwater face.
        This is what carries the fields.
    ``inlet``
        The slot mouth, the strip at the inlet radius running the full
        revolution, where the flow being collected enters.
    ``outlet``
        The flange face, or the collection plane itself without an exit cone.

    ``opening_color`` paints both openings, red by default. They carry no
    fields, so they keep that colour whichever field the scene is coloured by,
    rather than disappearing into the scale.

    ``results`` is a :class:`~pyskyfire.ducts.VoluteResult`, or a mapping of
    field names to point-associated FieldData with one value per azimuthal
    station. A supplied mapping holds its last value along the exit cone, which
    a result does only for quantities that have no variation to carry there.

    ``exit_points`` sets how finely the exit cone is swept along its length.

    Anchors mark the inlet, on the axis pointing in along -z, and the outlet:
    the flange face, oriented along its own normal, when there is an exit
    cone; otherwise the full-collection plane, oriented tangentially into
    where an exit cone would continue.
    """

    volute: object
    results: object = None
    n_points: int = 16
    exit_points: int = 17
    opening_color: str = "red"

    def build(self):
        if not isinstance(self.n_points, int) or self.n_points < 2:
            raise ValueError("n_points must be an integer of at least 2")
        volute = self.volute
        # polygon() repeats its first point to close the loop; a swept mesh
        # wraps that seam itself, so drop it rather than emit degenerate faces.
        rings = np.asarray(volute.sweep(n_points=self.n_points), dtype=float)[:, :-1]
        if rings.ndim != 3 or rings.shape[2] != 3 or not np.isfinite(rings).all():
            raise ValueError("Volute sweep must be finite with shape (N,M,3)")
        stations, width = rings.shape[:2]
        origins = _volute_fan_origins(volute)
        # The section's last boundary piece is the slot mouth, so its vertices
        # are the ones tracing the inlet opening around the volute.
        mouth = width - self.n_points

        exit_rings = np.zeros((0, width, 3))
        if volute.exit_geometry is not None:
            # The duct's first ring is the volute's own last one, so it is
            # dropped here and the two sweeps share that seam exactly.
            exit_rings = np.asarray(
                volute.sweep_exit_cone(n_points=self.n_points,
                                       n_stations=self.exit_points),
                dtype=float,
            )[1:]

        walls = [_surface(rings[:, :mouth], closed=False)]
        if len(exit_rings):
            # Every wall of the duct, mouth-side included: nothing enters it.
            walls.append(_surface(np.concatenate((rings[-1:], exit_rings)), closed=True))
            closing_ring, closing_origin = exit_rings[-1], volute.exit_geometry.end
        else:
            closing_ring, closing_origin = rings[-1], origins[-1]
        walls.append(_fan(rings[0], origins[0], outward=True))

        fluid = None if self.results is None or isinstance(self.results, Mapping) else self.results.fluid
        root = SceneNode("volute", anchors=_volute_anchors(volute),
                         metadata=_carrying({"kind": "volute"}, fluid))
        passage = SceneNode("passage", mesh=_merge(walls),
                            fields=self._wall_fields(stations, mouth, len(exit_rings)),
                            metadata={"kind": "volute-passage"})
        inlet = SceneNode("inlet", mesh=_surface(rings[:, mouth:], closed=False),
                          color=self.opening_color,
                          metadata={"kind": "volute-inlet"})
        outlet = SceneNode("outlet", mesh=_fan(closing_ring, closing_origin, outward=False),
                           color=self.opening_color,
                           metadata={"kind": "volute-outlet"})
        root.children = [passage, inlet, outlet]
        return root

    def _wall_fields(self, stations, mouth, exit_stations):
        """Field values laid out to match the merged wall mesh's vertices."""
        results = self.results
        if results is None:
            return {}
        if not isinstance(results, Mapping):
            results = volute_fields(results, exit_stations=exit_stations)
            expected = stations + exit_stations
        else:
            expected = stations
        fields = {}
        for name, field in results.items():
            if field.association != "point" or len(field.values) != expected:
                raise ValueError(
                    "Volute fields must be point-associated, one value per station")
            values = np.asarray(field.values, dtype=float)
            if len(values) == stations and exit_stations:
                values = np.concatenate((values, np.full(exit_stations, values[-1])))
            parts = [np.repeat(values[:stations], mouth)]
            if exit_stations:
                # The duct mesh repeats the shared seam ring, then its own.
                parts.append(np.repeat(values[stations - 1:], self.n_points * 6))
            # The cutwater cap: its ring, then its fan origin.
            parts.append(np.full(self.n_points * 6 + 1, values[0]))
            fields[name] = FieldData(np.concatenate(parts), field.units, field.label)
        return fields


def _volute_fan_origins(volute):
    """Star-shaped fan origins for the end sections, in machine coordinates."""
    geometry = volute.geometry
    local = volute.section.fan_origin(geometry.size)
    theta = geometry.theta
    radial = np.column_stack((np.cos(theta), np.sin(theta), np.zeros_like(theta)))
    axial = np.zeros_like(radial)
    axial[:, 2] = 1.0
    u_hat = np.cos(volute.tilt)*radial + np.sin(volute.tilt)*axial
    v_hat = -np.sin(volute.tilt)*radial + np.cos(volute.tilt)*axial
    return (volute.inlet.radius*radial
            + local[:, :1]*u_hat + local[:, 1:2]*v_hat)


def _volute_anchors(volute):
    """Inlet on the machine axis; outlet at the flange, or the collection plane."""
    inlet = Anchor.from_direction(np.zeros(3), np.array([0.0, 0.0, -1.0]))
    if volute.exit_geometry is not None:
        cone = volute.exit_geometry
        return {"inlet": inlet, "outlet": Anchor.from_direction(cone.end, cone.end_direction)}
    geometry = volute.geometry
    centre = _volute_fan_origins(volute)[-1]
    tangential = np.array([-np.sin(geometry.theta[-1]), np.cos(geometry.theta[-1]), 0.0])
    return {"inlet": inlet, "outlet": Anchor.from_direction(centre, tangential)}


@dataclass
class ManifoldComponent:
    """Display a resolved distributing or collecting manifold.

    The manifold supplies its own geometry through ``sweep(n_points=...)``,
    already placed on the wall it was built against, and the local frame is
    that frame: the machine axis is +z and the dead point sits at azimuth zero.
    A split manifold sweeps the whole ring, both branches together.

    Three child nodes come back, so the opening can be told apart from the
    walls at a glance:

    ``passage``
        Every wall of the cavity, and the cutwater face where a single-branch
        ring meets itself. This is what carries the fields.
    ``passage``
        Every wall: the dome, the stub's own tube, the chamber-facing face
        outside the slit, and the cutwater where a single-branch ring meets
        itself. This is what carries the fields.
    ``slit``
        The strip of the chamber-facing face reserved for the cooling channels
        to break out through. A reservation rather than modelled geometry --
        the breakout itself is CAD work -- so what is drawn is where it has
        been kept clear. A section that reserves nothing gives an ``opening``
        instead, its whole face being a real hole.
    ``outlet``
        The flange the stub ends on, where there is one.

    ``opening_color`` paints only those faces, red by default: the places flow
    crosses the boundary, not the walls that carry it there. They hold no
    fields, so they keep that colour whichever field the scene is coloured by,
    rather than disappearing into the scale.

    ``results`` is a mapping of field names to point-associated FieldData, one
    value per swept station of the resolved ring -- which for a split manifold
    is both branches, not one.

    ``n_points`` samples each boundary piece of the section. ``n_stations``
    refines the ring in azimuth for drawing, which is the direction a manifold
    facets in first, and leaves the resolved design alone; values on the
    resolved grid are read between its stations onto the refined mesh. ``None``
    draws the ring exactly as it was resolved.
    """

    manifold: object
    results: object = None
    n_points: int = 16
    port_points: int = 9
    n_stations: int | None = None
    opening_color: str = "red"

    def build(self):
        if not isinstance(self.n_points, int) or self.n_points < 2:
            raise ValueError("n_points must be an integer of at least 2")
        manifold = self.manifold
        # polygon() repeats its first point to close the loop; a swept mesh
        # wraps that seam itself, so drop it rather than emit degenerate faces.
        rings = np.asarray(manifold.sweep(n_points=self.n_points,
                                          n_stations=self.n_stations),
                           dtype=float)[:, :-1]
        if rings.ndim != 3 or rings.shape[2] != 3 or not np.isfinite(rings).all():
            raise ValueError("Manifold sweep must be finite with shape (N,M,3)")
        stations, width = rings.shape[:2]

        # The section says which of its boundary pieces faces the machine it
        # serves, rather than closing the passage.
        index = manifold.section.opening_index
        if index == 0:
            face, walls = rings[:, :self.n_points], rings[:, self.n_points:]
        elif index == -1:
            face, walls = rings[:, -self.n_points:], rings[:, :-self.n_points]
        else:
            raise ValueError(
                "a section's opening must be its first or last boundary piece "
                f"for a manifold to sweep it apart; this one reports {index}")

        meshes = [_surface(walls, closed=False)]
        # One station index per vertex of each mesh, so fields laid out per
        # station can be spread over a merged wall of mixed widths.
        order = [np.repeat(np.arange(stations), width - self.n_points)]
        openings = []

        if hasattr(manifold.section, "slit_boundary"):
            # A section measured against a wall meets it everywhere; only the
            # slit is open. The rest of that face is metal, so it joins the
            # walls rather than being painted as a hole.
            for strip in manifold.sweep_margins(self.n_points, self.n_stations):
                meshes.append(_surface(np.asarray(strip, dtype=float), closed=False))
                order.append(np.repeat(np.arange(stations), self.n_points))
            openings.append(SceneNode(
                "slit",
                mesh=_surface(np.asarray(manifold.sweep_slit(self.n_points,
                                                            self.n_stations),
                                         dtype=float), closed=False),
                color=self.opening_color, metadata={"kind": "manifold-slit"}))
        else:
            openings.append(SceneNode(
                "opening", mesh=_surface(face, closed=False),
                color=self.opening_color, metadata={"kind": "manifold-opening"}))

        root = SceneNode("manifold", metadata={"kind": "manifold"})
        if manifold.port_geometry is not None:
            duct = manifold.port_geometry
            stub = np.asarray(manifold.sweep_port(n_points=self.n_points,
                                                  n_stations=self.port_points),
                              dtype=float)
            # The stub's wall is wall; only the flange it ends on is an opening.
            meshes.append(_surface(stub, closed=True))
            order.append(np.full(stub.shape[0]*stub.shape[1], stations - 1))
            openings.append(SceneNode(
                "outlet", mesh=_fan(stub[-1], duct.end, outward=True),
                color=self.opening_color, metadata={"kind": "manifold-outlet"}))
            root.anchors = {"port": Anchor.from_direction(duct.end,
                                                          duct.end_direction)}

        if manifold.branches == 1:
            # A ring running one way round meets itself at full size, so the
            # step there is a wall: the cutwater.
            meshes.append(_fan(rings[0], manifold.fan_origins(self.n_stations)[0],
                               outward=True))
            order.append(np.full(width + 1, 0))

        root.children = [
            SceneNode("passage", mesh=_merge(meshes),
                      fields=self._wall_fields(stations, np.concatenate(order)),
                      metadata={"kind": "manifold-passage"}),
            *openings,
        ]
        return root

    def _wall_fields(self, stations, order):
        """Field values laid out to match the merged wall mesh's vertices."""
        results = self.results
        if results is None:
            return {}
        if not isinstance(results, Mapping):
            raise TypeError(
                "manifold results must be a mapping of field names to "
                "FieldData; there is no loss model to adapt yet")
        # Values belong to the ring the law was resolved on; a mesh refined for
        # drawing reads between them, the way a pipe's do along an opened bend.
        resolved = self.manifold.sweep_stations()
        fields = {}
        for name, field in results.items():
            if field.association != "point" or len(field.values) != resolved:
                raise ValueError(
                    "Manifold fields must be point-associated, one value per "
                    "swept station of the resolved ring")
            values = np.asarray(field.values, dtype=float)
            if stations != resolved:
                values = np.interp(np.linspace(0.0, 1.0, stations),
                                   np.linspace(0.0, 1.0, resolved), values)
            fields[name] = replace(field, values=values[order])
        return fields
