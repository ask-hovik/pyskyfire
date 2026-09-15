"""Engineering objects adapted into local-coordinate scene hierarchies."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from pyskyfire.ducts import DisplayPath, PipeResult
from pyskyfire.regen.channel_construction import build_channel_construction
from pyskyfire.regen.cad_geometry import transverse_construction
from pyskyfire.regen.channel_layout import ChannelLayout, resolve_analytical_layout, resolve_numerical_layout
from .scene import Anchor, FieldData, MeshData, SceneNode, Transform

__all__ = ["ThrustChamberComponent", "ImpellerComponent", "PipeComponent", "pipe_fields",
           "VoluteComponent", "volute_fields"]

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
        for family, blade, phase, color in (("blade", g.main_blade, 0., "tomato"),
                                             ("splitter", g.splitter, g.splitter_pitch_fraction, "goldenrod")):
            if blade is None:
                continue
            mesh = _impeller_blade_mesh(blade)
            for i in range(g.blade_count):
                root.children.append(SceneNode(
                    f"{family}_{i}", mesh=mesh, color=color,
                    transform=Transform(rotation=((i+phase)*2*np.pi/g.blade_count, 0, 0)),
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


def _impeller_blade_mesh(blade):
    """Join closed-edge sides with shared vertices and cap the hub and tip.

    Caps bridge corresponding pressure/suction stations, avoiding a fan across
    the concave planform of a wrapped blade. These are blade solids for viewing;
    endwalls remain flow-path surfaces, not a manufacturing model of the wheel.
    """
    from pyskyfire.pump.geometry import cylindrical_to_cartesian
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


# PipeResult attribute, field name, base SI units, label. Station arrays carry
# one value per PathStation; connection arrays one per connection. Names match
# the thrust chamber's where the quantity is the same, so one ``color_by`` can
# span a duct and a chamber. No unit conversion is performed.
_PIPE_STATION_FIELDS = {
    "p": ("p_static", "Pa", "Static pressure"),
    "T": ("T_static", "K", "Static temperature"),
    "density": ("density", "kg/m³", "Density"),
    "velocity": ("velocity", "m/s", "Velocity"),
}
_PIPE_CONNECTION_FIELDS = {
    "pressure_loss": ("pressure_loss", "Pa", "Connection pressure loss"),
    "reynolds": ("reynolds", "1", "Reynolds number"),
    "darcy_factor": ("darcy_factor", "1", "Darcy friction factor"),
    "loss_coefficient": ("loss_coefficient", "1", "Loss coefficient"),
}


def pipe_fields(result):
    """Adapt a :class:`~pyskyfire.ducts.PipeResult` into scene fields.

    Station quantities become point fields and connection quantities become
    cell fields, both in the base SI units the solver already works in. NaN is
    retained where a connection does not define a quantity, such as the Darcy
    factor of a lumped bend. :class:`PipeComponent` maps either association
    onto its swept mesh.
    """
    if not isinstance(result, PipeResult):
        raise TypeError("pipe_fields expects a ducts.PipeResult")
    fields = {name: FieldData(getattr(result, attribute), units, label)
              for attribute, (name, units, label) in _PIPE_STATION_FIELDS.items()}
    fields.update({name: FieldData(getattr(result, attribute), units, label, "cell")
                   for attribute, (name, units, label) in _PIPE_CONNECTION_FIELDS.items()})
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


@dataclass
class PipeComponent:
    """Display an explicitly defined pipe; never route or modify a solver model.

    The pipe supplies its own geometry. A :class:`~pyskyfire.ducts.Pipe` offers
    ``display_path(segments=...)``, which opens its lumped rounded bends into
    arcs for drawing without touching the solved model; ``segments`` sets that
    refinement. Any other object is read through ``centerline``, an (N,3)
    Cartesian polyline in metres, and ``radius``, a positive scalar or one
    radius per station, and is drawn exactly as supplied. The local frame is
    the frame of those coordinates (no implicit recentering).

    ``results`` is a :class:`~pyskyfire.ducts.PipeResult`, or a mapping of field
    names to FieldData: point-associated with one value per station, or
    cell-associated with one value per connection. Station values are
    interpolated along refined bends and repeated around each ring; connection
    values fill the rings spanning that connection.
    """

    pipe: object
    results: object = None
    n_theta: int = 16
    segments: int = 16

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
        results = pipe_fields(self.results) if isinstance(self.results, PipeResult) else self.results
        fields = {}
        if results is not None:
            for name, f in results.items():
                fields[name] = _mapped_pipe_field(f, path, self.n_theta)
        return SceneNode("pipe", mesh=_surface(rings, closed=True), fields=fields,
                         anchors=_end_anchors(p), metadata={"kind": "pipe"})


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

        root = SceneNode("volute", anchors=_volute_anchors(volute),
                         metadata={"kind": "volute"})
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
