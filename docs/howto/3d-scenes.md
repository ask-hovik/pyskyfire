# Assemble a 3D scene

`Scene` combines component geometry and simulation fields in a common coordinate
system. Use metres, radians, pascals, kelvin, and other base SI units. There are
no unit-conversion arguments. Scene geometry illustrates the simulated model;
CAD export remains a separate workflow.

```python
import pyskyfire as psf

scene = psf.viz.Scene(color_by="p_static")
chamber_view = scene.add(
    psf.viz.ThrustChamberComponent(
        thrust_chamber,
        cooling_data=cooling_data,
        layout="analytical",
        stride=3,
    ),
    name="chamber",
)
pump_view = scene.add(
    psf.viz.ImpellerComponent(impeller),
    name="fuel_pump",
    position=(-0.45, 0.30, 0.0),
    rotation=(0.0, 1.5707963267948966, 0.0),
)

scene.show()
scene.save_html("engine.html")
report.add_tab("Engine").add_iframe(scene.data_url)
```

The complete assembly example is `examples/visualisation/scene.py`. It reads
saved results without rerunning a simulation.

## Geometry and output

`Scene.add()` builds a component and returns its root `SceneNode`. Nodes hold
local transforms, children, anchors, optional meshes, and scalar fields. Their
path IDs distinguish circuits, channels, hubs, shrouds, and blades, providing
identity for future selection controls. Sibling names must be unique. A scene
can contain several instances of the same component, including overlapping
components at the origin. No fluid-network connectivity is checked.

The scene uses renderer-independent NumPy meshes. Desktop output uses PyVista;
standalone HTML uses bundled Three.js and OrbitControls with no network assets
or Python server. HTML carries geometry, hierarchy, transforms, and original
field values directly, rather than converting a VTK rendering window to glTF.
The HTML viewer currently supports orbit, pan, and zoom, with a fixed color
field chosen in Python. It has no picking or field-selection UI yet.

`show()` owns and closes its desktop window. For custom desktop integration,
`scene.to_plotter()` returns a PyVista plotter whose resources the caller must
close. `Scene` itself owns no rendering window and requires no `close()` call.
Each output reflects current scene properties; an already-open viewer is a
snapshot and is not synchronized with subsequent Python changes.

Cooling-channel construction runs once per resolved seed. Repeated channels
share its immutable mesh and field arrays and use rotation transforms. The
desktop renderer materializes rotated copies; the HTML payload references
shared seed geometry. Each physical channel retains its own node identity.
The numerical layout remains the default. Analytical layout retains its
existing restrictions; this change does not alter the layout solver, stride
behavior, or mesh resolution defaults. Hot-wall construction is not included.

## Cooling fields

`cooling_data` accepts one `RegenResult`, an iterable, or a mapping whose values
are results. Mapping keys have no geometric meaning. Results match by
`circuit_index`, with `circuit_name` checked as an additional safeguard.
Duplicate indices, invalid indices, and mismatched names raise errors. These
checks cannot prove that a result came from the identical geometry revision;
the caller must supply results belonging to the supplied chamber.

`Scene(color_by=None)` retains component colors while keeping attached data.
Set `scene.color_by = "p_static"` before display or export to select a field.
Supported cooling fields are:

| Names | Native result grid | SI units |
| --- | --- | --- |
| `p_static`, `p_stagnation` | `x_coolant` | Pa |
| `T_static`, `T_stagnation`, `coolant_T_sat` | `x_coolant` | K |
| `velocity` | `x_coolant` | m/s |
| `coolant_enthalpy` | `x_coolant` | J/kg |
| `coolant_quality` | `x_coolant` | 1 |
| `T_hot_wall`, `T_coolant_wall` | `x_wall` | K |
| `qpp_hot`, `dQ_dA` | `x_heat_flux` | W/m² |

The wall fields select the hot-wall and coolant-wall columns of `RegenResult.T`;
they color the channel surface as an illustration of those axial results.
Each quantity is linearly interpolated from its own grid onto channel sections.
Descending grids are reversed together with their values. Legacy results without
a native grid fall back to `result.x`. Values outside the result domain and NaN
gaps remain unavailable and appear grey. Categorical fields are not supported.

One common viridis scale covers all visible components providing the selected
field. Components without it retain their material color. A field missing from
every visible component raises an error. Set `scene.clim = (1e6, 5e6)` to fix the
range in SI values. Changing the selected field does not reconstruct geometry.

## Placement and anchors

All component geometry keeps its documented local coordinates. The common frame
is right-handed Cartesian `(x, y, z)`; the existing chamber and impeller axes
are local x. Neither component is automatically recentered. Positive rotations
follow the right-hand rule.

`position=(x,y,z)` and `rotation=(rx,ry,rz)` use translation in metres and
extrinsic xyz Euler angles in radians. A `scipy.spatial.transform.Rotation` is
also accepted. Alternatively supply `transform=psf.viz.Transform(...)`, or
construct one with `Transform(matrix=rigid_4x4_matrix)`. Scaling and reflection
are rejected. Explicit transform and position/rotation cannot be combined.

Use `parent=existing_node` to place a component relative to another node. World
transforms compose parent-to-child. Node visibility also respects ancestors.

```python
world_inlet = chamber_view.anchor("circuit_0.inlet")
pipe_view = scene.add(
    psf.viz.PipeComponent(pipe),
    name="fuel_feed",
    anchor="outlet",
    attach_to=world_inlet,
)
```

This is a one-time visual placement: the two anchor positions coincide and
their outward x directions oppose. The remaining frame axes determine roll.
`node.anchor(name)` returns a world-space snapshot. Attachments are not persistent
constraints: moving a target afterwards does not move the attached component
unless it is also its parent. Both `anchor` and `attach_to` must be supplied;
combining attachment with explicit placement raises an overconstraint error.

Thrust-chamber root anchors are `circuit_N.inlet` and `circuit_N.outlet`. They
refer to the first seed's first channel for that circuit, at its end-section
polygon centre, following coolant flow direction. They are representative
channel endpoints, **not manifold ports**. Individual channel nodes expose
`inlet` and `outlet` anchors too. Impellers expose `shaft` (+x) and `inlet` (-x)
at the axis and first hub-curve station. There is no discrete radial outlet
anchor because the current impeller model has no volute outlet geometry.

## Explicit pipe geometry and extension components

There is not yet a physical pipe class in pyskyfire. `PipeComponent` consumes an
object exposing `centerline` (an N×3 Cartesian array in metres) and `radius`
(a positive scalar or one radius per station). The frame is that of the supplied
coordinates. A swept open-ended polygonal tube follows those stations. Curved
bends must be sampled by the caller; no route or bend radius is inferred.

```python
from types import SimpleNamespace
import numpy as np

pipe = SimpleNamespace(
    centerline=np.array([[0., 0., 0.], [0.1, 0., 0.], [0.2, 0.03, 0.]]),
    radius=0.005,
)
```

Optional `results` maps names to `FieldData` with one value per pipe station:
`FieldData(pressure_values, "Pa", "Static pressure")`. Values repeat around the
section. These values must come from the caller's simulation; visualization
does not compute pressure loss. Translating or rotating the component changes
only its scene placement, never its length, bends, or engineering model.

Custom components implement `build() -> SceneNode`, returning a fresh hierarchy.
`MeshData(points, triangles)` stores local vertices and triangle connectivity.
`FieldData(values, units, label, association="point")` attaches immutable scalar
arrays; `association="cell"` means one value per triangle. Units are descriptive
SI metadata, never conversion inputs. Same-named fields must describe the same
physical quantity across components. Impeller `results` accepts mappings from
part names (`hub`, `shroud`, `blade_N`) to dictionaries of mesh-aligned fields.

## Migration

`make_engine_3d()` and `Engine3DViewer` emit `DeprecationWarning`. The old factory
delegates to the scene implementation and retains its viewer handle for existing
callers. Its `show` parameter now actually controls opening the desktop window.
Replace factory calls with `Scene()` plus `scene.add(ThrustChamberComponent(...))`.
Move `stride`, `n_arc_pts`, `layout`, and `circuit_colors` onto the component;
keep `save_html()` and `data_url` on the scene and remove `close()` calls on it.
