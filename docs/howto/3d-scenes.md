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
standalone HTML uses bundled Three.js and TrackballControls with no network assets
or Python server. HTML carries geometry, hierarchy, transforms, and original
field values directly, rather than converting a VTK rendering window to glTF.
The HTML viewer supports unconstrained trackball rotation, pan, zoom, and an
`R` key that restores the initial view.

Clicking a part selects the top-level node that `scene.add()` created, so a
click on one cooling channel identifies its whole thrust chamber. A box appears
beside the click showing the node's name, its `kind` metadata, the `position=`
and `rotation=` arguments that reproduce its placement, the fields it carries,
and its anchors in scene coordinates. The angles use the same extrinsic xyz
convention as `Transform`, so they can be pasted back into `scene.add()`
verbatim. Click empty space or press `Esc` to dismiss the box. Selection is
HTML-only; the desktop PyVista output has no picking.

The legend carries a dropdown that recolors the scene by any field its
components provide, plus a "None" entry that restores component colors. Every
field travels in the HTML already, so switching costs no extra payload and
needs no rerun. `color_by` chooses the field shown first, and `clim` pins that
field's range while the others range over their own values. A field whose units
disagree between components is left out of the dropdown.

`Scene(show_axes=True)` adds a corner orientation triad, labelled x, y, and z
and colored red, green, and blue, to both outputs. It is an orientation cue
that rotates with the view; it is not a measuring grid, and marks neither the
origin nor the scale.

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

## Pipes and ducts

`PipeComponent` draws a `ducts.Pipe` directly; the solved model is never routed
or modified.

```python
pipe = psf.ducts.Pipe(stations, connections, roughness=15e-6)
result = pipe.solve(inlet_state=inlet, mass_flow=0.1, transport=transport)

scene.add(psf.viz.PipeComponent(pipe, result), name="fuel_feed")
```

A lumped `RoundedBend` is one connection between its two tangent points, so
`Pipe.centerline` crosses the corner as a chord. For drawing, `PipeComponent`
calls `Pipe.display_path(segments=...)`, which reopens each such bend into the
arc its resolved `center` and `normal` describe. `segments` on the component
sets that refinement, defaulting to 16. A `RoundedBend` built by hand without
`center` and `normal` keeps its chord, because the bend plane is not
recoverable from the loss model alone. `Continuous` connections are already
physical arc stations, so their resolution stays with their own `segments`.

Refinement is display-only. Lengths, losses, and every other solver quantity
come from the stations and connections, which `display_path()` does not touch.

Passing a `PipeResult` attaches its fields through `psf.viz.pipe_fields()`:
`p_static`, `T_static`, `density`, and `velocity` per station, plus
`pressure_loss`, `reynolds`, `darcy_factor`, and `loss_coefficient` per
connection. Station names match the thrust chamber's, so one `color_by` spans a
duct and a chamber. Per-connection values become cell fields on the rings
spanning that connection; NaN where a connection does not define the quantity,
such as the Darcy factor of a lumped bend.

`results` also accepts a plain mapping of names to `FieldData`, either
point-associated with one value per station or cell-associated with one value
per connection:

```python
psf.viz.PipeComponent(pipe, {"p_static": FieldData(values, "Pa", "Static pressure")})
```

Any other object exposing `centerline` (an N×3 Cartesian array in metres) and
`radius` (a positive scalar or one radius per station) is accepted too, and is
swept exactly as supplied with no refinement:

```python
from types import SimpleNamespace
import numpy as np

pipe = SimpleNamespace(
    centerline=np.array([[0., 0., 0.], [0.1, 0., 0.], [0.2, 0.03, 0.]]),
    radius=0.005,
)
```

Translating or rotating the component changes only its scene placement, never
its length, bends, or engineering model.

## Extension components

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
