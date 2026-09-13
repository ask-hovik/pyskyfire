# Physical cooling-channel construction

`build_channel_construction(circuit, reference_path)` resolves a circuit's
nominal winding into the section it physically occupies. The same resolved
geometry feeds the simulation, CAD, and the native 3D viewer, so a channel is
sized once and drawn, meshed and solved from that one answer.

## The three paths

A circuit carries three cylindrical `(x, r, theta)` paths, and they answer
three different questions.

`reference_path` is the ideal contour carrying the nominal winding. Its radius
is `contour.r(x)` exactly, and its angular rate comes from the placement's
helix law. Nothing about the cross section or the wall stack enters it, which
is the point: two circuits asking for the same helix angle wind together and
cannot drift apart, whatever sections they carry.

`cross_section_base_point_path` is where the section is anchored. For a squared
section that is the midpoint of the coolant-facing floor; for a rounded one it
is the inner lobe centre. Section planes, drawing arcs and CAD sketches are all
built here.

`cooling_channel_centroid_path` is the coolant area centroid, computed from the
resolved section rather than assumed to be half-way up the channel. The
one-dimensional passage model takes its length, `ds/dx` and curvature from this
path, so a channel that sits further out correctly runs further.

`ThrustChamber` builds them in that order: the winding first, then the placed
sections, then the centroid the sections enclose. Swept volume pairs the
section area with the reference path's `ds/dx`, because that is the plane the
area was measured in; using the centroid there would make the volume of a
helically wound annulus depend on its helix angle, which it physically does
not.

Saved geometry from before these definitions is rejected on load rather than
renamed, because the stored coordinates meant something else. Regenerate it.

## Main walls and liners

`Wall` accepts `role="main"` (the default) or `role="liner"`:

```python
tube_wall = psf.regen.Wall(material=copper, thickness=0.6e-3)
liner = psf.regen.Wall(material=silver, thickness=insert_thickness, role="liner")
# Pass walls=[liner, tube_wall] to CoolingCircuit.
```

List layers from hot gas to coolant, with all liners before the main walls.
A circuit requires at least one main wall. Every main layer contributes to
the channel-wall thickness; for rounded channels these are the layers forming
the tubes and the material between neighboring coolant passages. A liner lies
only on the combustion side and does not add to the tube-wall thickness.
All layers still participate in through-wall thermal conduction.

`circuit.channel_wall_thickness(x)` and design input `t_wall` sum the main
layers. `circuit.hot_side_layer_thickness(x)` and `t_hot_layer` sum the liners.
`circuit.total_thickness(x)` and `t_wall_total` include both, describing the
complete stack between hot gas and coolant. The wall-layer plot includes all
layers on the hot side and assumes a closeout made from the main layers only.

Names and list positions do not assign roles. Existing calls and saved `Wall`
objects without a role default to main walls; explicitly mark existing liners
and regenerate saved simulations. The RL10 silver throat insert uses
`role="liner"`, while its stainless steel tube wall uses the main default.

## Rib heat transfer

Set `enable_fin` and `rib_material` on each `CoolingCircuit` independently:

```python
circuit = psf.regen.CoolingCircuit(
    # ... other circuit inputs ...
    enable_fin=True,
    rib_material=psf.common.solids.GRCop42,
)
```

`enable_fin` defaults to `True`. Set it to `False` to omit the fin heat-transfer
enhancement while retaining the rib geometry. `rib_material` defaults to `None`,
which uses the coolant-facing main-wall material. An explicit material supplies
`get_k(T)` to the fin calculation only; wall conduction still uses the individual
wall-layer materials. For rounded tubes an override represents an effective fin
material, not an additional wall layer.

`ThrustChamber` no longer accepts `enable_fin` or overwrites circuit settings.
Move any previous chamber-level setting onto each relevant circuit. Saved
circuits retain their existing enabled/disabled setting and default to the
wall material if they predate the material override.

## Squared channels

Squared sections use planes perpendicular to the engine axis, for both
straight and helical channels. Each origin sits on the base-point path, and the
first in-plane axis points radially outward. Consequently the engine axis is at
`(-r_base, 0)` in local plane coordinates. The channel floor and roof are exact
circles about that point.

`build_axis_normal_section_planes(path)` provides these planes directly.
They are also available through `CoolingCircuit.section_planes(axis_normal=True)`
and `ResolvedChannel.axis_normal_section_planes`. The original flow-normal
`section_planes` interfaces retain their behavior.

Let `alpha = atan(dr_ref/dx)` be the signed meridional slope angle, using the
existing `alpha` convention, and let `gamma` be the helix angle measured from
the meridian. The physical drawing dimensions are

```text
height_transverse = h / cos(alpha)
pitch_transverse  = theta / cos(gamma_nominal)
```

Here `theta` is the simulation's local angular pitch, already reduced by
`cos(gamma_nominal)`. Construction uses the stored nominal helix law. If a
numerical layout changes a seed's angular path, its derivative changes the
local `gamma` used for rib projection; the nominal pitch allocation stays the
same. A geometry-only circuit without a stored helix law uses the derivatives
of its own path. Cylindrical derivatives avoid interpreting a coarsely sampled
helix chord as an inward radial slope.

The complete wall stack also needs the meridional projection. The reference
path is the bare contour, so the whole stack is projected onto it once:

```text
r_base = r_ref + (t_hot_layer + t_wall) / cos(alpha)
r_roof = r_base + h / cos(alpha)
```

## Rib projection

A simple wedge would produce a rib that thickens toward the roof. The
construction instead intersects the floor and roof circles with the two faces
of each rib.

Rib faces run along the local radial direction, so a transverse section shows
ribs that point at the engine axis rather than leaning over. Writing `u` for
radial and `v` for circumferential offsets from a rib's own centre ray, its two
faces are the parallel lines

```text
v = +d
v = -d
```

A circle of radius `R` about the engine axis meets them at

```text
delta = +/- pitch_transverse/2 -/+ asin(d / R)
```

so the floor and roof arcs run between the inward-facing intersections at the
rib centres `-pitch_transverse/2` and `+pitch_transverse/2`, joined by the
straight faces. Invalid intersections and a closed channel opening raise
errors.

`d` is not simply half the rib width. Both faces of a rib are helicoids: their
separation inside a transverse plane is constant, but the material between them
is that separation times the cosine of the helix angle seen locally, and one
angular rate `dtheta/dx` is a steeper helix further out. A rib is therefore
thinnest at the channel roof, and

```text
d = rib_width/2 * sqrt(1 + (r_roof * dtheta/dx)**2)
```

puts that thinnest point at exactly `rib_width`, leaving the floor a little
thicker. The rate belongs to the reference path, where the helix angle is
measured: `dtheta/dx = tan(gamma) / (r_ref cos(alpha))`. With no helix this
reduces to the usual constant-width rib between annular channels.

The rib face normal in `(axial, radial, circumferential)` coordinates is

```text
n_rib = (r dtheta/dx, 0, -1) / sqrt(1 + (r dtheta/dx)**2)
```

at radius `r`. It has no radial component, because the face contains the radial
direction.

These are local tangent-plane projections: `alpha`, `gamma` and `dtheta/dx` are
held constant across a station's section, and the wall is treated as locally
straight. Contour curvature, taper and loft interpolation still introduce
differences between the local model and the resulting solid. They fall on the
safe side: on a tapering wall the neighbouring channel's roof edge has already
moved outward by the time two faces close on each other, so the solid keeps a
little more material than the model asks for, never less. The station spacing controls the axial
loft approximation. This construction avoids a CAD kernel in Python while
retaining chamber-centered section arcs and explicit rib faces.

## Choosing the rounded height input

Rounded sections accept the full internal coolant height by default. To use
its older base-height parameter, set the option on the cross section:

```python
cross_section = psf.regen.CrossSectionRounded(height_input="base")
# The same option is available on CrossSectionRoundedInternal.
```

Then `CoolingCircuit(channel_height=...)` supplies the base parameter. This
is the value used directly as `prof.h` in commit `ca6d36a`, before the switch
to full-height inputs. It excludes the rounded lobes and is distinct from the
straight-side length, which is `h_base*cos(theta/2)`.

Both lobes are tangent to the rays bounding the channel's sector, so raising
the outer lobe centre by `h_base` also grows its radius by
`h_base*sin(theta/2)`:

```text
q_inner = d_lobe*sin(theta/2) - t_wall
H_full  = 2*q_inner + h_base*(1 + sin(theta/2))
```

`d_lobe` is the solved distance from the engine axis to the inner lobe centre
(see below) and `theta` the sector angle in the section plane. The straight
side length is a separate quantity, `h_base*cos(theta/2)`.

Zero base height is valid: it produces the circular limiting section. Negative
or nonfinite heights are rejected. Full-height input must still meet the
geometric minimum `2*q_inner`.

The section exposes `full_height(profiles)` and `base_height(profiles)` in
either input mode. After chamber assembly, use
`circuit.channel_full_heights` and `circuit.channel_base_heights` for the
corresponding arrays on `circuit.x_domain`. The original `channel_height(x)`,
`channel_heights`, and design-input quantity `h` retain the supplied convention.

Height plots always default to full internal height, independently of the
input convention:

```python
psf.viz.PlotChannelHeight(thrust_chamber)                 # full height
psf.viz.PlotChannelHeight(thrust_chamber, height="base") # base height
psf.viz.PlotChannelHeight(thrust_chamber, height="both") # both, labeled

# The circuit selection option works with either representation.
plot = psf.viz.PlotChannelHeight(
    thrust_chamber, circuit_index=0, height="both"
)
report.add_tab("Channel geometry").add_figure(plot)
```

Squared channels have the same full and base height, so plotting both adds
only one trace for those circuits. The wall-layer plot also uses full height.
Native 3D and Onshape geometry use the section's selected interpretation.
The coolant centroid path uses the full height, after the local pitch and the
lobe placement have been resolved. Equivalent base and full inputs therefore
give the same section analytics, coolant path, and construction geometry when
the other resolved design inputs agree. A helix or width law that depends on the input `h` sees
the selected convention and may require a corresponding change to that law.

## Rounded channels and Onshape

`CrossSectionRounded` and `CrossSectionRoundedInternal` build their lobes in
planes normal to the nominal winding, and both lobes are tangent to the two
rays that separate a channel from its neighbours. The lobes are therefore sized
at the radius they end up at: there is no outward translation afterwards, which
is what used to leave the wrong amount of material between neighbouring tubes.

Where the lobes sit follows from the gas area. A locally prismatic channel
meets an engine-axis-normal plane in the image of its section swept along the
flow direction, so in a local radial/circumferential frame that projection is

```text
radial displacement from normal-plane u:          u / cos(alpha)
circumferential displacement from normal-plane u: u tan(alpha) tan(gamma)
circumferential displacement from normal-plane v: v / cos(gamma)
```

Mapping the two sector rays back through that projection gives the sector the
section really owns, with half-angle `theta/2`. The lobe-centre distance
`d_lobe` is then solved so that the gas remaining inside the sector, outside
the tube, equals that sector's share of the ideal contour's area,
`pitch * r_ref**2 / 2`. Every channel keeps the gas area of the contour it
replaces, and the tube exterior is tangent to the sector boundary, so
neighbouring tubes touch. A liner extends the gas-facing arc inward from the
tube and is included in the same solve.

Because the tube bulges inward between its tangent points and leaves a cusp of
material at them, `d_lobe` is not the contour radius plus a wall thickness. On
a strongly sloped wall it can even sit inside the contour, which is the correct
answer for a tube bundle whose area-equivalent contour is what was specified.

At the minimum rounded height, the two arcs form a circle and FeatureScript
omits the zero-length side lines.

The advanced example's exporter writes schema version 6. Both channel types
are exported in engine-axis-normal planes. The solver's rounded sections remain
flow-normal; `regen.cad_geometry.transverse_construction` provides a separate
geometry-only representation for CAD and the native viewer. It intersects the
locally prismatic section with a transverse plane, including both contour slope
and helix angle, the shear when both are nonzero, and the rotated lobe frame.
Circular arcs become elliptical arcs; straight sides stay straight. Zero helix
only gives circular arcs where the contour slope is also zero. A curved CAD loft
is still an approximation between those local sections.

Each channel carries `sectionType: "squared"` or `"rounded"`. Each station
contains a transverse plane, two affine arcs, and a circular-limit flag. An arc
is represented by `center + axisA*cos(t) + axisB*sin(t)`, with directed angle
limits. Its JSON also includes principal-axis ellipse parameters for Onshape's
`skEllipticalArc`. Circles use `skArc`; zero-length sides are omitted at the
rounded circular limit. FeatureScript lofts and patterns coolant seeds
independently. Regenerate older JSON with the updated exporter and update the
FeatureScript together.

The native viewer samples these same transverse arcs into polygons. Its
`n_arc_pts` setting controls circumferential tessellation and `stride` controls
the selected stations; the final station is always retained. Onshape retains
analytic arcs and uses all exported coolant stations. Neither Python tests nor
the native mesh renderer execute Onshape's loft kernel.

## Optional continuous hot-wall export

The advanced and RL10 export entry points share one exporter. To include the
hot wall in the same JSON:

```bash
python examples/advanced/export_onshape.py --hot-wall \
    --hot-wall-trim-fraction 0.1 --hot-wall-transition-fraction 0.02

python validation/RL10/export_onshape.py --hot-wall \
    --hot-wall-transition-length-m 0.01
```

Enable **Construct hot wall** in the included Onshape cooling-channel feature.
The optional `hotWall` block supplies a complete closed circumferential sketch
at each hot-wall station. The surface is open at its axial ends; it is not a
closed solid or a material-volume model. The native viewer currently shows the
transverse coolant geometry, not this additional hot-wall surface.

A hot-wall profile holds one curve per squared channel, or an arc and a blend
per rounded channel, so a finely divided engine reaches hundreds of similar,
evenly spaced curves. Left to pair those profiles by proximity, a loft fails
outright past roughly ninety curves, and once a helix advances more than half a
channel pitch between stations it pairs each curve with its neighbour, twisting
the wall against the channels. The export therefore also carries `runs`: spans
of consecutive sections whose curves are the same channels in the same order,
which `join_hot_wall` guarantees by pinning the profile seam to a stable
channel. FeatureScript lofts each curve along its run on its own, so every
profile is a single edge and correspondence is given rather than inferred, at
any channel count or helix angle.

A circuit hand-off replaces one channel pattern with another and starts a new
run, so the two profiles bounding it share no curves at all. Both still wind
once about the axis, so azimuth pairs them: the export cuts each at every curve
start of either profile, in `hotWall.bridges`, and the two then carry the same
number of curves spanning the same angular intervals in the same order.
Splitting is exact — circular arcs by parameter, Béziers by de Casteljau — and
a profile keeps its own curve starts wherever the two nearly coincide, so
neither side gains a sliver the other's cut would have made. Paired spans sit
at the same azimuth, which makes pairing them in order the nearest pairing too,
so a hand-off is lofted a couple of dozen spans at a time rather than one.

The resulting surface strips are knitted into one surface body, and are left
separate with a notice if the kernel declines to knit them. A hand-off is the
place that is most likely to decline: the run and the bridge meet on the same
curve, but the bridge has split it, so knitting there depends on the kernel
imprinting the extra vertices rather than on edges that already match.

One strip is one loft, so a hot wall covering several hundred channels is
several hundred operations and regenerates slowly. Exporting one circuit at a
time with `--circuit` is the lever that reduces it; station count does not.

Hot-wall geometry is resolved in Python, not reconstructed from coolant
boundaries in FeatureScript. Squared gas-facing pieces on the same contour
join directly, including the wall over cooling ribs. Squared/rounded and
rounded/rounded neighbors are trimmed and joined with cubic Bézier curves.
The endpoint handles follow the adjoining tangents and have length one third
of the endpoint separation. `--hot-wall-trim-fraction` is the fraction of the
original arc's parameter span removed from **each blended side**, in `(0, 0.5)`;
its default is 0.1. It is not a fraction of physical elliptical arc length.
Profiles are checked for closure and sampled angular folds. Tangency is exact
within each section, but a free loft does not guarantee 3D tangency or clearance
between stations; inspect the generated CAD before relying on wall thickness.

The hot wall starts at the first selected circuit and ends at the furthest
selected circuit, rather than extending to the contour ends. Internal gaps
are bridged by continuing the same loft onto the next circuit's profiles.
The coverage diagnostic uses the union of selected circuit **axial** spans,
including missing coverage at either contour end. It warns if more than 10%
of the full contour's axial length is uncovered; interleaving spans are not
counted twice. This is an axial coverage diagnostic, not angular coverage.

Hot-wall transitions have their own stations, independent of simulation
resolution. `--hot-wall-transition-fraction` defaults to 0.02 of the full
contour's axial length (`contour.xs[-1] - contour.xs[0]`). Alternatively,
`--hot-wall-transition-length-m` sets an absolute length in metres. The CLI
accepts one or the other; the Python API gives an explicit metre value precedence.
A transition is centered on each internal circuit hand-off. An uncovered gap
is one hand-off centered on that gap, with width at least the gap length.
The export inserts profiles at both transition limits and removes interior
hot-wall profiles. Coolant loft sampling is unchanged. Overlapping intervals
or intervals extending outside the overall circuit extent raise an error;
reduce the requested transition length instead of silently clipping it.

Python callers can pass `include_hot_wall=True`, `hot_wall_trim_fraction`,
`hot_wall_transition_fraction`, and `hot_wall_transition_length_m` to
`build_onshape_data` or `export_onshape_json` in the example exporter.

The simulation's flow areas, centroids, wall conduction, and effective hot-gas
area multiplier remain unchanged. CAD blending can change the actual gas area;
it does not reposition tubes or feed that change back into the solver. This
initial export follows the existing gas-facing section/liner geometry. Hot-side
ribs and an explicitly smooth silver-insert region are future extensions.
