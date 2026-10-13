# Build a path-based heat exchanger

`pyskyfire.heat_exchanger` supports both prescribed surrounding gas and two
fluids coupled through a wall. Channels use the same `ducts.Path`, `PathNode`,
and `PathStation` objects as pipes, and the existing rounded/squared sections.
All geometry angles are radians. Lengths are metres.

## Define and place a path

```python
import numpy as np
import pyskyfire as psf

path = psf.ducts.Path.from_nodes([
    psf.ducts.PathNode((0, 0, 0)),
    psf.ducts.PathNode((0.1, 0.002, 0)),
]).refined(max_step=0.005)

placements = [
    psf.ducts.PathPlacement(
        start_position=(0, 0.05*np.cos(phi), 0.05*np.sin(phi)),
        rotation=psf.ducts.rotation_about_axis((1, 0, 0), phi),
        initial_floor_direction=(0, -np.cos(phi), -np.sin(phi)),
    )
    for phi in np.linspace(0, 2*np.pi, 32, endpoint=False)
]
```

A placement maps the first template station to `start_position`, rotates the
centerline, then establishes the cross-section roll. The initial floor vector
is expressed in assembly coordinates and points from the passage centroid
toward its floor. It is projected perpendicular to the initial tangent. A
parallel or zero vector raises. The remaining frames use minimal-rotation
transport. No central reference line is needed.

The path passes through each passage's **area centroid**, including when its
section varies. Section dimensions are independent of placement. The user
controls whether adjacent copies meet without gaps or overlaps.

`Path.from_nodes` accepts continuous bends on interior control nodes, for
example `PathNode((0.1, 0, 0), bend=Continuous(radius=0.02))`, which becomes
one arc connection. Heat exchangers accept `Straight` and `Continuous`
connections. Lumped `RoundedBend` connections remain available to pipes only.
`Path.from_curve` works too: it samples a spline or any other curve into
exact arc connections (see [Route a pipe along a curve](pipes.md#route-a-pipe-along-a-curve)).
Unlike a pipe, a heat exchanger does not add stations by itself:
`refined(max_step)` sets its thermal integration resolution, and generated arcs
retain exact length, curvature, and bend-plane geometry when split. Manually supplied continuous paths without arc-plane
metadata should already resolve the intended geometry before refinement.

## Construct a channel bank

```python
hx = psf.heat_exchanger
wall = psf.regen.Wall(psf.common.solids.copper, thickness=0.0005)
section = psf.regen.CrossSectionSquared(
    sidewall_half_thickness=0.0005,
    floor="line", roof="line",
    exposure={"floor": 1.0, "roof": 1.0},
)
transport = psf.thermo.FluidTransport(
    psf.thermo.Fluid("fuel", ["Helium"], [1.0])
)
circuit = hx.HeatExchangerCircuit(
    name="helium", path=path, placements=placements,
    cross_section=section, wall=wall, transport=transport,
    dimensions={"height": 0.003, "floor_width": 0.005, "roof_width": 0.005},
    roughness=2e-6,
)
exchanger = hx.HeatExchanger([circuit])
result = exchanger.solve(
    {"helium": hx.Inlet(T=150, p=6e6, mass_flow=0.1)},
    boundaries={"helium": hx.PrescribedGas(temperature=1000, coefficient=500)},
)
```

`mass_flow` is the **total bank flow**. It is split equally among placements.
The solver represents one channel and verifies that all its copies see the
same prescribed gas conditions. Use separate circuits when the thermal
conditions differ. There is no hydraulic flow-balancing solve.

Dimensions accept constants, arrays with one entry per path station, or
functions of arc length `s`. Use `radius` and `angle` instead of the two widths
for sector-defined sections. `Wall.thickness(s)` defines floor thickness;
squared roof and sidewall settings still apply independently. A complete
`SectionProfiles` can instead be passed as `section_profiles=`.

All wetted surfaces count toward hydraulic diameter. Only selected exposure
fractions exchange heat. Partial exposure is centered on each surface. The
section's `enable_fin` switch adds the existing symmetric rib correction;
active sidewalls and rib enhancement are incompatible. Thermal conduction
uses each surface's own thickness and temperature-dependent wall conductivity.
The first solver has local through-wall conduction, without axial or general
circumferential conduction between wall patches.

## Connect two fluids

Both circuits have their own inlet, fluid, path, section and flow rate. A
`ThermalLink` joins selected **segments**, indexed by their upstream station:

```python
n = len(cold.path.stations) - 1
link = hx.ThermalLink(
    circuit_a="cold", circuit_b="hot",
    pairs=tuple((i, n-1-i) for i in range(n)),
    surface_a="floor", surface_b="floor", wall=shared_wall,
)
exchanger = hx.HeatExchanger([cold, hot], thermal_links=[link])
```

Reverse the pairing for counterflow; `(i, i)` gives cocurrent pairing. Paths
are always stored in their own flow direction. Pairing need not follow spatial
proximity: the user states which physical surfaces share a wall.

Each pair transfers exactly opposite heat into the two banks. The shared wall
is counted once. Its thickness must describe the complete physical separator,
including both wall layers if two tubes are bonded together. Contact area
is the smaller available exposed area, or may be supplied explicitly using
`areas=` (m² per pair per copied channel). Area budgets reject reused patches.
Unused exposure can exchange heat with a prescribed boundary. Linked banks
currently require equal copy counts and one-to-one copies.

`FlowConnection("down", "return")` separately connects outlet flow to the next
bank. Supply an inlet only for the first bank. The ideal header transfers all
flow without heat or fitting loss; banks require the same fluid, copy count,
and matching endpoint flow areas. Model area transitions on the actual paths.
Header turn geometry is not included in the swept envelope or pressure loss.

## Install in a thrust chamber

```python
installation = chamber.install_heat_exchanger(
    exchanger, position=(0, 0, 0), rotation=np.eye(3), coefficient=1000,
)
result = installation.solve(inlets)
print(installation.volume_change)
```

Placement converts every station into chamber coordinates. The chamber's
existing aerothermodynamics provides recovery temperature at its axial
position. Gas temperature, chemistry, and flow remain prescribed.

The exterior coefficient can be a constant or
`coefficient(chamber_position, wall_temperature)`. If omitted, the adapter
uses the existing Bartz-style reference-temperature correlation. This is an
engineering approximation for immersed hardware; select an appropriate
coefficient for the actual baffle/duct arrangement. There is no resolved
external blockage loss or radiation model in this first implementation.

`VolumeChange` reports total swept envelope, displaced volume inside the
contour, volume outside it, and `chamber_displaced` upstream of the throat.
It includes the enclosed coolant passages. The sum assumes nonoverlapping
copies; it is not a solid-union calculation. Axial integration uses station
areas and segment lengths. Clipped sections use triangle/axial quadrature;
refine the path/section resolution for tightly curved or intersecting geometry.
The chamber's nominal volume is not overwritten. A thermal solve requires
fully contained envelopes; clipped installations currently support volume
accounting only.

## Pressure and heat-transfer models

Shared scalar functions now live in `common.physics`, `common.friction`, and
`common.momentum`. Pressure includes Darcy friction, the existing distributed
curvature correction, and reversible acceleration. The energy march conserves
`h + velocity**2/2`, using static pressure for fluid-property queries.

Internal turbulent convection uses
[Gnielinski](https://ht.readthedocs.io/en/release/ht.conv_internal.html#ht.conv_internal.turbulent_Gnielinski).
The straight Darcy factor is shared with the friction calculation; the coil
friction multiplier is not automatically treated as a heat-transfer multiplier.
Gnielinski's range is enforced. Properties are bulk values, corrected for a
heated fluid by Gnielinski's gas factor `(T_b/T_w)^0.45`, where `T_w` is the
wall temperature at the fluid interface. It was fitted at moderate temperature
ratios; a cryogen against a hot wall is an extrapolation. A cooled fluid is not
corrected. Pass `property_correction=False` to a circuit for bulk properties.

Below Re=2300, the default is Nu=3.66, with a linear transition to Gnielinski
at Re=4000. **That laminar default is a circular,
fully developed, constant-wall-temperature approximation.** For other boundary
conditions, supply `laminar_nusselt=` or a `heat_transfer` callback accepting
`Re, Pr, k, Dh, darcy_factor`. Hydraulic-diameter use in noncircular passages
and the shared circular-duct laminar friction model are approximations.

The initial solver requires an enthalpy-capable CoolProp backend and rejects
two-phase states. Pressure/energy and coupled heat iterations raise on failure
to converge; they do not return unconverged results.

## Examples and reports

```bash
python examples/heat_exchanger/two_fluid.py
python examples/heat_exchanger/baffles.py
```

The first heats oxygen with hydrogen in counterflow. The second makes a shallow
V-shaped ring pair and six radial spokes with paired straight return channels.
Both produce interactive reports under `examples/heat_exchanger/output/`.

`psf.viz.heat_exchanger_report(exchanger, result, installation=installation)`
includes swept geometry, total/per-channel flow, areas, fluid temperature,
static pressure, velocity, both wall-interface temperatures, and energy balance.
`result.walls` identifies each wall's surfaces and segment indices.

## Put an exchanger and chamber in the same scene

```python
scene = psf.viz.Scene(color_by="T_coolant_wall", show_axes=True)
chamber_node = scene.add(
    psf.viz.ThrustChamberComponent(chamber, layout="analytical"), name="chamber",
)
for group in chamber_node.children:
    for channel in group.children:
        channel.opacity = 0.16
scene.add(psf.viz.HeatExchangerComponent(installation, result), name="baffles")
scene.save_html("baffle-scene.html")
```

Passing a `ChamberInstallation` uses its simulation transform automatically.
A raw `HeatExchanger` instead uses assembly coordinates and supports ordinary
`Scene.add(position=..., rotation=...)` placement. If the chamber itself has a
scene transform, add the installed exchanger with `parent=chamber_node` so
both move together. Explicit scene placement overrides the installation's root
transform, and affects the drawing only.

Each channel is drawn as its wetted passage, the surface the coolant touches,
not the outer envelope of its walls. The adapter exposes fluid temperature,
static pressure, velocity, enthalpy, and `T_coolant_wall`, the wall
temperature at the wetted interface, in the field selector. `T_coolant_wall`
uses contact-area averages on each exposed surface segment; inactive surfaces
stay uncoloured. The opposite interface is not drawn; it remains in
`result.walls`. Each copied channel has inlet/outlet anchors; root
`circuit_N.inlet/outlet` anchors refer to the first copy.

The baffle example writes `baffle-scene.html` and embeds the same scene in its
report. It includes 96 unhelical squared chamber channels from injector to
nozzle exit, shown transparently, and places the baffles 2 mm downstream of the
injector datum in both the simulation and scene. Chamber cooling geometry is
constructed for this example; the existing baffle thermal solve supplies the
wall-temperature colours.
