# Define and solve a pipe

A `Pipe` holds a shared `Path` plus its diameter profile. The path contains public `PathStation` objects with one explicit
connection between every adjacent pair. The solver uses these objects directly;
it does not infer whether a corner is a continuous curve or a lumped bend.

## Generate stations from human-friendly corner nodes

`generate_path_stations()` turns control vertices into tangent points. Here the
middle node is the intersection of the two straight control legs, not a
physical fluid station:

```python
import pyskyfire as psf

nodes = [
    psf.ducts.PathNode((0.0, 0.0, 0.0)),
    psf.ducts.PathNode(
        (0.0, 0.0, 0.5),
        bend=psf.ducts.RoundedBend(radius=0.050),
    ),
    psf.ducts.PathNode((0.0, 0.5, 0.75)),
]

path = psf.ducts.Path.from_nodes(nodes)
pipe = psf.ducts.Pipe(path, diameter=0.010, roughness=15e-6)
```

The generated `RoundedBend` is one complete lumped loss evaluated with
`fluids.fittings.bend_rounded`; its arc is not also integrated as distributed
friction.

To use the cooling-channel-style continuous-curvature model instead, replace
the middle-node bend instruction:

```python
bend=psf.ducts.Continuous(radius=0.050, segments=20)
```

This generates twenty physical arc segments with exact radius and arc length.
Resolution is therefore explicit and controlled by the user.

The generator rejects impossible radii, reversing legs, unresolved corners,
and neighboring bends whose tangent trims overlap. Start and end nodes cannot
carry bends.

## Supply the calculation objects directly

The generator is optional:

```python
stations = [
    psf.ducts.PathStation((0.0, 0.0, 0.0)),
    psf.ducts.PathStation((0.0, 0.0, 0.5)),
]
connections = [psf.ducts.Straight()]
path = psf.ducts.Path(stations, connections)
pipe = psf.ducts.Pipe(path, diameter=[0.010, 0.012], roughness=15e-6)
```

Use `Continuous(radius=...)` when the curvature is known. If its radius is
omitted, the solver estimates unsigned Menger curvature from neighboring
station triples. This estimation is only performed because the user explicitly
selected `Continuous`; the solver never switches connection models itself.

Lumped rounded-bend correlations currently require equal diameters at their two
tangent stations. Distributed straight and continuous connections support
changing diameter. Supply a scalar, a function of path distance `s`, or one value
per resolved station on `Pipe(diameter=...)`. Nodes and stations carry geometry only.

## Solve independently

```python
fluid = psf.thermo.Fluid("fuel", ["Water"], [1.0])
transport = psf.thermo.FluidTransport(fluid)
resolver = psf.thermo.FluidStateResolver(transport, p_ref=2.0e6)
inlet = resolver.from_tp(300.0, 2.0e6)

result = pipe.solve(
    inlet_state=inlet,
    mass_flow=0.10,
    transport=transport,
    friction_model="fluids",
)

print(result.total_pressure_loss)
print(result.states[-1])
```

`result` contains pressure, temperature, density, velocity, and fluid state at
every `PathStation`, plus loss, Reynolds number, Darcy factor, and loss
coefficient for every connection.

This first pipe solver is adiabatic and loss-only. Temperature remains constant,
and reversible acceleration from changing area or density is not included.
Diameter variation affects local velocity, Reynolds number, and irreversible
loss, but `result.p` is not yet a complete variable-area static-pressure
solution.

## Use the same pipe in an engine network

```python
block = psf.network.PipeBlock(
    "fuel_feed",
    "tank_out",
    "pump_in",
    pipe,
    transport,
)
```

`PipeBlock` returns the outlet network station, emits `dp_fuel_feed`, and stores
the detailed `PipeResult` in `network.block_results["fuel_feed"]` after network
convergence. Valves and other discrete devices remain separate network blocks.

## Draw it in a 3D scene

```python
scene = psf.viz.Scene(color_by="p_static")
scene.add(psf.viz.PipeComponent(pipe, result), name="fuel_feed")
scene.save_html("duct.html")
```

`PipeComponent` reads the pipe's own geometry and reopens lumped `RoundedBend`
connections into arcs for drawing only. Passing the `PipeResult` colours the
tube by any solved quantity. See [Assemble a 3D scene](3d-scenes.md) for
placement against a chamber or a pump.

## A worked duct

`validation/RL10/plumbing/` builds real RL10 tubing from the fabrication
drawings, one module per component. `bending_table.py` is the reusable part:
it converts a tube-bending table -- straight lengths, bend angles, bend radii,
and the dihedral rotation between consecutive bend planes -- into `PathNode`
control vertices, absorbing the tangent distance `R * tan(angle / 2)` into
both legs at each corner and carrying bend planes by rolling the bend binormal
about the intervening straight. `turbine_discharge_duct.py` is one such table,
with the bore taken from the drawing's OD and wall callouts.

Those pipes are built but not yet wired into the cycle in `sizer_sim.py`,
which still represents each duct as a `SimpleDuctBlock` on a calibrated
pressure ratio. A ratio covers a whole network leg, so it can only be replaced
once every component along that leg exists.
