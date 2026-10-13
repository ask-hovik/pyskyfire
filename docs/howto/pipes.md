# Define and solve a pipe

A `Pipe` in Pyskyfire is an object that transports fluid along a `Path`. A path
is a sequence of `PathStation` points joined by connections, and each
connection says how the pressure loss between its two stations is computed.
The program is set up this way so that a host of different 1D pipe models can
be used interchangeably along the same pipe.

All the code on this page comes from `examples/pipe/pipe_paths.py`, which also
writes every figure shown here.

## Nodes and stations

There are two kinds of point in a pipe definition, and it is worth keeping them
apart:

`PathNode`
: A control vertex: a point you place to describe the route, like the corners
  on a tube-bending drawing. A node may carry a `bend` saying how its corner is
  to be rounded. Nodes are input only. The solver never sees them, and a
  corner node with a bend does not even lie on the pipe, because the bend cuts
  the corner.

`PathStation`
: A point on the pipe's actual centerline. The solver evaluates the fluid state
  at every station, and each pair of neighbouring stations is joined by one
  connection (`Straight`, `RoundedBend` or `Continuous`). A `Path` is made of
  stations and connections only.

There are three ways to lay out a path, and the next three sections take them
in turn:

- [Route a pipe through corners](#route-a-pipe-through-corners): place a few
  `PathNode`s and let `Path.from_nodes()` round their corners into stations
  and connections. This is the easiest way for straight legs and bends, as on a
  tube-bending drawing.
- [Route a pipe along a curve](#route-a-pipe-along-a-curve): let
  `Path.from_curve()` sample a spline or any other curve into stations.
- [Route a pipe station by station](#route-a-pipe-station-by-station): build
  the `PathStation`s and connections yourself.

## Route a pipe through corners

`Path.from_nodes()` turns a list of `PathNode`s into stations and connections
with `generate_path_stations()`. Each corner is rounded by the bend its node
carries, and the bend can book its loss in one of two ways: as one lumped loss
for the whole bend, or as wall friction distributed along its arc.

### Lumped bends: `RoundedBend`

Here is how to define the nodes and use the `bend` attribute:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:rounded
:end-before: tutorial:end:rounded
:dedent: 4
```

So what happened above? We placed three `PathNode`s, at $(0, 0, 0)$,
$(0, 0, 0.5)$ and $(0, 0.5, 0.75)$. On the second node we specified a
`RoundedBend` with a radius of 0.25 m, so the two straight legs that meet at
$(0, 0, 0.5)$ are joined by an arc of that radius, tangent to both. The legs
turn through $\theta = 63.4°$, so the generator trims each of them by the
tangent distance $R\tan(\theta/2) = 0.155$ m and places a `PathStation` at
each tangent point. The corner node does not become a station. The resolved
path has four stations, the inlet, the two tangent points and the outlet,
joined by a `Straight`, the `RoundedBend` and another `Straight`. The figure
below shows the three nodes, the four stations generated from them, and the
arc the bend follows as a dashed line:

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/rounded-bend-path.html"
    title="PathNodes, the PathStations generated from them and the rounded-bend arc"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

The generated `RoundedBend` is one complete lumped loss evaluated with
`fluids.fittings.bend_rounded`; its arc is not also integrated as distributed
friction.

Once the `Pipe` is defined, `PipeComponent` draws its geometry in a 3D
[scene](3d-scenes.md):

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:scene-rounded
:end-before: tutorial:end:scene-rounded
:dedent: 4
```

```{raw} html
<div class="psf-demo-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/rounded-bend-pipe.html"
    title="Pipe with a rounded bend"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

Drag to rotate, scroll to zoom and right-drag to pan. The tube is drawn through
the bend's arc, but that arc is for display only: the solver still sees one
lumped connection between the two tangent stations.

### Distributed bends: `Continuous`

The rounded-bend correlations are well established, but they book the whole
bend as one loss coefficient at one diameter, so they require the same bore at
both tangent stations. Say the bore also flares out along the pipe, here from
25 mm to 90 mm. A `RoundedBend` cannot represent that, so we use a
`Continuous` bend instead:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:flared
:end-before: tutorial:end:flared
:dedent: 4
```

On a `PathNode`, a `Continuous` bend is still a circular arc of the given
radius; what changes is how its loss is booked. The generator emits the arc as
one connection that carries its exact radius, arc length and plane, and the
solver charges it distributed wall friction corrected for curvature
(`common.friction.curvature_multiplier`). Because the loss is distributed, the
diameter may change along the bend.

You never choose how many stations a pipe needs. The solver evaluates each
connection at one diameter, the mean of its two ends, so wherever the bore
varies, `Pipe` splits straight and continuous connections by itself until the
diameter changes by no more than 2 % across each piece. That keeps the error
of the one-diameter evaluation below 0.05 %. The extra stations go where the
bore changes fastest, and the arcs stay exact. A constant bore needs none: a
constant-bore arc loses exactly the same as one connection as it does split
into many. The tolerance can be changed with `Pipe(max_diameter_change=...)`.

The figure below puts the two pipes side by side. On the left, the rounded bend
is inferred from the corner: the generator places only the two tangent
stations, and the arc between them is one lumped connection. On the right, the
flaring pipe has been given 67 stations. They are densest at the narrow inlet,
where a millimetre of bore is the largest relative change, and they continue
along the straight legs, because the bore flares there too. Hover over a
station to read its diameter.

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/bend-comparison.html"
    title="Rounded bend and flaring continuous bend compared"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

The flaring pipe in 3D:

```{raw} html
<div class="psf-demo-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/flared-bend-pipe.html"
    title="Flaring pipe with a continuous bend"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

The generator rejects impossible radii, reversing legs, unresolved corners,
and neighboring bends whose tangent trims overlap. Start and end nodes cannot
carry bends.

## Route a pipe along a curve

`Path.from_curve()` makes a path from a curve: a spline, a helix, or any shape
you can write as a function. Here a SciPy `CubicSpline` routes a pipe between
two flanges, one at the origin facing up and one at $(0.5, 0.4, 0.6)$ facing
along $+x$, through one waypoint between them. Its end conditions fix the
direction at each end, so the pipe leaves and enters both flanges square to
them:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:spline
:end-before: tutorial:end:spline
:dedent: 4
```

`from_curve` accepts anything *curve-like*, in the same spirit as NumPy's
"array-like": anything that maps a parameter `t` to a point. Called with an
array of parameters, it returns points of shape `(len(t), 3)`, or `(len(t), 2)`
for a curve in the x-y plane. A function that only takes one `t` at a time and
returns one point also works. SciPy's `CubicSpline`, `make_interp_spline` and
`PchipInterpolator` qualify as they are, and so does a plain function such as a
helix:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:helix
:end-before: tutorial:end:helix
:dedent: 4
```

Output in other forms needs a one-line wrapper. `scipy.interpolate.splprep`,
for example, returns spline data rather than a callable, so pass
`lambda t: np.array(splev(t, tck)).T`. `t_span` is the parameter range to
follow, which you state because every library keeps its range differently.

You do not choose how many stations to take. `from_curve` starts from 16 even
intervals and halves any whose curve turns by more than `max_turn`, 5° by
default, so stations crowd into tight bends and stay sparse along straight
stretches. It also checks the turn across each station, so a sharp corner
cannot hide on a sample point: a curve with a true corner raises an error,
because a pipe cannot follow one without a bend. Each interval then becomes the
circular arc through its two stations and the curve's own midpoint between
them, a `Continuous` connection carrying that arc's radius, length and plane,
or a `Straight` where the three points line up.

The samples become stations, not nodes. A `PathNode` route is straight legs
meeting at corners, each rounded by a bend you specify, so curve samples used
as nodes would need a bend radius at every one and would give straight legs
with short arcs rather than the curve. A station lies on the centerline
itself, so the pipe passes through every sample.

The spline above becomes 39 stations. The figure compares each connection's
curvature with the spline's own, $|\mathbf{r}' \times \mathbf{r}''| / |\mathbf{r}'|^3$.
The pipe turns hardest just after it leaves each flange, with a tightest radius
of 0.26 m, and most gently near the waypoint. The connections follow the
spline to within 0.3 % of the peak curvature, and 0.02 % on average:

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/spline-curvature.html"
    title="Curvature along the spline pipe, fitted and exact"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

```{raw} html
<div class="psf-demo-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/spline-pipe.html"
    title="Pipe following a cubic spline"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

Solved with the water flow of [Solve a pipe](#solve-a-pipe), the default 39
stations give a loss of 7955.0 Pa. Tightening `max_turn` to 1° (172 stations)
or to 0.25° (677 stations) gives 7955.5 Pa, so the default is within 0.01 %.
The length, 0.96172 m, is the same in all three, because each connection
carries its arc length rather than its chord.

## Route a pipe station by station

The two routes above are conveniences. You can also skip them and supply the
`PathStation`s and the connections between them yourself: the exact objects
the solver marches over. Here seven stations are placed by hand on a quarter
circle of 0.2 m radius, between a straight leg up and a straight leg across,
and the arc's six intervals are joined by `Continuous()` connections:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:stations
:end-before: tutorial:end:stations
:dedent: 4
```

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/manual-path.html"
    title="PathStations placed by hand on a quarter circle"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

A `Continuous()` without a radius has its curvature estimated from its
stations and their neighbours, as the circle through each three in a row
(Menger curvature), averaged over the triples the connection belongs to. The
solver only does this because `Continuous` was chosen; it never switches
connection models itself. In the middle of the arc the estimate is exact. At
either end of it, though, one of the triples reaches back onto a straight leg,
and the estimate there comes out far too gentle, 0.77 m and 0.97 m instead of
0.2 m. Where the radius is known, give it:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:stations-radius
:end-before: tutorial:end:stations-radius
:dedent: 4
```

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/manual-curvature.html"
    title="Curvature of each hand-placed connection, estimated and given"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

With the water flow of [Solve a pipe](#solve-a-pipe), the estimated radii give
a loss of 7492 Pa and the given ones 7640 Pa, 2 % more, all of it in the two
end connections of the arc. A `Continuous` without a `length` also takes the
straight distance between its stations, so the arc here is 0.09 % short of its
true length. Pass `length=` as well where that matters, or sample the curve
with `Path.from_curve()`, which gives every connection its exact arc.

The bore is given on `Pipe(diameter=...)` the same way whichever route made
the path: a scalar, a function of arc length `s`, or one value per supplied
station, taken as linear in `s` between them. Nodes and stations carry
geometry only. Where the bore varies, `Pipe` adds the stations it needs, as it
did for the flaring bend. This pipe widens by 20 % along its one connection:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:tapered
:end-before: tutorial:end:tapered
:dedent: 4
```

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/tapered-stations.html"
    title="The two supplied stations and the eleven resolved ones"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

`Pipe` splits it into ten, so `tapered_pipe.stations` holds eleven stations
rather than the two supplied, each piece widening by under 2 %. A
`RoundedBend` correlation books its whole bend at one bore, so it needs equal
diameters at its two tangent stations; `Straight` and `Continuous`
connections can change diameter.

## Solve a pipe

Here the rounded-bend pipe carries 2 kg/s of water, which is 4.1 m/s in the
25 mm bore and a Reynolds number of about $1.2 \times 10^5$:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:solve
:end-before: tutorial:end:solve
:dedent: 4
```

`result` is a `DuctResult`. It contains pressure, temperature, density,
velocity, flow area and fluid state at every `PathStation`, plus loss, Reynolds
number, Darcy factor, and loss coefficient for every connection.

This first pipe solver is adiabatic and loss-only. Temperature remains constant,
and reversible acceleration from changing area or density is not included.
Diameter variation affects local velocity, Reynolds number, and irreversible
loss.

The state it marches is a stagnation state: `inlet_state.p` and
`inlet_state.T` are taken as the inlet stagnation pressure and temperature.
`result.p_stagnation` is that pressure less the irreversible losses up to each
station, and `result.T_stagnation` is the inlet temperature throughout, since
the pipe is adiabatic. `result.p_static` subtracts the local dynamic pressure
`density*velocity**2/2`, and `result.T_static` the dynamic temperature
`velocity**2/(2*cp)`, so the transport must provide `cp` as well as `rho` and
`mu`. Both are first-order corrections, not a complete variable-area
static-state solution.

### Compare a lumped and a distributed bend

To see how the two bend models differ, build the same corner again, with the
same constant 25 mm bore, but as a `Continuous` bend, and solve it with the
same inlet and flow:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:continuous
:end-before: tutorial:end:continuous
:dedent: 4
```

`viz.PlotPipePressureLoss` plots the cumulative stagnation pressure loss of one
or more solved pipes against arc length:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:plot-loss
:end-before: tutorial:end:plot-loss
:dedent: 4
```

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/pressure-loss-comparison.html"
    title="Cumulative pressure loss of the rounded and continuous bends"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

Both pipes have the same four stations and lose the same along the two
straights; they differ only across the bend, where each has one connection.
The totals are close: the rounded bend loses 2576 Pa ($K = 0.310$) and the
continuous bend 2485 Pa ($K = 0.299$). These come from two different
correlations, Rennels' bend formula and curvature-corrected wall friction, and
for a gentle bend like this one, with $R/D = 10$, they agree to within 4 %. The
difference that remains is mostly the extra loss at a bend's entry and exit,
which Rennels' formula includes and the curvature correction, made for flow
that has settled into a long coil, does not.

### Losses where the area changes

Where the flow area changes, the flow loses not only to wall friction but also
to the change. Where the area grows, it loses to its own slowing, as a diffuser
does; above about a 7° full cone angle that loss exceeds the friction. Where
the area shrinks, it loses to its speeding up, as a contraction does: little
for a gentle taper, but about half the downstream velocity head at a sharp
step. `solve` adds both on top of the friction, and
`pipe.area_change_loss_coefficient` holds them for every connection.

A stretch over which the area keeps growing is treated as one diffuser, and
one over which it keeps shrinking as one contraction. Each ends where the area
changes the other way, or at a constant-section stretch longer than four
hydraulic diameters, since the flow recovers there; a shorter one, such as a
bend between two cones, leaves it one diffuser. A diffuser loses the expansion
term of Rennels' conical-diffuser correlation
(`fluids.fittings.diffuser_conical`) for its two end diameters, and a
contraction the convergence term of his conical-contraction correlation
(`contraction_conical`). Each is taken at the cone angle averaged with weight
where the flow slows or speeds up. For a true cone that is Rennels' value, and
for a step his sudden-expansion or sharp-contraction loss. It does not depend
on how finely the pipe is split, and a nearly parallel stretch at either end
of a steep cone does not water it down. `pipe.cone_angle` gives the local full
cone angle of every connection.

This is a one-dimensional estimate. The conical correlations themselves differ
by 30 to 50 %, and a long bore that widens slowly at first and steeply at the
end, like a trumpet bell, can be underestimated: the flow separates at the steep end,
where little of the slowing happens. Where the flow is evidently separated,
use a more detailed method.

### Follow the pressure along a flaring pipe

The flared pipe is solved exactly like the others; its 67 stations are already
in place. `viz.PlotPipePressure` shows its static and stagnation pressure:

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:plot-pressure
:end-before: tutorial:end:plot-pressure
:dedent: 4
```

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/flared-pressure.html"
    title="Static and stagnation pressure along the flared pipe"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

The stagnation pressure only falls, by 827 Pa in all: 687 Pa of wall friction
and 140 Pa of diffuser loss, at the pipe's gentle 3.6° cone angle. The static pressure rises: as the bore widens the flow slows, and its
dynamic pressure, 8.3 kPa at the inlet, is converted back into static pressure.

## Give a pipe a core

A `Pipe` can carry a concentric core inside its bore: a centrebody, a turbine's
nose, or the inner line of coaxial tubing. Give its diameter as
`core_diameter`, the same ways as `diameter`: a constant, one value per
station, or a function of arc length. Zero, the default, is no core.

```python
nodes = [
    psf.ducts.PathNode((0.0, 0.0, 0.0)),
    psf.ducts.PathNode((0.0, 0.0, 0.5), bend=psf.ducts.RoundedBend(radius=0.1)),
    psf.ducts.PathNode((0.0, 0.5, 0.5)),
]
coaxial = psf.ducts.Pipe(psf.ducts.Path.from_nodes(nodes),
                         diameter=0.04, core_diameter=0.02)
```

The section normal to the path is then an annulus. Wall friction acts on its
hydraulic diameter, `D - d`, and laminar friction uses the exact
concentric-annulus `f*Re`, which runs from 64 for a pipe to 96 for parallel
plates (`common.friction.annulus_laminar_factor`). Turbulent friction is left
uncorrected; on the hydraulic diameter it differs by a few percent.
`coaxial.area`, `coaxial.hydraulic_diameter` and `coaxial.r_mean` give the
section at every station.

The core need not run the whole length. It may start from a point and grow,
as a nose does, and it may end, tapering or bluntly. The flow area changes
wherever the core does, so a growing core is a contraction and a core that
ends is a diffuser, and the area-change loss above applies to both: a blunt
core end loses as a sudden expansion. Where the section varies, the stations
are refined on the area and the hydraulic diameter, which both stay finite as
the core shrinks to nothing, so a pointed core needs no special treatment:

```python
path = psf.ducts.Path.from_nodes([psf.ducts.PathNode((0.0, 0.0, 0.0)),
                                  psf.ducts.PathNode((0.2, 0.0, 0.0))])
nosed = psf.ducts.Pipe(path, diameter=0.04,
                       core_diameter=lambda s: 0.03*min(s, 0.1)/0.1)
```

A bend's centreline radius must exceed the bore radius, with or without a
core, or its inner side folds. A `RoundedBend` round a core loads one side of
the annulus and drives secondary flow round the core, and no annular-bend
correlation is applied: the circular-pipe correlation is charged, on the
hydraulic diameter or, with `solve(bend_length_scale="outer")`, on the bore,
with a `PipeCorrelationWarning`. Solving with both brackets the bend loss. The
flow leaving such a bend is uneven round the circumference, which a 1D model
cannot show.

## Give a pipe swirl

A pipe, with or without a core, can carry swirl about its centreline. Give the
circumferential velocity at the inlet's mean radius, `pipe.r_mean`, which is
`(D + d)/4`, and a quarter of the bore without a core:

```python
result = coaxial.solve(inlet_state=inlet, mass_flow=2.0, transport=transport,
                       inlet_swirl=3.0)
print(result.swirl_velocity[-1], result.swirl_angle[-1])
```

Swirl is carried as angular momentum `r*c_u`. A change of radius leaves it
unchanged, because that is reversible; wall friction is what takes it away.
The wall shear lies along the total velocity, so a swirling flow rubs along a
helix that is longer than the pipe by `c/c_m`, and its friction loss grows as
`(c/c_m)**3` on the through-flow head. `result.meridional_velocity` is the
through-flow velocity, `result.swirl_velocity` the swirl at the mean radius,
and `result.velocity` their magnitude, which `p_static` and `T_static` use.
Every `loss_coefficient` is on the through-flow head `rho*c_m**2/2`. The bend
correlations are for flow without swirl, so swirl entering a `RoundedBend`
warns.

## Put a pipe in an engine network

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:network
:end-before: tutorial:end:network
:dedent: 4
```

`PipeBlock` returns the outlet network station, emits `dp_fuel_feed`, and stores
the detailed `DuctResult` in `network.block_results["fuel_feed"]` after network
convergence. Valves and other discrete devices remain separate network blocks.

## Draw a pipe in a 3D scene

```{literalinclude} ../../examples/pipe/pipe_paths.py
:language: python
:start-after: tutorial:start:scene-pressure
:end-before: tutorial:end:scene-pressure
:dedent: 4
```

```{raw} html
<div class="psf-demo-frame">
  <iframe
    src="../_static/howto-artifacts/pipes/flared-pressure-3d.html"
    title="Flared pipe coloured by static pressure"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

`PipeComponent` reads the pipe's own geometry and reopens lumped `RoundedBend`
connections into arcs for drawing only. Passing the `DuctResult` colours the
tube by any solved quantity, here the static pressure. A pipe with a core is
drawn as the core inside a translucent bore, the core only where it runs. See
[Assemble a 3D scene](3d-scenes.md) for placement against a chamber or a pump.

## Worked example: RL10 tubing

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
