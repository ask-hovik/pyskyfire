# Design and solve a volute

A `Volute` collects a swirling flow around a circumference and delivers it to a
round pipe. It is the shaped manifold on a pump discharge, and — through a
different area law — the plain constant-area toroidal manifold sometimes used
on a cooling-circuit inlet.

Three pieces are combined: a **cross-section family**, the **inlet conditions**,
and an **area law** that decides how big the section is at each azimuth.

## Define the cross-section

The passage cross-section is an elliptical bore fed through a parallel-sided
slot, with the wall-to-bore junction rounded on both sides. The slot width is
fixed by whatever feeds the volute, so the section grows only through its
bore size:

```python
import pyskyfire as psf

section = psf.ducts.EllipticalSlottedSection(
    inlet_width=0.012,    # slot width b5, fixed by the pump discharge
    fillet_radius=0.003,
    throat_length=0.002,  # straight wall before the fillet starts
)
```

Each fillet is tangent to its slot wall and to the bore, which is enough to
place the bore centre for any bore size. The only requirement is a bore at
least as wide as the slot feeding it, so `section.minimum_size` is
`inlet_width / 2`.

### Aspect ratio

`aspect_ratio` is the bore's semi-axis across the slot over its semi-axis
outward from the inlet plane. It defaults to one, which is a circle. Above one
the bore is stretched across the slot and squashed radially, keeping the volute
radially compact; below one it reaches further out than it is wide:

```python
section = psf.ducts.EllipticalSlottedSection(
    inlet_width=0.012, fillet_radius=0.003, throat_length=0.002,
    aspect_ratio=1.6,       # oblong across the slot
    max_aspect_ratio=5.0,   # how far it may be flattened at the tongue
)
```

`size` is always the semi-axis *across* the slot, whatever the aspect ratio, so
it plays exactly the role the bore radius does for a circle and the slot still
fixes its floor at `inlet_width / 2`.

The bore can reach inboard of the inlet plane when it is small, which makes an
intermediate rather than a fully external volute. That is valid geometry that
may simply foul neighbouring hardware, so it is reported rather than rejected:

```python
section.inboard_overhang(0.008)   # metres inboard of the inlet plane, 0.0 if clear
```

Lengthening `throat_length` pushes the bore back outboard.

## Define the inlet and build the volute

The inlet carries the radius at which flow is collected and its two velocity
components. These are design inputs: run the pump sizing first and pass its
outputs in here.

```python
inlet = psf.ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)

volute = psf.ducts.Volute(section, inlet)
```

An inlet plane need not be radial. A pump whose discharge points a little
downward, or a manifold sitting against a sloping nozzle wall, has a tilt:

```python
volute = psf.ducts.Volute(section, inlet, tilt=math.radians(45.0))
```

Tilt is constant around the circumference and rotates the section plane
rigidly, so it decides where the passage sits in space and changes nothing in
the loss analysis.

The default area law conserves the angular momentum being collected, which is
Aungier's method of choice. Geometry is resolved when the volute is built, so
one volute can afterwards be solved at several operating points:

```python
geometry = volute.geometry
geometry.area[-1]              # full-collection-plane area [m^2]
geometry.centroid_radius[-1]   # its centroid radius r6 [m]
```

### The tongue

A shaped volute wants zero area at the tongue, and no bore that has to span the
slot can shrink that far. Its *shape* is not stuck, though, and the aspect
ratio is the way past the floor: holding the size at its minimum and flattening
the bore shrinks the area roughly as one over the aspect ratio, without ever
narrowing it below the slot. `Volute` does that automatically, up to
`max_aspect_ratio`, so the passage tracks the area law considerably further in
than a fixed circular bore could:

```python
geometry.aspect               # resolved aspect ratio, per station
geometry.flattened            # boolean: where it was pushed past the design value
geometry.clipped              # boolean: where even the flattest section was too big
geometry.tongue_clip_angle    # azimuth up to which the law was clipped [rad]
```

Only where the flattest allowed section is still too big does the law get
clipped, and it is clipped to that flattest section rather than to the
design-aspect one. The passage right at the tongue is therefore still larger
than conservation of angular momentum asks for, which is also what real
hardware does, since a cutwater has finite thickness.

Setting `max_aspect_ratio` equal to `aspect_ratio` disables flattening, leaving
the fixed-shape behaviour of a family with a hard area floor.

## Solve for pressure loss

```python
result = volute.solve(
    inlet_state=state,
    mass_flow=27.14,
    transport=transport,
    exit_diameter=0.050,
    roughness=15e-6,
)

result.SP                 # realised sizing parameter
result.loss_breakdown()   # coefficients by mechanism
result.pressure_loss      # [Pa]
```

The four loss coefficients are all referred to the inlet velocity head, so the
total pressure loss is `omega_total * rho * C_5**2 / 2`:

| Mechanism | What it charges |
| --- | --- |
| `meridional` | The through-slot velocity head, which becomes swirl in the passage and is never recovered |
| `tangential` | Departure from conserved angular momentum between inlet and collection plane |
| `friction` | Wall shear over the mean spiral path |
| `exit_cone` | The change of velocity delivering to the round pipe |

The meridional term is usually the largest, and nothing the volute does can
change it — it is set by `C_m` at the inlet. Reducing it is a job for the pump
and diffuser design upstream.

## Sizing parameter

`SP` is the ratio of actual passage area to the area that exactly conserves
angular momentum. It is both the design knob and the parameter the tangential
loss is computed from, so a designed volute realises the `SP` it was built at.

The default `SP = 1.0` zeroes the tangential loss. That is not quite the
minimum of the *total* loss: a slightly larger volute slows the flow at the
collection plane, which buys back more in friction and exit-cone loss than the
tangential term costs. Aungier recommends 1.0 to 1.2, biased upward because
undersizing is penalised far more steeply than oversizing.

```python
volute = psf.ducts.Volute(
    section, inlet, psf.ducts.ConstantAngularMomentum(SP=1.15)
)
```

## Outlet ducts

A volute's full-collection plane connects to a pipe flange whose diameter is
usually fixed by a piping standard, in a direction that depends on how the
machine is packaged: some outlets run straight on, some curve round until they
point radially outward, and real hardware sits anywhere between. `ExitCone`
describes that duct:

```python
volute = psf.ducts.Volute(
    section, inlet,
    exit_cone=psf.ducts.ExitCone(
        final_diameter=0.050,
        first_angle="radial",       # or "tangential", or a number
        second_angle=0.0,           # elevation out of the volute's plane
        arc_length_parameter=2.0,   # duct length, in outlet diameters
    ),
)
```

`first_angle` turns within the volute's own plane, measured from *tangential*
-- the direction the flow is already travelling as it leaves the
full-collection plane. That reference is fixed at the volute outlet, since by
the time the duct reaches its flange it has moved round the volute and there
is no tangential direction to measure against there anymore. `"radial"` is
different: it always means straight out from the machine axis, evaluated at
wherever the flange ends up, so it is closer than a right angle away from
tangential -- exactly by however far the duct has travelled round. A numeric
first angle beyond that bound would point the duct back inward, and is
rejected.

`second_angle` lifts the outlet out of the volute's own plane, on top of
whatever the first angle already set. A volute with no exit cone stops at the
full-collection plane, and `Volute.solve()` then needs an explicit
`exit_diameter`; with one attached, the flange diameter is stated once and used
for both the geometry and the loss analysis's station 7.

`ExitCone` geometry does not feed back into the loss analysis beyond that
diameter. Aungier's exit-cone loss (Eq. 5-110) is already a lumped allowance
for the whole transition, including the turn out of the volute's log-spiral
path, so a tangential duct and a radial one at the same flange diameter are
charged the same loss here, whatever their length or how sharply they turn.

## Visualise it

`Volute.sweep()` places each station's section in machine coordinates, with the
axis along +z and the tongue at azimuth zero. `VoluteComponent` wraps that into
a scene. With no exit cone the sweep is capped at both ends, so the step
between the two faces is the cutwater. With one attached, `Volute.sweep_exit_cone()`
continues the same tube seamlessly out to the flange -- its first ring is
bit-identical to the volute's own last one -- and the cap moves there instead:

```python
scene = psf.viz.Scene(title="Volute", color_by="velocity")
scene.add(psf.viz.VoluteComponent(volute, result), name="volute")
scene.save_html("volute.html")
```

Passing the result colours the passage by any per-station field: through-flow
velocity, passage area, bore radius, centroid radius, or the fraction of mass
flow collected so far. The exit cone is not resolved station by station, so its
rings simply carry the full-collection-plane value forward.

Note that drawing a section in its own `(u, v)` coordinates hides the tilt,
because tilt *is* a rotation of that frame. To see it, map into the meridional
plane — `examples/volute/volute_report.py` does this.

## Constant-area collectors

Swap the area law for a uniform passage — the unshaped manifold:

```python
volute = psf.ducts.Volute(section, inlet, psf.ducts.ConstantArea(area=7.8e-4))
```

A collector tries to diffuse the entering flow toward rest where collection
begins, so it is charged an extra half of the tangential velocity head that a
shaped volute avoids. At equal collection-plane area it is substantially lossier
than the shaped volute, which is the trade against its lower manufacturing cost.

## In an engine network

`VoluteBlock` adapts the solver to the network, exactly as `PipeBlock` does:

```python
block = psf.network.VoluteBlock(
    "discharge", "pump_out", "jacket_in", volute, transport, roughness=15e-6,
)
```

`exit_diameter` is optional here too, and defaults from an attached `ExitCone`
the same way `Volute.solve()` does; pass it explicitly to override, or when the
volute has no exit cone.

## Scope

The model is Aungier's one-dimensional inlet-to-collection-plane analysis. It
does **not** resolve how mass flow distributes around the circumference, so it
cannot yet tell you whether cooling channels far from the feed point are
starved relative to near ones. It also requires incoming swirl: a manifold fed
axially through a pipe, with no tangential velocity to conserve, is outside
what these correlations describe.

## References

Aungier, R. H., *Centrifugal Compressors: A Strategy for Aerodynamic Design and
Analysis*, ASME Press, 2000. Sections 5.7 and 11.2.

Weber, C. R. and Koronowski, M. E., "Meanline Performance Prediction of Volutes
in Centrifugal Compressors," ASME 86-GT-216, 1986.
