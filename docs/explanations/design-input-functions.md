# Design-Input Functions

Most of the geometry you hand to Pyskyfire can be either a number or a
function. A wall can be 1 mm thick everywhere, or it can thin out towards the
nozzle exit. A rib can be 1 mm wide everywhere, or it can always take up a
fifth of the local channel pitch.

A function used this way is called a **design input**. It declares the
quantities it needs *by parameter name*, and Pyskyfire supplies exactly those:

```python
sidewall_half_thickness = lambda: 0.5e-3                    # a 1 mm shared rib
sidewall_half_thickness = lambda x: 0.5e-3 + 0.5e-5 * x     # grows downstream
sidewall_half_thickness = lambda pitch_mid: 0.1 * pitch_mid # shared rib takes 20% of pitch
```

The parameter names are the interface. Order does not matter, and a design
input that does not depend on position takes no parameters at all rather than
an ignored `x`.

## Where design inputs are accepted

| Argument | Passed to | Tiers available |
| --- | --- | --- |
| `thickness` | {py:class}`~pyskyfire.regen.thrust_chamber.Wall` | A |
| `channel_height` | {py:class}`~pyskyfire.regen.thrust_chamber.CoolingCircuit` | A, B |
| `helix_angle` | {py:class}`~pyskyfire.regen.channel_placement.SurfacePlacement`, {py:class}`~pyskyfire.regen.channel_placement.InternalPlacement` | A, B, and `h` from C |
| `channel_angular_width` | {py:class}`~pyskyfire.regen.channel_placement.InternalPlacement` | A, B, and `h` and `gamma` from C |
| `sidewall_half_thickness` | {py:class}`~pyskyfire.regen.cross_section.CrossSectionSquared` | A, B, C |
| `roof_thickness` (`section_thickness` slot) | {py:class}`~pyskyfire.regen.cross_section.CrossSectionSquared` | A, B, C |
| `exposure`, `exposure_angles` values (`surface_exposure` slot) | Cross sections; angle inputs only on rounded floor/roof | A, B, C |
| `roughness` | {py:class}`~pyskyfire.regen.thrust_chamber.CoolingCircuit` | A, B, C |
| `hot_gas_surface_area_multiplier` | {py:class}`~pyskyfire.regen.thrust_chamber.CoolingCircuit` | A, B, C |

Every one of these also accepts a plain `float`, which is the same as a
function that takes no parameters and returns that number.

`sidewall_half_thickness` means one complete sidewall of an isolated tube, or
half of the shared rib between neighbouring chamber channels. Two 0.5 mm
sidewalls form a 1 mm shared rib. Exposure fractions range from 0 to 1;
exposure angles use radians. Local exchanger profiles expose their path
coordinate as `x`; chamber-only quantities require a chamber context.

## Why the tiers exist

The available quantities differ from one argument to the next because the
geometry is built in a fixed order, and each step consumes what the previous
ones produced:

```
contour  →  wall stack  →  channel heights  →  centrelines  →  angular widths  →  section
```

A design input can only be handed quantities that already exist when its own
step runs. Asking for one from further down the chain would be circular: a
helix angle that depended on the channel pitch could not be evaluated, because
the pitch is `2π·cos(γ)/n` and so depends on the helix angle. 

The quantities are grouped into three tiers accordingly.

**Tier A — the contour.** Where you are on the chamber wall and what the wall
is doing there. Known before any channel exists, so every design input can use
tier A.

**Tier B — layout and walls.** How many channels there are and how thick the
wall stack is. Known once the circuits are configured. Available everywhere
except `Wall.thickness`, which is what defines the wall stack in the first
place.

**Tier C — the resolved section.** The channel's height, helix angle, pitch and
radii. Known only after heights, centrelines and angular widths have been
built, so it is available to `sidewall_half_thickness`, `roughness` and
`hot_gas_surface_area_multiplier`. The two slots evaluated part of the way
through see part of it: `helix_angle` gets `h`, because heights are built
first, and `channel_angular_width` gets `h` and `gamma`. Neither can see the
pitch quantities they are themselves about to determine.

A placement may be shared by several circuits, and each is evaluated against
its own context. Two circuits sharing one placement but carrying different
channel heights will therefore get different helix angles from the same law —
which the analytical channel layout rejects, since it needs the circuits it
interleaves to follow one angular evolution. If you share a placement, keep its
design inputs to quantities the circuits agree on.

## Examples

A wall that thins from 1.2 mm in the chamber to 0.8 mm at the exit, using the
normalised axial coordinate rather than absolute metres:

```python
psf.regen.Wall(
    material=material,
    thickness=lambda x_frac: 1.2e-3 - 0.4e-3 * max(x_frac, 0.0),
)
```

A rib of constant thickness, which is what `CrossSectionSquared` is built
around — the rib keeps this width from the channel floor to the channel roof,
so the channel itself widens with radius while the rib does not:

```python
psf.regen.CrossSectionSquared(sidewall_half_thickness=0.5e-3)
```

A rib that instead always takes a fixed share of the pitch, so that the channel
and the rib scale together as the chamber radius changes:

```python
psf.regen.CrossSectionSquared(sidewall_half_thickness=lambda pitch_mid: 0.05 * pitch_mid)
```

This last form gives a shared rib occupying 10 % of the mid-height pitch.
It removes approximately a tenth of the coolant area; the section integrates
the actual curved boundary rather than approximating it as a sector minus a rectangle.

A rib sized from the raw circumference and the total channel count instead,
which gives the same answer for a non-helical channel and a slightly different
one for a helical one, because `pitch` is measured normal to the channel axis
while `circumference` is not:

```python
psf.regen.CrossSectionSquared(
    sidewall_half_thickness=lambda circumference, n_channels_all_circuits: (
        0.1 * circumference / n_channels_all_circuits
    ),
)
```

A channel height keyed to the wall it sits on:

```python
psf.regen.CoolingCircuit(
    ...,
    channel_height=lambda t_wall: 3.0 * t_wall,
)
```

A constant helix angle, which needs no parameters:

```python
psf.regen.SurfacePlacement(
    n_channel_positions=90,
    helix_angle=lambda: np.radians(30.0),
)
```

An **area-preserving helix**. A squared section holds

```
A ≈ h * [(pi * cos(gamma) / n) * (2 * r_base + h) - 2 * sidewall_half_thickness]
```

so with the channel height, rib width and channel count fixed, the coolant area
is constant exactly when `cos(gamma) * (2 * r_base + h)` is. Anchoring that
product at one station pins the helix there and opens it out as the radius
grows, instead of letting the passage swell with the wall it wraps:

```python
reference = np.cos(np.radians(20.0)) * (2.0 * (r_start + t_wall) + height)

def area_preserving_helix(r, t_wall):
    return np.arccos(
        np.clip(reference / (2.0 * (r + t_wall) + height), -1.0, 1.0)
    )
```

`examples/advanced/sizer_sim.py` uses exactly this on its nozzle circuits. Note
that it closes over a fixed `height` rather than taking `h`, because the
placement is shared; a single-circuit placement can take `h` directly. Where
the channel count changes — under an overlapping return pass, say — the pitch
changes with it and no helix law can hold the area across that step.

A roughness that grows with the channel it roughens, which is a tier C
dependency and therefore only legal in a slot that is evaluated last:

```python
psf.regen.CoolingCircuit(..., roughness=lambda h: 0.002 * h)
```

## Errors

Naming a quantity that does not exist, or one that is not available in that
slot, raises {py:class}`~pyskyfire.regen.design_inputs.DesignInputError` when
the owning object is constructed. The message names the offending parameter,
says why it was refused, and lists what is legal there:

```text
CrossSectionSquared.sidewall_half_thickness: 'radius' is not a design-input quantity.
Available here: area_ratio, circumference, contour, dr_dx, gamma, h, ...
```

```text
Wall.thickness: 't_wall' is a tier B quantity describing the wall stack,
which is what this function defines; requesting it would be circular.
Available here: area_ratio, circumference, contour, dr_dx, ...
```

Design inputs are evaluated one station at a time and must return a scalar. A
function is therefore free to branch on position — `if x > 0.1:` — without
worrying about array semantics.

Two signature forms are rejected outright. `*args` is refused because design
inputs are called with keyword arguments, so a variadic function cannot say
what it wants; write out the names instead. Positional-only parameters are
refused for the same reason. A `**kwargs` catch-all *is* accepted and receives
every quantity available in that slot, which is occasionally useful for
debugging but makes for an opaque input.

## Glossary

Every quantity a design input may name, exhaustively.

### Tier A — contour

| Name | Meaning |
| --- | --- |
| `x` | Axial coordinate measured along the engine axis [m]. |
| `r` | Hot-wall (contour) radius at `x` [m]. |
| `circumference` | Hot-wall circumference `2*pi*r` at `x` [m]. |
| `dr_dx` | Contour slope `dr/dx` at `x` [-]. |
| `normal_angle` | Angle between the outward wall normal and the plane normal to the engine axis [rad]. |
| `x_frac` | Signed normalised axial coordinate: −1 at the chamber start, 0 at the throat, +1 at the nozzle exit [-]. |
| `area_ratio` | Local area ratio `(r/r_throat)**2` [-]. |
| `x_throat` | Axial coordinate of the throat [m]. |
| `r_throat` | Throat radius [m]. |
| `x_start` | Axial coordinate of the contour start [m]. |
| `x_end` | Axial coordinate of the contour end [m]. |
| `r_chamber` | Contour radius at the chamber start [m]. |
| `r_exit` | Contour radius at the nozzle exit [m]. |
| `contour` | The {py:class}`~pyskyfire.regen.contour.Contour` object itself, as an escape hatch for quantities not listed here. |

### Tier B — layout and walls

| Name | Meaning |
| --- | --- |
| `n_channel_positions` | Number of circumferential channel positions (leaves) in this circuit [-]. |
| `n_channels` | Total number of channels in this circuit, i.e. leaves times channels per leaf [-]. |
| `n_channels_all_circuits` | Number of surface-occluding channels present at `x`, summed over every circuit [-]. |
| `t_wall` | Combined thickness of all main wall layers forming the channel walls, excluding liners [m]. |
| `t_wall_total` | Combined thickness of the whole wall stack [m]. |
| `t_hot_layer` | Combined thickness of layers explicitly marked as liners on the combustion side of the main walls [m]. |

`n_channels_all_circuits` is the one quantity a circuit cannot work out alone,
because it counts channels belonging to other circuits. It is available only
once the circuit has been given to a
{py:class}`~pyskyfire.regen.thrust_chamber.ThrustChamber`.

### Tier C — the resolved section

| Name | Meaning |
| --- | --- |
| `h` | Supplied channel height [m], following the cross section's `height_input` convention: full height by default, optionally rounded base height. |
| `gamma` | Local helix angle of the channel, measured from the meridional direction [rad]. |
| `theta` | Angular pitch of one channel measured in the plane normal to the channel axis [rad]. Already reduced by `cos(gamma)`. |
| `r_channel_base` | Radius of the cross-section base point: the channel floor for a squared section, the inner lobe centre for a rounded one [m]. |
| `r_channel_top` | Radius of the channel roof, `r_channel_base + h` [m]. |
| `pitch` | Channel-to-channel pitch at the hot wall, `r * theta` [m]. |
| `pitch_base` | Channel-to-channel pitch at the channel floor, `r_channel_base * theta` [m]. |
| `pitch_mid` | Channel-to-channel pitch at mid channel height, `(r_channel_base + h/2) * theta` [m]. |

For rounded sections, `height_input="base"` restores the older base-height
parameter. Design-input functions requesting `h` receive that supplied value;
the algebraic `r_channel_top` and `pitch_mid` quantities above use it too.
The height report defaults to the converted full internal height. See
{doc}`channel-construction` for conversion and plotting examples. Keeping `h`
in its input convention also lets helix laws use it before the angular pitch
needed for full-height conversion has been resolved.

`r_channel_base` is resolved geometry, not `r + t_wall`: the section stands off
the contour along its normal, so a sloping wall projects the stack outward, and
a rounded section's base point is the lobe centre its gas-area constraint
solves for. Before the section has been placed — while a height or helix law is
being evaluated — it falls back to the contour radius.

The three pitch quantities differ only in the radius at which the pitch is
measured, and all three are measured in the plane normal to the channel axis,
so they shrink with `cos(gamma)` on a helical channel. Which one to use depends
on what the rib is being sized against: `pitch_base` is the material available
at the rib root, while `pitch_mid` is the one that makes a proportional rib
remove exactly that proportion of the coolant area.
