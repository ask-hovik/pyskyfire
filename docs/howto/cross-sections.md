# Define chamber and heat-exchanger cross sections

`CrossSectionRounded` and `CrossSectionSquared` support both chamber sectors
and independently dimensioned tubes. The geometry, hydraulic properties, and
thermal surface selection use the same section objects.

## Dimensions and wall thickness

All widths describe the **internal fluid passage**: the straight-line distance
between the two sidewall endpoints at the floor or roof. A rounded width is a
lobe chord, not its diameter or arc length. Every angle in these inputs uses
radians.

| Input | Meaning |
| --- | --- |
| `floor_width`, `roof_width` | Internal endpoint-to-endpoint widths, metres |
| `height` | Rounded bottom-to-top height, or squared floor-to-roof reference spacing |
| `wall_thickness` | Rounded uniform wall thickness; squared floor thickness |
| `sidewall_half_thickness` | Squared thickness of **one** complete sidewall; half the full shared rib between adjacent channels |
| `roof_thickness` | Squared roof thickness; defaults to the floor thickness |
| `radius`, `angle` | Alternative occupied-sector construction; angle includes the material allotted to the sidewalls |

An existing 1 mm full squared rib becomes
`CrossSectionSquared(sidewall_half_thickness=0.5e-3)`. The old full-rib keyword
has been removed; convert constants and callable return values by dividing by
two. The separate internal rounded class has also been removed. Select the
desired thermal surfaces on `CrossSectionRounded` instead.

## Local exchanger profiles

Use `SectionProfiles.from_dimensions` without a thrust chamber or placed
path. Scalars and one-dimensional arrays broadcast over the supplied stations:

```python
import numpy as np
from pyskyfire.regen import CrossSectionRounded, CrossSectionSquared, SectionProfiles

profiles = SectionProfiles.from_dimensions(
    stations=np.linspace(0, 0.2, 11),
    floor_width=0.006,
    roof_width=0.006,
    height=0.010,
    wall_thickness=0.0005,
)

rounded = CrossSectionRounded(exposure={"floor": 1, "roof": 1})
rectangle = CrossSectionSquared(
    sidewall_half_thickness=0.0004,
    roof_thickness=0.0008,
    floor="line",
    roof="line",
    exposure={"floor": 1, "roof": 1},
)

flow_area = rounded.A_coolant(profiles)
hydraulic_diameter = rounded.Dh_coolant(profiles)
outer_area = rectangle.A_exterior(profiles)
floor_width, roof_width = rounded.widths(profiles)
```

`A_exterior` includes both wall material and enclosed fluid. Inner and outer
curves are available through `boundary`, `outer_boundary`, and the named
`surfaces` mapping. Their coordinates are in the local section plane; placement
onto a path is a separate operation.

Equal rounded widths give parallel sidewalls. At full height equal to width,
the passage becomes a circle. Unequal widths give different circular lobes
joined tangentially. Some width/height combinations are impossible and raise
`ValueError`. If two tangent constructions are possible, full-height input
chooses the one with the smaller sidewall angle. `height_input="base"` instead
interprets height as the distance between lobe centres and determines that
angle directly.

For squared profiles, choose `floor="arc"` or `"line"` independently of `roof`.
With direct widths and finite `radius`, the floor arc radius is `radius` and
the roof arc radius is `radius + height`; their midpoints are separated by
`height`. Straight ends use the same widths and midpoint spacing. Omitting
radius, or setting it to infinity, produces exact lines. No large artificial
radius is needed to approximate a rectangle. A batch of arc profiles must use
either all finite or all infinite radii.

Squared exterior corners are intersections of the normally offset floor,
roof, and sidewalls. Rounded exterior surfaces all use the same uniform wall
offset, retaining tangency.

## Occupied sectors and existing chamber construction

Alternatively, supply `radius`, `angle`, and `height`, without either width:

```python
sector = SectionProfiles.from_dimensions(
    radius=0.10,
    angle=0.10,
    height=0.012,
    wall_thickness=0.001,
)
section = CrossSectionSquared(sidewall_half_thickness=0.0005)
floor_width, roof_width = section.widths(sector)
```

The sector construction subtracts the sidewall material before resolving the
internal widths. Rounded channels use the uniform wall thickness; squared
channels use `sidewall_half_thickness`. Both routes then construct the same
endpoint-based section geometry as direct widths.

As in chamber profiles, sector `radius` locates the internal floor midpoint
for squared sections and the lower lobe centre for rounded sections. It is
independent of where an exchanger will eventually be placed in world space.

For sector-based squared channels, switching an end from arc to line joins
the same sidewall endpoints with a chord. Thus it preserves the existing
sidewall locations and material allocation. The input height remains the
spacing between the floor and roof reference radii; `full_height(profiles)`
reports the actual spacing between the selected curves' midpoints.

`CoolingCircuit` continues to supply the section profiles from its contour,
placement, channel height, and wall stack. Rounded profiles retain the exact
sidewall angle specified by this construction, including circular limits.
The established projection of helical squared channels is retained; see
{doc}`../explanations/channel-construction` for its material-thickness convention.

## Select thermal surfaces

The four surface names are `floor`, `roof`, `side_left`, and `side_right`.
Each accepts a centred exposed fraction from 0 to 1:

```python
section = CrossSectionSquared(
    sidewall_half_thickness=0.0005,
    exposure={"floor": 0.75, "roof": 0.5, "side_left": 1.0},
)
```

Unspecified surfaces are inactive. With no exposure arguments the complete
floor is active. An explicit empty mapping disables all direct heat transfer.
These choices change thermal perimeters; flow area and the full wetted
perimeter used in hydraulic diameter are independent of exposure.

Rounded floor and roof arcs may instead specify a centred angle in radians:

```python
chamber_section = CrossSectionRounded(
    exposure_angles={"floor": 112 * np.pi / 180},
)
```

This explicitly reproduces the former 112-degree contact model. The angle is
now user input, and the examples and RL10 validation supply it themselves.
An angle larger than its available lobe arc raises. A surface cannot specify
both a fraction and an angle.

`thermal_perimeters(profiles)` returns each active internal length;
`exterior=True` selects the outer side. `P_coolant` and `P_thermal` provide
the totals used by the existing regen interface. The established squared
chamber floor-area approximation is retained for sector-based squared floor
arcs; the other surfaces use their geometric lengths.

Exposure values and squared sidewall/roof thicknesses can use design-input
callables. For local profiles, `x` is the supplied station coordinate; it is
not a world or chamber position. For example, `exposure={"floor": lambda x:
0.5 if x < 0.1 else 1.0}` changes the active fraction along a tube.

## Symmetric ribs and fins

Set `enable_fin=True` on an independent section to include the existing
symmetric sidewall rib model in `R_coolant_per_len`. It represents half a shared
rib on each side of a channel. The conduction calculation therefore uses twice
the individual sidewall thickness as the full rib thickness.

```python
section = CrossSectionSquared(sidewall_half_thickness=0.0005, enable_fin=True)
resistance = section.R_coolant_per_len(profiles, h_c=1500.0, k_wall=300.0)
```

There is one fin switch for both sides. Any nonzero direct exposure on either
sidewall conflicts with fin modeling and raises. A `CoolingCircuit` supplies
its own existing `enable_fin` setting when calculating resistance, allowing
several circuits to share a section definition with different fin settings.

This change supplies section geometry and thermal surface accounting. A general
two-fluid exchanger solver can consume these interfaces in a subsequent step.
