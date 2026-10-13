# Design an axial turbine

Blade generation needs pyTurbo, which comes with the impeller extra, and the
scene needs the viz extra:

```sh
pip install 'pyskyfire[turbine,viz]'
```

A turbine is designed in three steps, each with its own object:

1. `TurbineDuty`: what the turbine is asked to do.
2. `AxialTurbine`: a meanline design that meets the duty and draws the
   hardware as a `TurbineDefinition`.
3. `TurbineDefinition.build_geometry()`: the blade surfaces, and what can
   only be measured on them.

```python
from pyskyfire.turbine import AxialTurbine, RowSettings, Stage, TurbineDuty
from pyskyfire.viz import Scene, TurbineComponent

duty = TurbineDuty.from_states(fluid="Methane", mdot=22.8, p0_in=160e5, T0_in=633.,
                               p_out=115.5e5, n=40_000)
turbine = AxialTurbine(duty, efficiency=0.75,
                       stages=[Stage(flow_coefficient=0.6, reaction=0.5)],
                       mean_diameter=0.16,
                       rows=RowSettings(aspect_ratio=1.5, zweifel=0.8, tip_clearance=3e-4))
geometry = turbine.definition.build_geometry()

scene = Scene(title="Turbine")
scene.add(TurbineComponent(geometry))
scene.save_html("turbine.html")
```

## The duty

`TurbineDuty.from_states` takes a medium, a mass flow, the inlet total state,
the outlet pressure and the shaft speed. `TurbineDuty.from_stations` takes
the same from two network stations. `fluid` is a CoolProp name, a
`thermo.Fluid`, or any `thermo.EquationOfState`.

The duty reports the isentropic head `isentropic_head`, `pressure_ratio`, the
spouting velocity `c0 = sqrt(2*Δh_s)`, the inlet and isentropic-outlet flows
and densities, the specific speed `Ns`, and `blade_jet_ratio(D)` and
`specific_diameter(D)` for a candidate diameter.

`outlet_basis` says whether the outlet pressure is total (`'total'`, the
default, and what network stations carry) or static (`'static'`). The
efficiency a design is given is read on the same basis. Total-to-static, the
exit kinetic energy counts as loss, so a design refuses an efficiency that
leaves less loss than the exit velocity carries away.

## The meanline

The design runs at one mean radius through the whole turbine, so the blade
speed `u` is the same at every station. Give either `mean_diameter` or
`blade_jet_ratio` (`u/c0`).

Each `Stage` is a stator and a rotor. Its `flow_coefficient` (`c_x/u`) sets
the axial velocity through it. It also takes one of two things: `reaction`
(the rotor's share of the stage's static enthalpy drop) or `exit_swirl` (the
absolute swirl leaving the rotor, in m/s). The other follows from it.

Each stage does its `work_split` share of `efficiency * Δh_s`. The
efficiency is an input. Where it comes from is up to you, and the design does
not check it against anything. It enters the design twice:

- in the work, and so in how far the rotor turns the flow;
- through the entropy rise it implies, which the design shares among the
  rows (evenly, or by `loss_split`). Each station's static state follows
  from its enthalpy and that entropy, and continuity at that density gives
  the blade height.

`turbine.stations` holds the flow at every row edge: axial and swirl
velocity, absolute and relative angles, static state, total pressure, and
absolute and relative Mach numbers. `turbine.stages` gives each stage's
reaction, flow coefficient, loading and Euler work. `turbine.exit_swirl` and
`turbine.exit_flow_angle` describe the flow a downstream duct or collector
receives.

Angles are radians from axial, positive in the direction the rotor turns.

## From triangles to blades

`RowSettings` applies to every row:

| Setting | Meaning |
|---|---|
| `aspect_ratio` | Mean blade height over true chord. |
| `zweifel` or `pitch_chord` | Blade spacing: from Zweifel's loading coefficient (default 0.8) or directly. |
| `incidence_deg`, `deviation_deg` | Metal angle minus flow angle at the inlet; how far the leaving flow falls short of the outlet metal. |
| `stator`, `rotor` | A `SectionProfile` for each kind of row: see below. |
| `axial_gap` | Space between rows, over the upstream row's axial chord. |
| `inlet_length`, `outlet_length` | Straight annulus before the first row and after the last, over that row's axial chord. |
| `tip_clearance` | Rotor tip gap [m]. Zero is a shrouded rotor. |
| `hub_clearance` | Gap under a cantilevered stator, over a rotating drum [m]. |

Blade heights follow continuity at constant mean radius, and every row twists
as a free vortex from hub to tip.

A `SectionProfile` is a section's shape, in fractions of its chord:

| Setting | Meaning |
|---|---|
| `max_thickness`, `leading_edge_diameter`, `trailing_edge_thickness` | Thicknesses, over chord. |
| `trailing_edge_wedge_deg` | Suction- and pressure-side wedge angles at the trailing edge. |
| `distribution` | A `ThicknessDistribution`: suction- and pressure-side control-point offsets along the camber line, how far along it they reach, and how they are spaced. The default is front-loaded, as a reaction blade or nozzle vane is; `ThicknessDistribution.symmetric()` peaks at mid-chord. |
| `stagger_fraction` | Where the chord line sits between the inlet (0) and outlet (1) metal angles. Default 0.5, their mean. |
| `stagger_deg` | Stagger at mean radius, instead of `stagger_fraction`. |

`SectionProfile.impulse()` is a starting point for impulse rotor blades:
symmetric thickness and equally thin edges.

## The definition

`turbine.definition` is a `TurbineDefinition`: hardware only. It holds a flat
tuple of `BladeRow` (count, axial position and chord, annulus radii, metal
angles and stagger at mean radius, thicknesses, clearances, admission, and a
`twist` law), plus the `MeridionalPassage` through them. The passage records
which walls turn and which blade gaps exist:

- the hub turns under each rotor;
- the casing turns over a shrouded rotor;
- a rotor tip gap is recorded against the casing, and a cantilevered
  stator's gap against the drum under it.

A definition can also be written by hand, with no meanline behind it:

```python
from pyskyfire.turbine import BladeRow, TurbineDefinition

definition = TurbineDefinition(rows=(BladeRow(...), BladeRow(...)),
                               rotation="clockwise", outlet_length=0.004)
```

`definition.exit_annulus` is the passage's outlet: radii, area and wall
directions, for whatever joins it there.

## Built geometry and what is measured on it

`build_geometry(span_points=9, chord_points=129)` draws each row's sections
with pyTurbo at their span radii. It stacks them on their centroids and wraps
them onto the annulus. A rotor's tips stop `tip_clearance` short of the
casing. Each `RowGeometry` holds one blade's suction and pressure surfaces, in
cylindrical `(x, r, theta)`.

These are measured on the built sections rather than given, at every span:

- `throat`: the shortest line from the trailing edge to the neighbouring
  blade. `gauging_angle`, `acos(throat/pitch)`, is roughly the angle the flow
  leaves at. Set beside the outlet metal angle, it shows how far the drawn
  section turns the flow from what the triangles asked for;
- `narrowest` and `narrowest_at`: the narrowest gap anywhere between
  neighbours, and where along the axial chord it is. `inlet_width` is the
  width the flow enters through, `pitch*cos(inlet metal angle)`. `pinch` is
  the narrowest gap over the smaller of the inlet width and the throat: below
  one, the flow speeds up into the pinch and slows again after it, which
  costs loss. `narrowest_height` is the blade height at the pinch, for
  comparing areas rather than widths;
- `max_thickness`: the section's actual thickness.

Stagger decides most of the difference between the gauging and metal angles.
A row that turns the flow far from axial needs its chord leaned well towards
the outlet angle. With the mean-angle default, a vane turning 70° leaves its
flow more than 10° short. Raise `stagger_fraction` until the measured gauging
angles are where you want them. An impulse rotor takes its flow in as
steeply as it lets it out, so its inlet is as narrow as its throat. A
front-loaded thickness or a thick leading edge then narrows the passage
further, just inside the inlet; `SectionProfile.impulse()` avoids that.

`TurbineComponent` shows the passage walls and repeats each row's blade round
the wheel, with rotors and stators in different colours.

## Plots

`pyskyfire.viz` has four Plotly figures for a turbine, ready for
`Report.add_figure`:

- `plot_velocity_triangles(turbine)`: each stage's rotor inlet and exit
  triangles.
- `plot_turbine_meridional(geometry)`: the passage, the blade outlines and
  the wall stretches that turn with the shaft.
- `plot_turbine_cascade(geometry, span_fraction=0.5)`: every row's section at
  one span, unrolled.
- `plot_turbine_spanwise(geometry)`: metal and gauging angles, throat and
  pitch, hub to tip.

`examples/turbine/turbine_report.py` puts all of them and the 3D scene in one
report, for a single-stage impulse oxygen turbine with a shrouded rotor.
