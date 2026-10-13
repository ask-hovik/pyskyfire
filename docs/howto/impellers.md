# Create an impeller scene and report

Install the geometry and visualization extras:

```sh
pip install 'pyskyfire[impeller,viz]'
```

`Impeller` performs preliminary sizing and resolves physical blade inputs.
`build_geometry()` calls pyTurbo Aero to construct the pressure and suction
surfaces. `ImpellerComponent` places these surfaces in the regular scene API:

```python
from pyskyfire.pump import Impeller
from pyskyfire.viz import ImpellerComponent, Report, Scene

impeller = Impeller(
    Q=0.01, H=100., n=10_000,
    shaft_diameter=0.020, inlet_hub_diameter=0.030,
    blade_count=6, blade_thickness=0.001, hydraulic_efficiency=0.90,
    phi1=0.20, incidence_deg=2., rho=1000.,
)
impeller.build_geometry()

scene = Scene(title="Impeller", show_axes=True)
scene.add(ImpellerComponent(impeller), name="impeller")
scene.save_html("impeller-scene.html")

report = Report("Impeller design")
report.add_tab("Geometry").add_iframe(scene.data_url)
report.save_html("impeller-report.html")
```

The complete runnable example is `examples/impeller/impeller_report.py`.
It adds a design table, meridional plot, and blade-angle table to the report.
It is designed at low specific speed with four full blades and one splitter
per passage. Run it with `--screenshot` for a PNG preview. HTML generation
does not need a display.

## Duty, sizing and supplied efficiencies

Use SI units with rotational speed in rpm and blade angles in degrees.
`Q` is delivered volume flow and `H` is the total head at the impeller
outlet: the head the pump delivers plus whatever the volute or diffuser after
it loses. Alternatively use
`Impeller.from_pressure_rise(mdot=..., rho=..., dp=..., n=..., **design_inputs)`,
or size for a `PumpDuty` with `Impeller.from_duty(duty, **design_inputs)`.

A `PumpDuty` states what the pump is asked to do, separately from how it is
built. It can be written the way pumps are specified or the way an engine
presents them, and stores the same thing either way:

```python
from pyskyfire.pump import PumpDuty
from pyskyfire.thermo import Incompressible

PumpDuty.from_head(0.01, 100., 10_000, rho=1000.)          # Q, H, n
PumpDuty.from_pressure_rise(mdot=10., rho=1000., dp=9.8e5, n=10_000)
PumpDuty.from_states(fluid="Methane", mdot=22.8, p_in=5e5, T_in=135.3,
                     p_out=230e5, n=40_000)                  # real fluid
PumpDuty.from_stations(tank, pump_out, fluid="Methane", n=40_000)
```

`fluid` is a CoolProp name, a `thermo.Fluid`, or any `thermo.EquationOfState`,
such as `Incompressible(rho)`. `H` is the isentropic head: the enthalpy an
ideal pump adds, over `g`. At constant density that is `Δp/(rho*g)`. For a
real liquid it is a little less, because the liquid compresses on the way
through; for methane pumped from 5 to 230 bar the difference is about 2 %. The
duty also reports the isentropic outlet flow `Q_out` and the specific speed.

Eye sizing retains `phi1 = Cm1 / Ueye`, solving annular continuity with a
fixed wetted hub diameter. `d1` can override the computed eye diameter. The
shaft diameter and wetted hub diameter are separate inputs; no shaft stress
sizing is performed. `b2` defaults to Gülich's mean outlet-width curve,
Eq. (7.1).

### Slip and the head closure

The flow does not leave along the blades. Gülich's Table 3.2 predicts the
swirl it does leave with: a slip factor from the blade outlet angle, the blade
count, the blade blockage and the inlet diameter, a modified form of Wiesner's
correlation. `outlet_slip()` evaluates it on its own. The outlet triangles
carry that swirl, and `beta_flow_deg` is the flow angle it implies, below the
blade angle `beta_metal_deg`.

The design is closed on it. Exactly one quantity is solved so the predicted
swirl is the swirl the head needs, `c_u2 = g*H/(eta_h*u2) + u1*c_u1/u2`
(Gülich Eq. T3.2.8):

- By default the blade outlet angle is solved, which is the order of Gülich's
  own procedure. The outlet diameter comes from one of `psi`, `d2`, `u2` or a
  supplied `meridional`. With none of them it comes from Gülich's head
  coefficient curve, Eq. (3.26), at `fT` (default 1.0; he finds 1.0 describes
  multistage pumps and 1.1 low-head single-stage volute pumps).
- Give `beta2_deg` and the outlet diameter is solved instead.
- Giving both over-specifies the wheel and is refused.

Eq. (3.26) describes whole pumps: their delivered head at a stage efficiency.
`H` and `hydraulic_efficiency` here are the impeller's, so the same wheel sits
higher in the band by roughly the ratio of the two efficiencies. That is the
reason to move `fT` or give `psi` outright.

At low specific speed the swirl a wheel can deliver has a broad, flat peak
against the blade angle. A head above it cannot be had at any angle, and the
error says what head coefficient can. The solved angle is always the one below
the peak.

`blade_count` is an input. Gülich rejects empirical blade-count formulas as
likely to violate the criteria that actually set it. Only full blades count
toward slip; splitters are not in the correlation.

`impeller.diagnostics` reports what was solved, the resolved and correlated
head coefficients, the mean-streamline slip prediction, the deviation angle
(Gülich suggests at most 10 to 14 degrees) and the deceleration ratio
`w2/w1` at the shroud (de Haller: at least 0.7). None of them is enforced.

### Efficiencies

`hydraulic_efficiency` is the impeller's own, Gülich's eta_h,La: the losses
between its inlet and outlet. Each other loss is charged by the component it
happens in, never lumped into one pump efficiency:

- the volute or diffuser, by its own block, from the head it is handed;
- disk friction on the wheel's outer faces, by `pump.DiskFriction` (Gülich
  Table 3.6);
- bearings and shaft seals, by the pump block's `mechanical_efficiency`;
- leakage, by `volumetric_efficiency`, which sets the internal flow to
  `Q / volumetric_efficiency`. Its default of 1 means no leakage.

`blade_power` is `rho * g * Q_internal * H / eta_h`, the work the blades do,
or `None` without `rho`. `network.PumpBlock` combines the parts: its outlet is
the impeller outlet, its shaft power is blade power plus disk friction over
the mechanical efficiency, and the first two heat the fluid.

For a typical efficiency, Gülich's Fig. 3.33 splits the losses of pumps near
`nq = 15` into about 2 to 3.5 percent of the head in the impeller against 13 to
18 percent in the collector, which puts eta_h,La near 0.97 where the whole
stage is near 0.85.

The efficiencies are supplied assumptions, not predictions. There is no
cavitation assessment or impeller loss prediction. The generated geometry is a
preliminary design for inspection.

## Design controls and conventions

- `blade_count`, `blade_thickness`, `hydraulic_efficiency`, `shaft_diameter`
  and `inlet_hub_diameter` are required design inputs. `beta2_deg` is solved
  unless given.
- `beta1_deg=None` derives hub, midspan and tip inlet metal angles from the
  unblocked annular velocity triangles, with `incidence_deg` added. An explicit
  `beta1_deg` overrides this calculation. `inlet_prewhirl` specifies absolute
  circumferential inlet velocity in m/s, positive with rotation.
- Metal angles are measured from the backward circumferential direction toward
  the meridional direction: 90° is meridional. This implementation supports
  backward-swept blades. Incidence means metal angle minus relative flow angle.
- Span runs from hub (`span_fraction=0`) to blade tip (`span_fraction=1`).
  The three defining profiles are interpolated by pyTurbo across that span.
- `wrap_deg` is total circumferential sweep from leading to trailing edge,
  positive backward for rotation about +x. If omitted, it is estimated by
  integrating an interpolated metal-angle law along each meridional section.
  This integral sets the total wrap; pyTurbo constructs its own camber curve
  satisfying the endpoint controls.
- Sections are stacked on the trailing edge: it lies in one axial plane at
  `theta = 0`, and the leading edge carries the spanwise wrap difference as
  sweep. Stacking rotates each section rigidly about the shaft, leaving metal
  angles, wrap and thickness unchanged.
- `BladeSettings.start` slants the leading edge when given a `(hub, tip)` pair
  of stations. Slant it: a hub streamline is much longer in `m'` than a tip
  streamline, because `dm' = dm/r` and the hub runs close to the axis, so a
  leading edge on one station makes the hub section wrap far further than the
  tip section. That difference ends up as leading-edge sweep. Starting the
  blade further downstream on the hub, as a blade development is normally
  drawn, evens the wrap out. In the example, `start=(0.2, 0.01)` takes the
  spread from 65° to 6°.
- The inlet velocity triangle and the inlet metal angles are evaluated on the
  leading edge, not at the inlet plane, with `cm` from unblocked continuity
  through the leading edge. The component of the through-flow crossing that
  line is what counts, so a slanted edge passes less than the area it sweeps.
  Eye sizing still uses the inlet plane and `phi1`.
- `tip_clearance` is a distance in metres along the paired hub-to-shroud span.
  It trims the blade-tip boundary while preserving the displayed flow shroud.

`blade_thickness` is the total physical section-normal **body** thickness in
metres. Supply a float or a function requesting any subset of the named scalar
arguments `chord_fraction`, `span_fraction`, `x`, `r`:

```python
def thickness(chord_fraction, span_fraction):
    return 0.0015 * (1 - 0.30*span_fraction) * (1 - 0.20*chord_fraction)
```

`chord_fraction` is normalized conformal meridional distance, with `dm′ = dm/r`,
from main-blade leading edge (0) to trailing edge (1). `x` and `r` are in metres.
Angle, incidence, prewhirl and wrap inputs also accept functions, with arguments
`span_fraction`, `x`, `r`. Functions return a finite scalar and should be smooth.
Angle and wrap functions are sampled at the three defining spans.

Thickness is fitted to pyTurbo's surface controls. Leading and trailing edges
close, so the thickness tapers there. Body thickness is checked at the defining
profiles; fits exceeding 20% relative error away from the edge closures are
rejected. `geometry.thickness_max_relative_error` reports this fit error. A
sharp or rapidly varying function is unsuitable for the current fitter.

Optional `BladeSettings(start=0.05, end=0.98, trailing_edge="rounded",
angle_law="cosine")` controls blade placement and wrap estimation. Start and
end are fractions of the shared meridional station coordinate. `start` also
accepts a `(hub, tip)` pair, which slants the leading edge; `start_hub` and
`start_tip` always report both. A section can only begin on a span line, so a
slant is a statement about which hub and shroud stations pair, and the passage
is re-paired to carry that line: the hub keeps its stations and the shroud is
resampled, then relaxed back to the original pairing downstream. A supplied
`meridional` passage is re-paired the same way. Pairings that would push the
shroud below the hub are rejected. Trailing edge
can be `"rounded"` or `"cut"`; independent edge-radius inputs are not exposed
because the supported upstream version does not reliably honor them.

Optional `SplitterSettings(start=0.5, pitch_fraction=0.5, thickness_scale=1.0)`
adds one splitter per main-blade passage. It follows the downstream main-blade
camber, starting at the specified main-blade chord fraction. `pitch_fraction`
places it circumferentially between main blades. Its thickness function uses
the main-blade chord coordinate, multiplied by `thickness_scale`.

## Meridional passage and generated representation

The default passage uses smooth quarter-ellipse hub and shroud curves with an
axial inlet and radial outlet. `axial_length` sets the shroud's axial extent;
otherwise a preliminary dimension-based default is used. Supply
`meridional=MeridionalContours(hub=hub_xr, shroud=shroud_xr)` for custom paired
`(N, 2)` arrays of `(x, r)` stations. Both arrays must have the same length,
positive passage width and monotone axial/radial coordinates. Their endpoints
must match the sized inlet radii and, if `b2` is given, the outlet width; the
outlet radius is taken from it and the blade angle solved. Pair stations
across the passage; the shared station coordinate is normalized hub arc length.

`MeridionalContours` is the radial-impeller form of the shared
`pyskyfire.common.MeridionalPassage`, so a passage also reports its `area()`
and `hydraulic_diameter()` at each station and can record which wall stretches
rotate (`WallMotion`) and where blades stop short of a wall (`Clearance`).
The impeller does not fill those in yet.

`impeller.dimensions`, `impeller.velocity_triangles` and
`impeller.blade_definition` are available before importing pyTurbo. Geometry
generation is explicit and can be repeated at another resolution:

```python
geometry = impeller.build_geometry(span_points=33, chord_points=129)
geometry.export_npz("impeller-geometry.npz")
```

The geometry contains the hub and shroud curves plus one main blade and,
optionally, one splitter. Each blade has pressure and suction arrays shaped
`(span_points, chord_points, 3)` in cylindrical `(x, r, theta)` coordinates;
theta is in radians. Cartesian conversion uses `(x, r*cos(theta), r*sin(theta))`.
The visualizer caps the blade ends and instances the blade mesh around the
axis. Hub and shroud are flow boundaries, not a solid wheel or casing.

The adapter targets `pyturbo-aero==1.3.8` and translates its angle, wrap,
thickness and clearance conventions. This dependency is optional and imported
only by `build_geometry()`. Native blade surface generation and the two old
experimental impeller modules have been removed. Existing users should import
`Impeller` from `pyskyfire.pump`, supply the required design inputs, and call
`build_geometry()` before adding an impeller to a scene. Create a new `Impeller`
when changing design inputs.
