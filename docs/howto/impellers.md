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
    blade_count=6, beta2_deg=25., blade_thickness=0.001,
    phi1=0.20, incidence_deg=2.,
    volumetric_efficiency=0.97, hydraulic_efficiency=0.85,
    mechanical_efficiency=0.98, rho=1000.,
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
`Q` is delivered volume flow and `H` is delivered head. Alternatively use
`Impeller.from_pressure_rise(mdot=..., rho=..., dp=..., n=..., **design_inputs)`.

The outlet diameter and width use the existing preliminary pressure-coefficient
and width correlations. Eye sizing retains `phi1 = Cm1 / Ueye`, solving annular
continuity with a fixed wetted hub diameter. `d1` can override the computed eye
diameter. The shaft diameter and wetted hub diameter are separate inputs; no
shaft stress sizing is performed.

`volumetric_efficiency` sets the internal flow to `Q / volumetric_efficiency`.
Its default of 1 means no leakage. Hydraulic and mechanical efficiencies are
optional supplied values. With both and `rho`, `shaft_power` is
`rho * g * Q * H / (eta_v * eta_h * eta_m)`. Otherwise it is `None`.

These efficiencies do not constitute a performance solver. The outlet velocity
triangle assumes blade-congruent flow without slip. Its Euler head and the
hydraulic-efficiency-adjusted head are diagnostics in `velocity_triangles`;
they are not forced to equal the requested head. There is currently no head
closure iteration, blade blockage correction, cavitation assessment, or loss
prediction. The generated geometry is a preliminary design for inspection.

## Design controls and conventions

- `blade_count`, `beta2_deg`, `blade_thickness`, `shaft_diameter` and
  `inlet_hub_diameter` are required design inputs.
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
must match the sized inlet radii, outlet radius and outlet width. Pair stations
across the passage; the shared station coordinate is normalized hub arc length.

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
