# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added
- Modular 3D scenes with thrust-chamber, impeller, and explicitly shaped pipe
  components; SI rigid placement and visual inlet/outlet anchors; optional
  cooling-result scalar fields in desktop and offline HTML output. Cooling
  channel instances share geometry built once per resolved seed. See the
  3D scene guide and `examples/visualisation/scene.py`.
- `make_engine_3d` and `Engine3DViewer` are deprecated in favour of `Scene`
  and `ThrustChamberComponent`; examples and RL10 post-processing use the scene API.
- Explicit wall roles: `Wall(..., role="main")` is the default;
  `role="liner"` identifies combustion-side layers. Channel thickness sums
  all main layers, while liners contribute to hot-side offset and conduction.
  The RL10 silver throat insert is explicitly marked as a liner.
- Shared lazy cooling-channel section planes plus analytical and numerical
  layout resolvers, including validation for exactly patternable interleaves.
- An advanced-example Onshape JSON exporter and matching FeatureScript that
  loft unique channel seeds and circularly patterns their instances. The
  FeatureScript builds the true constant-width-rib section, whose flat rib
  faces make the channel widen with radius while the rib does not; the JSON
  schema is at version 6 and carries a `sectionType` per channel and the two
  resolved construction arcs per station.
- Optional continuous hot-wall profiles in the Onshape export, with circuit
  hand-off transitions and Bézier blends between unlike neighbours. The
  FeatureScript lofts the wall one exported curve at a time along each `run`
  of sections with constant profile topology, then knits the strips. Pairing
  whole profiles by proximity instead failed past roughly ninety curves per
  profile, and twisted the wall against the channels once a helix advanced
  more than half a channel pitch between stations. The two profiles either
  side of a hand-off share no curves, so the export cuts both onto one set of
  azimuths and they pair up as well.
- Design-input functions: geometry inputs now declare the quantities they need
  by parameter name, drawn from a tiered vocabulary of contour, layout and
  channel-section quantities. See the new *Design-Input Functions* explanation
  page for the slot table and the full glossary.

### Changed
- **Breaking:** `pyskyfire.skycea` is now `pyskyfire.thermo`, which also owns
  `Fluid` (formerly in `common`) and fluid-state resolution (formerly in
  `regen.coolant_state`). `CoolantTransport` becomes `FluidTransport`,
  `CoolantState` becomes `FluidStateResolver`, `BulkState` becomes `FluidState`,
  `CoolantCapability` becomes `FluidCapability`, and `GasState` becomes
  `CombustionGasState`. `TransportProperties`, `SaturationState`,
  `Aerothermodynamics`, and `NozzleLoss` keep their names under `thermo`.
  `probe_coolant` becomes `thermo.state.probe_fluid`; film-property policy moves
  to `regen.coolant_properties.film_properties(resolver, T_film, p, getters)`.
  No backwards aliases are provided. Old stored pickle files may require
  regenerating due to module and class-name changes; compatibility is not guaranteed.
- **Breaking:** `enable_fin` now belongs to `CoolingCircuit` (default `True`),
  not `ThrustChamber`. Each circuit can also set `rib_material`; omitting it
  retains the coolant-facing main-wall material for the fin calculation.
  The RL10 case sets fin enhancement on its two circuits individually.
- **Breaking:** `CrossSectionSquared` now takes `rib_width` instead of
  `blockage_ratio`, and models the rib as having the same thickness from the
  channel floor to the channel roof rather than as a wedge. Pass
  `rib_width=lambda pitch_mid: f * pitch_mid` to reproduce a blockage ratio of
  `f` exactly. The 3D viewer draws the area-equivalent wedge, since it cannot
  yet draw the true section.
- **Breaking:** design-input functions are now called with keyword arguments
  and must name what they need. `lambda x: ...` keeps working, a function that
  ignores position becomes `lambda: ...`, and `lambda *args: ...` is rejected
  with an explanatory error. This covers `Wall(thickness=...)`,
  `CoolingCircuit(channel_height=..., roughness=...,
  hot_gas_surface_area_multiplier=...)`, and the placement classes'
  `helix_angle` and `channel_angular_width`.
- `SurfacePlacement.helix_angle` and `InternalPlacement.channel_angular_width`
  are now methods that evaluate the design input, rather than stored callables.
  They, and `dtheta_dx`, take the circuit as an optional trailing argument so
  that a placement shared by several circuits evaluates against each in turn.
- The advanced example's nozzle circuits now use an area-preserving helix that
  starts at 20 deg, holding the coolant section constant instead of letting it
  swell with the nozzle radius.
- Renamed angular cooling-channel width APIs from `channel_width` to
  `channel_angular_width`; previously saved result files are migrated when loaded.
- Refactored the 3D channel visualizer to consume the shared layout and plane
  helpers; numerical redistribution remains its default.
- **Breaking:** a cooling circuit now carries three explicitly named paths in
  place of `centerlines`/`coolant_centerlines`. `reference_path` is the ideal
  contour carrying the nominal winding, `cross_section_base_point_path` is
  where the cross section is anchored (the coolant-facing floor for a squared
  section, the inner lobe centre for a rounded one), and
  `cooling_channel_centroid_path` is the coolant-area centroid that sets
  passage length and curvature. Circuits with the same helix angle now share
  one winding and cannot drift apart. Saved geometry using the old definitions
  is rejected on load rather than renamed; regenerate it.
- **Breaking:** `ChannelPlacement.compute_centerline_radius` and
  `compute_flow_centerline_radius` are gone. A placement owns the
  circumferential packing and the winding; where the section sits radially is
  resolved from the cross section and the wall stack.
- Design input `r_channel_base` is now the resolved cross-section base point
  rather than `r + t_wall`, so it accounts for the wall stack's projection onto
  a sloping contour and for a rounded section's solved lobe placement.
- The Onshape exporter rejects cross sections its FeatureScript cannot sketch
  instead of exporting arcs that would be mis-built.
- Removed the superseded drawing helpers in `pyskyfire.viz.engine_viz`
  (`get_section_scalars`, `rounded_lobe_arcs_yz`, `rounded_slice_arcs_uv`,
  `rounded_base_height` and their local helpers). The viewer and the exporter
  both build sections from `build_channel_construction`.

### Fixed
- Resolved channel geometry is built once per geometry instead of once per
  solver residual. A single-station section has to be interpolated out of the
  whole-circuit bundle, since it carries resolved fields such as the lobe
  sector radius and the winding rate that nothing can reconstruct from a
  radius alone. Assembling that bundle walks every station of every wall and
  evaluates each design-input function a station at a time, and the wall
  balance reached it on every residual it formed. The bundle is now cached
  against the path it was built from, single stations against the station they
  describe, and each section's boundary integral against the profile it
  measures; the profile setters and the construction steps invalidate all
  three. Results are bit-identical. The RL10 validation case goes from 30 s to
  4 s, and the advanced example's full cycle from 46 s to 36 s. The worst of it
  was an engine with a callable wall thickness, such as the RL10 silver
  insert: its thickness function was called 677,000 times per run, and
  two-thirds of those were for tier B quantities that a rounded section never
  reads.
- Squared channel ribs now point at the engine axis. Their faces used to be
  ruled by the wall normal swept along the channel, which showed up in a
  transverse CAD section as ribs leaning by up to 20 degrees on a sloping
  helical wall. They are ruled by the radial direction instead.
- A helical rib is no longer thinner than asked for. Both of its faces are
  helicoids, and one angular rate is a steeper helix further out, so the
  material between neighbouring channels closes in above the channel floor --
  about 4 percent over the channel height in the advanced example. The
  azimuthal offset is now sized so that the tightest point, at the channel
  roof, is exactly ``rib_width``. ``rib_width`` therefore means the minimum
  material through the solid, which is what a CAD measurement reports.
- Rounded channel lobes are now sized at the radius they actually end up at.
  Each tube is tangent to the sector that separates it from its neighbours and
  keeps that sector's share of the combustion-gas area, instead of being sized
  on the contour and then translated outward. Neighbouring tubes therefore meet
  with the intended amount of material between them in CAD. The RL10 numerical
  regression baseline moved by about one percent as a result; the experimental
  validation cases were unaffected.
- Corrected the rounded full-height conversion, which used
  `h_base*(sin(theta/2) + cos(theta/2))` where the geometry gives
  `h_base*(1 + sin(theta/2))`.
- Corrected squared-channel visualization to account for the ribs and to
  avoid counting hot-side wall layers twice when locating the wetted passage.

---

## [0.3.0] - 30-08-2026

### Added
- A new coupled cooling solver that balances hot-gas convection, conduction
  through multi-material walls, and coolant heat pickup. It supports separate
  wall, heat-flux, and coolant grids, reports convergence at every wall node,
  and returns a richer `RegenResult` for post-processing.
- Enthalpy-based coolant marching for pure fluids, including saturation-state,
  phase, vapour-quality, and boiling-region diagnostics. Mixtures and fluids
  without saturation data retain a temperature-based fallback.
- Standalone and regeneratively coupled subcritical film-cooling analysis based
  on the Grisson model, including liquid-film evaporation and dryout, gaseous
  film entrainment and radiation, and result plots and examples.
- A redesigned cooling-channel geometry system supporting straight, slanted,
  helical, and interleaved circuits; multiple channels per placement leaf;
  variable channel dimensions, and uneven number of channels in each interleaved circuit. 
- Interactive 3D views of the chamber and its actual cooling
  passages.
- Interactive engine-cycle network diagrams, mass-flow Sankey diagrams, and a
  much broader plotting suite for regenerative cooling, film cooling, wall
  layers, coolant phase, hot-gas fields, and engineering reference charts.
- Markdown blocks, responsive layouts, improved tables, and embedded
  interactive content in standalone HTML engineering reports.
- A full Sphinx documentation site with API reference, tutorials, design
  examples, theory explanations, and project contribution guidance.
- Automated tests and CI coverage for package imports, examples, solver
  physics, geometry, serialization, plotting, reports, and validation data.
- End-to-end minimal, advanced engine, mixture-ratio optimization, boiling,
  curvature-correction, and film-cooling examples.
- A comprehensively rebuilt RL10A-3-3A validation case covering engine sizing,
  cycle layout, regenerative cooling, channel geometry, published reference data, and a single consolidated report.
- Citation metadata in `CITATION.cff`, so releases are archived on Zenodo under
  a citable DOI.

### Changed
- Regenerative cooling was overhauled from the legacy heat-exchanger ODE into
  a station-marching coupled solution. That is, film cooling and regenerative cooling is available in the same engine. Coolant energy is advanced on enthalpy
  where possible, wall temperatures come from local nonlinear heat balances,
  and pressure loss follows the full friction-and-acceleration momentum balance. Performance was improved by utilising caching during solve. 
- Hot-side, coolant-side, curvature, channel-rib, and multi-layer wall models
  were revised and made consistent with the physical channel count used by the
  solver, geometry, volume calculations, and visualizations.
- `Aerothermodynamics` was rebuilt around the native `cea` package. Gas states
  are now evaluated on demand along the attached contour, cached by station,
  serializable without native solver handles, and available for equilibrium or
  frozen transport calculations.
- Thrust-chamber contours, circuit placement, cross-sections, and 3D geometry
  generation were substantially reworked for stable interpolation and a common
  representation of real channel centerlines and dimensions.
- Hot-gas plot resolution is now selected when plotting rather than when
  attaching a contour. Plots can use a solved result grid, a uniform node
  count, or explicit coordinates and reuse cached CEA states across properties.
- Visualization dependencies moved to the optional `pyskyfire[viz]` extra, so
  the core installation no longer requires the plotting and meshing stack.
- The supported runtime is now Python 3.12 or newer, and development,
  documentation, and test dependencies are organized into dedicated groups.

### Fixed
- Corrected coolant pressure-drop calculations to include acceleration from
  heating and area changes with the correct sign and magnitude.
- Corrected physical channel counting and coolant-volume calculations when a
  placement leaf represents more than one channel.
- Fixed geometry errors affecting ribs, rounded cross-sections,
  helical passages, contour splines and curvature, wall layers, and their 3D
  representations.
- Improved hot-gas property caching and low-temperature behavior, result
  serialization, report sizing and table layout, and solver diagnostics for
  unconverged wall nodes.

### Removed
- The legacy `regen.solver` implementation and its `steady_heating_analysis`
  entry point; use `coupled_steady_heating_analysis` and the new result model.
- Deprecated `CoolingCircuitGroup` and `WallGroup` wrappers; cooling circuits
  are passed directly to `ThrustChamber`, and wall layers are defined per
  circuit.
- The Cantera combustion path and obsolete combustion/nozzle-solver modules.
  Chemical-equilibrium and hot-gas calculations now use `Aerothermodynamics`
  backed by `cea`.
- The `npts` sizing-constructor argument, stored aerothermodynamic plotting
  grid, and deprecated `Nt` argument. Resolution now belongs to the requested
  solve or visualization.

---

## [0.2.1] - 08-09-2025
## Added
- Various bug fixes

---

## [0.2.0] - 08-09-2025
## Added
- New combustion properties driver CEA_Wrap, replaces RocketCEA and Cantera. The old interface to RocketCEA and Cantera, CombustionTransport, is replaced with Aerothermodynamics. Leverages the large amount of work that has gone into the NASA CEA database. This means that all of NASA CEAs propellants are now supported! Legacy CombustionTransport depricated. 
- New method of hot gas heat transfer estimation in the regen solver based on the new Aerothermodynamics class. 
- Revamped cross sections to use improved class structure.
- Fixed an issue in the CrossSectionSquared class where the rib width was impossible to control. 
- Added a new fluids interface that allows for fuel and oxidizer blends. An example would be an ethanol-water mix.
- Rewrote entire visualisation pipeline. HTML Reports are now the predominant way of viewing results. 
- Added minimal example
- Rewrote RL10 validation case

---

## [0.1.2] - 23-08-2025
### Added
- Packaging and distribution support: `pyskyfire` is now available via pip and PyPI.
- An updated `README.md` explaning the package.
- New logo, which is currently just "pyskyfire" written in the orbitron open-source font. 
- This `CHANGELOG.md` file to track notable changes.

---

## [0.1.0] - 29-06-2025
### Added
- First public release of **pyskyfire**.
- Rocket engine **multi-pass regenerative cooling solver**.
- **Full cycle solver** for expander-cycle engines and related architectures.
- **Visualisation tools** for engine performance and cooling analysis.
