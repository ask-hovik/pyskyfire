# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

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
