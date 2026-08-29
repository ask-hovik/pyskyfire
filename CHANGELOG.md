# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added
- Coupled regenerative-cooling solves now accept independent wall,
  heat-flux, and lumped coolant node grids through the `nodes` argument.
- Wall-temperature plots mark nodes that exceed the solver's scaled residual
  tolerance with red `x` symbols.
- `ThrustChamber.coolant_volume` reports the wetted coolant volume of all
  cooling circuits combined, and `CoolingCircuit.n_channels` gives a circuit's
  total physical channel count.

### Fixed
- `CoolingCircuit.compute_volume` counted circumferential leaves only, so a
  placement with `n_channels_per_leaf > 1` under-reported its coolant volume by
  that factor. Circuit and solver now share one channel count,
  `CoolingCircuit.n_channels`.
- The coolant march now integrates static pressure against the full momentum
  equation, `dp/dx = -friction - rho*u*du/dx`, instead of marching stagnation
  pressure against friction alone and recovering static pressure as
  `p0 - rho*u^2/2`. The old scheme counted heating-driven acceleration at half
  strength, because at constant area the acceleration loss is `G^2 d(1/rho)`,
  which is twice the change in dynamic head. On the RL10 validation case the
  predicted jacket static-pressure drop rises from 0.50 to 0.63 MPa against a
  reference range of 0.94 to 1.39 MPa.

### Changed
- `CoupledHeatExchangerPhysics.coolant_pressure_rate` is replaced by
  `coolant_friction_rate`, which returns only the friction gradient. The
  acceleration term is applied by the march station to station. The removed
  method also had a sign error in its static-pressure return value, which
  predicted a static-pressure *rise* through a contraction; that return value
  was unused by the solver.
- `CoupledHeatExchangerPhysics.bulk_velocity` returns the bulk density and
  single-channel velocity at a station.
- Visualization now chooses the resolution of hot-gas property plots.
  `PlotTransportProperty` takes exactly one of `results=` (reuse the axial
  stations of a solved run, the recommended default), `nodes=` (uniform grid
  over the contour), or `x=` (explicit coordinates). Repeated plots over the
  same grid reuse the transport object's station cache, so CEA runs once per
  station no matter how many properties are plotted.

### Removed
- The `npts` argument of every `Aerothermodynamics` sizing constructor, and the
  uniform `x_nodes` grid that `compute_aerothermodynamics` built when a contour
  was attached. Gas states are solved live at the requested coordinate, so the
  simulation no longer carries a visualization resolution. Pass the grid at
  plotting time instead. `compute_aerothermodynamics` also no longer accepts
  the deprecated, ignored `Nt` argument.

---

## [0.2.2] - 08-09-25
## Added
- Adding Cantera back in as a possible combustion backend model. CEA_Wrap will live on as an option when propellants not supported by cantera is needed. This will be used for all advanced implementations, including actual delivered isp estimates. 

---

## [0.2.1] - 08-09-25
## Added
- Various bug fixes

---

## [0.2.0] - 08-09-25
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
