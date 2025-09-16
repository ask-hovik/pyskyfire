# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

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
