# Start Here

Pyskyfire is a Python library for analysing liquid-propellant rocket engines.

**Project links:** [GitHub](https://github.com/ask-hovik/pyskyfire) · [PyPI](https://pypi.org/project/pyskyfire/)

It supports:

* combustion-performance and nozzle-flow calculations using NASA CEA;
* thrust-chamber and nozzle contour generation;
* regenerative-cooling analysis of chamber walls and coolant channels;
* cooling-channel geometry, wall-material, and coolant-property models;
* pump, turbine, and engine-network utilities;
* plots, interactive 3D geometry views, and standalone HTML reports.

Pyskyfire is intended for preliminary engine design, thermal analysis, and engine-cycle studies.

## Start with the minimal simulation

{doc}`Minimal Simulation <tutorials/minimal-simulation>` is the recommended introduction to Pyskyfire.

The tutorial builds and analyses a 5 kN nitrous-oxide/ethanol rocket engine. It shows how to:

1. define an engine design point;
2. calculate combustion-gas properties;
3. generate a thrust-chamber contour;
4. define walls and cooling channels;
5. run a regenerative-cooling simulation; and
6. generate plots, a 3D model, and an HTML report.

The complete source code is available in `examples/minimal/minimal_sim.py`. Use it as a starting point for new engine cases.

## Documentation structure

The documentation is organised around the [Diátaxis](https://diataxis.fr/)
framework, with Validation added as a first-class section. Each section serves
a different purpose.

### Tutorials

Tutorials are guided, runnable introductions to a complete workflow.

* {doc}`Minimal Simulation <tutorials/minimal-simulation>`
  Build and analyse a regeneratively cooled thrust chamber from a defined engine design point.

### Howto

Howtos guide you through achieving a specific goal. For example hoe to use pyskyfire to find an optimum mixture ratio: {doc}`Mixture ratio optimisation <howto/mixture-ratio-optimisation>`  

* {doc}`Define and solve a pipe <howto/pipes>`
* {doc}`Define and solve an annular duct <howto/annular-ducts>`
* {doc}`Size and solve orifices <howto/orifices>`
* {doc}`Design and solve an injector <howto/injectors>`
* {doc}`Design and solve a volute <howto/volutes>`
* {doc}`Design a vaned diffuser <howto/diffusers>`
* {doc}`Design an axial turbine <howto/turbines>`
* {doc}`Hold the coolant area constant with a helix <howto/area-preserving-helix>`
  Construct explicit pipe stations and connections, solve pressure loss independently, or use the same pipe in an engine network.

### Validation

Validation reports compare Pyskyfire analyses with published engine data.

* {doc}`RL10A-3-3A <validation/rl10a-3-3a>`
  Reconstruction and validation of the RL10A-3-3A thrust chamber,
  regenerative cooling system, and expander cycle.

### Explanations

Explanations describe the engineering models, code structure, and analysis methods used by Pyskyfire.

* {doc}`Capabilities <explanations/capabilities>`
  Overview of the analyses, outputs, and visualisation tools provided by Pyskyfire.

* {doc}`Package Structure <explanations/package-structure>`
  Overview of the repository, source packages, examples, validation cases, and their roles.

* {doc}`Design-Input Functions <explanations/design-input-functions>`
  Explanation of the quantities a geometry input function may depend on, which are available where, and why.

* {doc}`Regenerative Cooling <explanations/regenerative-cooling-explanation>`
  Explanation of the regenerative-cooling model, solver structure, heat-transfer paths, and governing quantities.

* {doc}`Specific Impulse and Thrust Coefficient <explanations/specific-impulse>`
  Explanation of the vacuum, optimum-expansion, ambient, and sea-level performance figures, and the conditions under which each applies.

* {doc}`The Moody Diagram <explanations/moody-diagram>`
  Explanation of the friction-factor model used for cooling-channel pressure loss.

* {doc}`Rao Nozzle Contour Angles <explanations/rao-nozzle>`
  Explanation of the angle data used to construct Rao-style bell nozzles.

### Reference

The reference section is generated from the source code and documents the available Python interfaces.

| Package                                                    | Purpose                                                                                                                   |
| ---------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| {doc}`pyskyfire.common <autoapi/pyskyfire/common/index>`   | Shared building blocks including solid-material models, result containers, friction correlations, and blade and section geometry. |
| {doc}`pyskyfire.thermo <autoapi/pyskyfire/thermo/index>`   | Fluid definitions, state resolution, transport properties, NASA CEA-based combustion and aerothermodynamics, and nozzle losses.         |
| {doc}`pyskyfire.regen <autoapi/pyskyfire/regen/index>`     | Thrust-chamber geometry, cooling channels, wall models, regenerative cooling, and film cooling.                           |
| {doc}`pyskyfire.ducts <autoapi/pyskyfire/ducts/index>`     | Shared duct friction, explicit pipe and annular-duct geometry, pressure-loss solving, and pipe-path generation.          |
| {doc}`pyskyfire.heat_exchanger <autoapi/pyskyfire/heat_exchanger/index>` | Path-based heat exchangers with prescribed gas or coupled fluid circuits. |
| {doc}`pyskyfire.injector <autoapi/pyskyfire/injector/index>` | Injector element sets and the flow from each feed through them. |
| {doc}`pyskyfire.pump <autoapi/pyskyfire/pump/index>`       | Pump and turbopump analysis utilities.                                                                                    |
| {doc}`pyskyfire.turbine <autoapi/pyskyfire/turbine/index>` | Turbine models and utilities for engine-cycle calculations.                                                               |
| {doc}`pyskyfire.network <autoapi/pyskyfire/network/index>` | Engine-network stations and blocks that connect components into a solvable cycle. |
| {doc}`pyskyfire.viz <autoapi/pyskyfire/viz/index>`         | Plotting, reporting, pressure-temperature diagrams, and 3D engine visualisation.                                          |

Use the reference when you need the exact inputs, outputs, methods, or inheritance structure of a class or function.

## Documentation and source code

Tutorials use runnable scripts from the repository. The minimal simulation tutorial includes selected sections of `examples/minimal/minimal_sim.py`, so the documented workflow matches the executed example.

The API reference is generated from docstrings in `src/pyskyfire/`.


```{toctree}
:maxdepth: 1
:titlesonly:
:hidden:

self
```

```{toctree}
:caption: Tutorials
:maxdepth: 1
:titlesonly:
:hidden:

tutorials/minimal-simulation
tutorials/advanced-simulation
```

```{toctree}
:caption: How-to Guides
:maxdepth: 1
:titlesonly:
:hidden:

howto/mixture-ratio-optimisation
howto/pipes
howto/annular-ducts
howto/orifices
howto/injectors
howto/heat-exchangers
howto/volutes
howto/cross-sections
howto/area-preserving-helix
howto/impellers
howto/diffusers
howto/turbines
howto/3d-scenes
```

```{toctree}
:caption: Validation
:maxdepth: 1
:titlesonly:
:hidden:

validation/rl10a-3-3a
```

```{toctree}
:caption: Explanations
:maxdepth: 1
:titlesonly:
:hidden:

explanations/capabilities
explanations/package-structure
explanations/design-input-functions
explanations/channel-construction
explanations/regenerative-cooling-explanation
explanations/specific-impulse
explanations/moody-diagram
explanations/rao-nozzle
```

```{toctree}
:caption: Reference
:maxdepth: 1
:titlesonly:
:hidden:

autoapi/pyskyfire/common/index
autoapi/pyskyfire/ducts/index
autoapi/pyskyfire/heat_exchanger/index
autoapi/pyskyfire/injector/index
autoapi/pyskyfire/regen/index
autoapi/pyskyfire/pump/index
autoapi/pyskyfire/turbine/index
autoapi/pyskyfire/thermo/index
autoapi/pyskyfire/network/index
autoapi/pyskyfire/viz/index

```
