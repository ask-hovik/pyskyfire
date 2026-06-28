# Start Here

Pyskyfire is a Python library for analysing liquid-propellant rocket engines.

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

The documentation follows the [Diátaxis](https://diataxis.fr/) framework. Each section serves a different purpose.

### Tutorials

Tutorials are guided, runnable introductions to a complete workflow.

* {doc}`Minimal Simulation <tutorials/minimal-simulation>`
  Build and analyse a regeneratively cooled thrust chamber from a defined engine design point.

### Howto

Howtos guide you through achieving a specific goal. For example hoe to use pyskyfire to find an optimum mixture ratio: {doc}`Mixture ratio optimisation <howto/mixture-ratio-optimisation>`  

### Explanations

Explanations describe the engineering models, code structure, and analysis methods used by Pyskyfire.

* {doc}`Capabilities <explanations/capabilities>`
  Overview of the analyses, outputs, and visualisation tools provided by Pyskyfire.

* {doc}`Package Structure <explanations/package-structure>`
  Overview of the repository, source packages, examples, validation cases, and their roles.

* {doc}`Regenerative Cooling <explanations/regenerative-cooling-explanation>`
  Explanation of the regenerative-cooling model, solver structure, heat-transfer paths, and governing quantities.

### Reference

The reference section is generated from the source code and documents the available Python interfaces.

| Package                                                    | Purpose                                                                                                                   |
| ---------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| {doc}`pyskyfire.common <autoapi/pyskyfire/common/index>`   | Shared objects including fluids, material models, result containers, stations, blocks, and engine-network infrastructure. |
| {doc}`pyskyfire.skycea <autoapi/pyskyfire/skycea/index>`   | NASA CEA-based combustion and aerothermodynamic calculations, nozzle utilities, and coolant transport properties.         |
| {doc}`pyskyfire.regen <autoapi/pyskyfire/regen/index>`     | Thrust-chamber geometry, cooling channels, wall models, regenerative cooling, and film cooling.                           |
| {doc}`pyskyfire.pump <autoapi/pyskyfire/pump/index>`       | Pump and turbopump analysis utilities.                                                                                    |
| {doc}`pyskyfire.turbine <autoapi/pyskyfire/turbine/index>` | Turbine models and utilities for engine-cycle calculations.                                                               |
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
```

```{toctree}
:caption: Explanations
:maxdepth: 1
:titlesonly:
:hidden:

explanations/capabilities
explanations/package-structure
explanations/regenerative-cooling-explanation
```

```{toctree}
:caption: Reference
:maxdepth: 1
:titlesonly:
:hidden:

autoapi/pyskyfire/common/index
autoapi/pyskyfire/regen/index
autoapi/pyskyfire/pump/index
autoapi/pyskyfire/turbine/index
autoapi/pyskyfire/skycea/index
autoapi/pyskyfire/viz/index

```
