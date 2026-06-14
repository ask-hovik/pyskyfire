# Package Structure

This page explains how the Pyskyfire repository is organised. It is written for users who are comfortable with engineering and Python scripts, but who may not have spent much time inside a Python package before.

The goal of this page is to give you a mental map of where things live, what each part is responsible for, and how the different parts of the project fit together.

## 1. The repository as a whole

When you open the Pyskyfire repository, you are looking at more than the Python code itself. A software repository usually contains:

- the source code,
- examples showing how the code is used,
- validation cases showing whether the code reproduces known results,
- documentation,
- configuration files telling Python and Git how to treat the project,
- images and generated content

A simplified view of the root directory is:

```text
pyskyfire/
├── .gitattributes
├── .gitignore
├── LICENSE
├── README.md
├── pyproject.toml
├── docs/
├── examples/
├── images/
├── src/
└── validation/
```

The most important folder for the actual library code is `src/`. The other folders explain, demonstrate, test, package, or present that code.

| File or folder | What it is for |
|---|---|
| `.gitignore` | Tells Git which files should not be tracked. In practice, this keeps local caches, virtual environments, build products, generated outputs, private data, and temporary files out of the repository history. |
| `.gitattributes` | Gives GitHub extra instructions about how to treat certain files. In Pyskyfire it marks generated HTML reports in `examples/` and `validation/` as generated, so they do not dominate GitHub's language statistics. |
| `LICENSE` | States the legal terms for using and modifying the project. Pyskyfire uses the MIT License, which is permissive, but also makes clear that the software is provided without warranty. |
| `README.md` | The front page of the GitHub repository. It introduces what Pyskyfire is, shows the main capabilities, gives a minimal installation command, and points users toward examples and validation cases. |
| `pyproject.toml` | The main configuration file for the Python package. It defines the package name, version, dependencies, build system, optional documentation dependencies, and tells the build system that the importable package lives in `src/pyskyfire`. |
| `docs/` | Contains the documentation website. This includes tutorials, how-to guides, explanations such as this page, and generated API reference pages. |
| `examples/` | Contains scripts that demonstrate how to use the package. These are primarily teaching material: they show workflows and interfaces that users can copy and adapt. |
| `images/` | Stores figures used in the README and documentation, such as plots, screenshots, and rendered engine or cooling-channel images. |
| `src/` | Contains the actual installable Python package. When a user writes `import pyskyfire`, Python imports code from `src/pyskyfire`. |
| `validation/` | Contains cases that compare Pyskyfire against reference data, published results, known solutions, or other tools. These are mainly about credibility rather than teaching the interface. |


## 2. The `src/pyskyfire/` package

The `src/pyskyfire/` directory contains the code that becomes the `pyskyfire` Python package after installation.

The package content is as follows:

```text
src/pyskyfire/
├── __init__.py
├── common/
├── regen/
├── skycea/
├── pump/
├── turbine/
└── viz/
```

Each subfolder is a *subpackage*. A subpackage is a folder containing related Python modules. A module is a single `.py` file.

For example:

```text
src/pyskyfire/regen/solver.py
```

is a module, and it belongs to the subpackage:

```python
pyskyfire.regen
```

### Top-level `__init__.py`

The file `src/pyskyfire/__init__.py` marks `pyskyfire` as a Python package and exposes the main subpackages.

### `common`

`common` contains objects and helper classes that are used by several other parts of Pyskyfire. It is where broadly useful concepts live: fluids, solids, results, engine stations, and engine-cycle blocks.

| Module | Role |
|---|---|
| `fluids.py` | Defines a `Fluid` helper for propellant and coolant mixtures. It stores component names, composition fractions, and whether the mixture is given on a mass or mole basis. |
| `solids.py` | Defines material-property models for solids. These can represent constant, polynomial, tabulated, or piecewise temperature-dependent properties such as thermal conductivity. |
| `results.py` | Defines a dictionary-like `Results` container with dot access and save/load functionality. This is useful for storing simulation outputs in a convenient way. |
| `engine_network.py` | Defines `Station` objects and the `EngineNetwork` fixed-point framework. This is used to connect pumps, turbines, cooling circuits, and other components into a larger engine-cycle calculation. |
| `blocks.py` | Defines the block interface used by the engine network. A block reads stations and/or scalar signals, performs a calculation, and returns updated stations and signals. |
| `constants.py` | Stores shared constants that may be used across modules. Keeping constants in one place avoids repeating the same numbers throughout the code. |
| `__init__.py` | Re-exports selected objects from the `common` subpackage so they are easier to import. |

### `regen`

`regen` contains the thrust-chamber and regenerative-cooling model. This is one of the central subpackages in Pyskyfire.

A useful mental model is:

```text
geometry + materials + coolant + hot-gas properties
                  ↓
          regenerative cooling solver
                  ↓
wall temperatures, coolant temperatures, pressure losses, heat fluxes
```

| Module | Role |
|---|---|
| `contour.py` | Defines the hot-gas contour of the chamber and nozzle. It provides local radius, area, throat geometry, expansion ratio, contraction ratio, and wall slope. |
| `cross_section.py` | Defines cooling-channel cross-section models. These provide flow area, hydraulic diameter, thermal perimeter, coolant wetted perimeter, and thermal resistance information. |
| `channel_height.py` | Contains logic for channel-height distributions along the chamber. This lets cooling-channel geometry vary with axial position. |
| `channel_placement.py` | Describes how cooling channels are placed around or inside the chamber wall. This includes channel counts and placement logic that affect per-channel mass flow. |
| `thrust_chamber.py` | Ties the physical thrust-chamber model together. It contains wall layers, cooling circuits, film-cooling inputs, and the chamber object used by the solver. |
| `physics.py` | Contains local engineering correlations used by the solver, such as Bartz-style hot-gas heat transfer, Colburn coolant-side heat transfer, Reynolds number, Darcy friction factor, and adiabatic wall temperature. |
| `solver.py` | Runs the regenerative-cooling calculation. It marches along the cooling circuit, solves local wall heat balances, and updates coolant temperature and pressure. |
| `film_solver.py` | Contains film-cooling model functionality. Film cooling is separate from regenerative cooling because coolant is injected into the hot-gas side rather than only flowing inside wall channels. |
| `film_solver_2.py` | Contains a newer or alternative film-cooling implementation. New users should normally start from documented examples or validation cases before relying on experimental modules. |
| `__init__.py` | Re-exports selected objects from the `regen` subpackage so common classes and functions are easier to import. |

### `skycea`

`skycea` handles chemical-equilibrium and transport-property calculations. It connects engine geometry and propellant choices to the gas and coolant properties needed by the rest of the package.

| Module | Role |
|---|---|
| `aerothermodynamics.py` | Computes and stores hot-gas properties along the chamber/nozzle contour using CEA-based calculations. The regenerative-cooling solver uses this to get gas temperature, pressure, Mach number, enthalpy, viscosity, thermal conductivity, Prandtl number, and other properties at each axial location. |
| `coolant_transport.py` | Defines transport-property containers for coolants and other fluids. Properties such as density, heat capacity, viscosity, conductivity, and Prandtl number can be constants or functions of temperature and pressure. |
| `nozzle_solver.py` | Contains nozzle and CEA-based performance utilities. It supports the broader role of connecting propellants, chamber conditions, nozzle geometry, and CEA outputs. |
| `data/` | Contains data files needed by the CEA wrapper, such as transport-property library data. |
| `__init__.py` | Re-exports selected skycea functionality for easier imports. |

### `pump` and `turbine`

The `pump` and `turbine` subpackages contain tools for turbopumps. They are both under development, and currently the pyskyfire package works fine without them. In the future, these packages will generate pump, turbine and manifold geometries, and expose other turbopump-relevant properties such as eigenfrequencies of the turbopump assembly. 

### `viz`: plotting, reporting, and geometry visualisation

`viz` turns simulation data into plots, reports, and geometry views. It is where the package moves from numerical results to things a user can inspect visually.

| Module | Role |
|---|---|
| `core.py` | Defines common Plotly helper functionality used by several plotting classes. It provides shared behaviour for configuring, showing, and saving figures. |
| `plot_regen.py` | Plotting tools for regenerative-cooling results, such as wall temperature, coolant temperature, pressure, heat flux, and heat-transfer coefficients. |
| `plot_film_cooling.py` | Plotting tools for film-cooling results. |
| `plot_skycea.py` | Plotting tools for aerothermodynamic and CEA-derived results. |
| `plot_common.py` | Shared plotting utilities used by several plotting modules. |
| `engine_viz.py` | 3D engine and cooling-channel visualisation tools. This includes geometry-oriented views of cooling circuits and engine structure. |
| `impeller_viz.py` | Visualisation tools for pump and impeller geometry. |
| `report.py` | Report-generation helpers for collecting plots, tables, and simulation summaries. |
| `lookup_tables.py` | Lookup-table support used by plotting or reporting workflows. |
| `__init__.py` | Re-exports selected visualisation functions and classes for easier imports. |

## 3. How the subpackages fit together

The subpackages each have a responsibility, but in order to run a regenerative cooling analysis or a full cycle analysis, the packages work together. A typical regenerative-cooling analysis might work like this:

1. `common` defines shared concepts such as fluids, solids, stations, blocks, and results.
2. `skycea` provides gas and coolant property models.
3. `regen` defines the thrust chamber, cooling channels, wall layers, heat-transfer correlations, and regenerative-cooling solver.
4. `pump` and `turbine` provide machinery-related utilities and can be connected into an engine-cycle model.
5. `common.engine_network` and `common.blocks` allow component models to be connected into a full engine network.
6. `viz` turns the resulting data into plots, reports, and geometry visualisations.

Geometry, material data, gas properties, heat-transfer correlations, numerical solvers, engine-cycle logic, and visualisation are different responsibilities. The idea is that keeping these responsibilities in different subpackages makes the code easier to test, replace, and understand.

## 4. Example versus validation entry

The differences between examples and validation cases are important.

An **example** primarily teaches usage. A good example says:

> “Here is how you use this feature.”

It may use simplified input data. It may prioritise readability over physical correctness. It should be easy to copy, modify, and run.

Examples belong in the `examples/` folder in the root directory.

A **validation entry** primarily tests credibility. A good validation case says:

> “Here is evidence that this feature gives reasonable results for a known case.”

It should explain:

- what reference data is used,
- where the reference data came from,
- what assumptions were made,
- what quantities are compared,
- where the model agrees and where it does not,
- an interpretation of the result.

Validation cases belong in the `validation/` folder in thr root directory:

A validation case can also be educational. After all, they show how the program is used to represent real hardware. But their primary purposes are different:

| Folder | Main question |
|---|---|---|
| `examples/` | “How do I use this?” |
| `validation/` | “What does the result reveal” |

## 5. How to read the code without getting lost

Pyskyfire is mostly written in an **object-oriented** style. This means that the code is organised around objects that represent things in the problem.

An object combines data and behaviour. For example, a thrust chamber is not only a list of numbers. It has geometry, wall layers, cooling circuits, coolant properties, and hot-gas properties. It also has methods that answer questions such as “what is the local channel area?” or “what is the wall thickness here?”

In object-oriented code, you can build larger objects from smaller objects:

```text
Material models
      ↓
Wall layers
      ↓
Cooling-channel cross-section + channel placement + channel height
      ↓
Cooling circuit + Contour
      ↓
Thrust chamber
```

A `ThrustChamber` object therefore acts as a container for many other objects. It may contain a `Contour`, one or more `CoolingCircuit` objects, `Wall` objects, `CoolantTransport` objects, and an aerothermodynamic model for the combustion gas. The solver therefore does not need every detail written directly inside `solver.py`. Instead, it receives the thrust-chamber object and asks it for the geometry, material, gas-property, and coolant-property information it needs.

A value used in the solver may come from another object several layers away. For example, the solver may ask a cooling circuit for its hydraulic diameter, and that cooling circuit may compute it from a cross-section object, a channel-height function, and a channel-placement object.

For regenerative cooling, start with:

| Start here | Why |
|---|---|
| `examples/` or `validation/regen_vali/...` | Shows a complete setup script where the objects are actually created and connected. |
| `regen/thrust_chamber.py` | Shows how the physical thrust chamber, wall layers, and cooling circuits are represented. |
| `regen/solver.py` | Shows the main regenerative-cooling solve procedure. |
| `regen/physics.py` | Shows the local equations and correlations used by the solver. |
| `skycea/aerothermodynamics.py` | Shows where the hot-gas properties come from. |
| `viz/plot_regen.py` | Shows how regenerative-cooling results are plotted. |

For engine-cycle work, start with:

| Start here | Why |
|---|---|
| `common/blocks.py` | Defines the component blocks that read and write stations/signals. |
| `common/engine_network.py` | The underlying orchestrator that connects the blocks. |

## 6. Summary

Pyskyfire is organised like a typical Python engineering package:

- `pyproject.toml` defines how the package is built and installed.
- `.gitignore` keeps local, generated, private, and temporary files out of version control.
- `.gitattributes` helps GitHub classify generated outputs correctly.
- `README.md` introduces the project.
- `docs/` contains the documentation website.
- `examples/` teaches usage.
- `validation/` demonstrates credibility against reference cases.
- `images/` stores figures used in documentation and presentation.
- `src/pyskyfire/` contains the actual importable Python package.

Inside `src/pyskyfire/`:

- `common` contains shared data structures, material/fluid helpers, results containers, and engine-network logic.
- `regen` contains thrust-chamber geometry, cooling-channel definitions, heat-transfer correlations, and the regenerative-cooling solver.
- `skycea` provides CEA-based aerothermodynamics and transport-property interfaces.
- `pump` contains preliminary centrifugal pump and impeller design tools.
- `turbine` contains early turbine-related utilities.
- `viz` contains plotting, reporting, and visualisation tools.

Pyskyfire is a set of cooperating objects and modules. The package separates geometry, properties, physics, solvers, engine-cycle blocks, and visualisation, and orchestrates them together to run a regenerative cooling or full cycle simulation. 
