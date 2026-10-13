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
| `pyproject.toml` | The main configuration file for the Python package. It defines the package name, version, core and optional dependencies, developer dependency groups, build system, and tells the build system that the importable package lives in `src/pyskyfire`. |
| `docs/` | Contains the documentation website. This includes tutorials, how-to guides, explanations such as this page, and generated API reference pages. |
| `examples/` | Contains scripts that demonstrate how to use the package. These are primarily teaching material: they show workflows and interfaces that users can copy and adapt. |
| `images/` | Stores figures used in the README and documentation, such as plots, screenshots, and rendered engine or cooling-channel images. |
| `src/` | Contains the actual installable Python package. When a user writes `import pyskyfire`, Python imports code from `src/pyskyfire`. |
| `validation/` | Contains cases that compare Pyskyfire against reference data, published results, known solutions, or other tools. These are mainly about credibility rather than teaching the interface. |


## 2. `src/pyskyfire/`

The `src/pyskyfire/` directory contains the code that becomes the `pyskyfire` Python package after installation.

The package content is as follows:

```text
src/pyskyfire/
├── __init__.py
├── common/
├── network/
├── regen/
├── thermo/
├── pump/
├── turbine/
├── heat_exchanger/
├── ducts/
├── injector/
└── viz/
```

Each subfolder is a *subpackage*. A subpackage is a folder containing related Python modules. A module is a single `.py` file.

For example:

```text
src/pyskyfire/regen/coupled_solver.py
```

is a module, and it belongs to the subpackage:

```python
pyskyfire.regen
```

### `__init__.py`

The top-level file `src/pyskyfire/__init__.py` marks `pyskyfire` as a Python package and exposes the main subpackages.

### `pyskyfire.common`

`common` holds dependency-free primitives that any other subpackage can use. It does not import from `regen`, `ducts`, `pump`, `turbine`, `network`, or `viz` — this keeps it a true leaf package and avoids circular imports. Fluid definitions and property models belong to `thermo`.

| Module | Role |
|---|---|
| `solids.py` | Defines material-property models for solids. These can represent constant, polynomial, tabulated, or piecewise temperature-dependent properties such as thermal conductivity. |
| `results.py` | Defines a dictionary-like `Results` container with dot access and save/load functionality. This is useful for storing simulation outputs in a convenient way. |
| `section_geometry.py` | Defines 2D cross-section primitives (circular and elliptical arcs, lines, boundary construction) shared by cooling-channel geometry in `regen` and duct/volute geometry in `ducts`. |
| `physics.py` | Shared heat-transfer functions: Bartz, Colburn, Gnielinski and recovery temperature. |
| `friction.py` | Shared Darcy factors and continuous-curvature corrections. |
| `momentum.py` | Shared reversible acceleration pressure change. |
| `blade_surface.py` | Defines `BladeSurface`, one blade's suction and pressure sides in cylindrical coordinates, shared by impellers and turbines. |
| `meridional.py` | Defines `MeridionalPassage`, the wetted flow path of a turbomachine or annular duct as paired hub and shroud curves, with `WallMotion` (which wall stretches rotate) and `Clearance` (where blades stand off a wall). |
| `constants.py` | Stores shared constants that may be used across modules. Keeping constants in one place avoids repeating the same numbers throughout the code. |
| `__init__.py` | Re-exports selected objects from the `common` subpackage so they are easier to import. |

### `pyskyfire.network`

`network` connects component models from other subpackages into a full engine-cycle calculation. Unlike `common`, it sits above the other subpackages: it is free to import from `regen`, `ducts`, `thermo`, and `common`.

| Module | Role |
|---|---|
| `engine_network.py` | Defines `Station` objects and the `EngineNetwork` fixed-point framework. This is used to connect pumps, turbines, cooling circuits, and other components into a larger engine-cycle calculation. |
| `blocks.py` | Defines the block interface used by the engine network. A block reads stations and/or scalar signals, performs a calculation, and returns updated stations and signals. |
| `__init__.py` | Re-exports selected objects from the `network` subpackage so they are easier to import. |

### `pyskyfire.regen`

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
| `common.physics` (imported) | Shared engineering correlations used by the solver. |
| `solver.py` | Runs the regenerative-cooling calculation. It marches along the cooling circuit, solves local wall heat balances, and updates coolant temperature and pressure. |
| `coolant_properties.py` | Defines the regenerative-cooling film-property policy through `film_properties(resolver, T_film, p, getters)`, separate from the general fluid-state resolver in `thermo.state`. |
| `film_solver.py` | Contains film-cooling model functionality. Film cooling is separate from regenerative cooling because coolant is injected into the hot-gas side rather than only flowing inside wall channels. |
| `film_solver_2.py` | Contains a newer or alternative film-cooling implementation. New users should normally start from documented examples or validation cases before relying on experimental modules. |
| `__init__.py` | Re-exports selected objects from the `regen` subpackage so common classes and functions are easier to import. |

### `pyskyfire.thermo`

`thermo` owns fluid definitions, thermodynamic state resolution, transport properties, and combustion-gas aerothermodynamics. It connects engine geometry and propellant choices to the fluid properties needed by the rest of the package; coolant-specific heat-transfer policy remains in `regen`.

| Module | Role |
|---|---|
| `fluid.py` | Defines `Fluid` for propellant and coolant mixtures, including component names, composition fractions, and mass or mole basis. |
| `transport.py` | Defines `FluidTransport` and `TransportProperties`. Properties such as density, heat capacity, viscosity, conductivity, and Prandtl number can be constants or functions of temperature and pressure. |
| `state.py` | Defines `FluidStateResolver`, the resolved `FluidState`, `FluidCapability`, and `SaturationState`; `probe_fluid` detects the available fluid capabilities. |
| `eos.py` | Defines the `EquationOfState` interface a machine design reads a fluid through, with `Incompressible` (constant density) and `CoolPropEOS` (a real fluid). |
| `duty.py` | Defines `MachineDuty`, the head, flows and speed asked of a turbomachine, which `pump.PumpDuty` and `turbine.TurbineDuty` build on. |
| `aerothermodynamics.py` | Defines `Aerothermodynamics` and `CombustionGasState`. It computes and stores CEA-based hot-gas properties along the chamber/nozzle contour for the regenerative-cooling solver. |
| `nozzle_loss.py` | Defines `NozzleLoss` for nozzle performance-loss calculations. |
| `__init__.py` | Publicly exports `Fluid`, `FluidTransport`, `TransportProperties`, `FluidStateResolver`, `FluidState`, `FluidCapability`, `SaturationState`, `probe_fluid`, `EquationOfState`, `Incompressible`, `CoolPropEOS`, `MachineDuty`, `Aerothermodynamics`, `CombustionGasState`, and `NozzleLoss`. |

#### Thermo migration (unreleased)

Update imports directly; no backwards aliases are provided. Paths below are
relative to `pyskyfire`:

| Previous import | Current import |
|---|---|
| `skycea` | `thermo` |
| `common.fluids.Fluid` / `common.Fluid` | `thermo.Fluid` (defined in `thermo.fluid`) |
| `skycea.coolant_transport.CoolantTransport` | `thermo.transport.FluidTransport` |
| `skycea.coolant_transport.TransportProperties` | `thermo.transport.TransportProperties` |
| `regen.coolant_state.CoolantState` | `thermo.state.FluidStateResolver` |
| `regen.coolant_state.BulkState` | `thermo.state.FluidState` |
| `regen.coolant_state.CoolantCapability` | `thermo.state.FluidCapability` |
| `regen.coolant_state.SaturationState` | `thermo.state.SaturationState` |
| `regen.coolant_state.probe_coolant` | `thermo.state.probe_fluid` |
| `skycea.aerothermodynamics.GasState` | `thermo.aerothermodynamics.CombustionGasState` |
| `skycea.aerothermodynamics.Aerothermodynamics` | `thermo.aerothermodynamics.Aerothermodynamics` |
| `skycea.nozzle_loss.NozzleLoss` | `thermo.nozzle_loss.NozzleLoss` |

All classes in the table are also available directly from `pyskyfire.thermo`.
Film-property evaluation is no longer a method on the state resolver: use
`pyskyfire.regen.coolant_properties.film_properties(resolver, T_film, p, getters)`.
The `coolant_transport` argument on cooling circuits and film inputs retains its
name; pass a `thermo.FluidTransport` instance.

Old stored pickle files may require regenerating because modules and class names
have changed. Compatibility with those files is not guaranteed.

### `pyskyfire.ducts`

`ducts` holds duct hydraulics that no single component owns: the friction
physics shared by cooling channels and piping, explicit pipe geometry,
annular ducts, and volute and collector manifolds. It depends only on `common` for planar section
primitives and on `thermo` for fluid state records, and it imports from neither
`regen` nor `network`, so cooling channels, turbopump hardware and engine-cycle
adapters can all build on it without a circular import.

| Module | Role |
|---|---|
| `path.py` | Shared path nodes, stations, connections, sampling of curves into paths (`Path.from_curve`), rigid placement, and transported frames. |
| `common.friction` (imported) | Geometry-agnostic friction correlations shared with regen and heat exchangers. |
| `pipe.py` | Explicit pipe geometry, with an optional core, and the loss march every duct shares: friction, bends, area change, and swirl. Returns a `DuctResult`. |
| `annular.py` | Annular ducts about the shaft (`AnnularDuct`, on `common.MeridionalPassage`), solved by the shared march. |
| `orifice.py` | Orifices (round, annular, slot) discharging into a volume, with isentropic flow that chokes and cavitates, and ISO 5167 orifice plates in a pipe. |
| `volute_section.py` | Volute cross-section families: an elliptical bore at fixed inlet slot width, sized by one parameter and reshaped by a second. |
| `area_law.py` | Rules relating azimuth to passage size: conservation of angular momentum for a shaped volute, uniform area for a collector. |
| `exit_cone.py` | The duct carrying a volute's collected flow to its pipe flange, at any angle between tangential and radial and any elevation out of the volute's own plane. |
| `volute.py` | Volute geometry generation and Aungier's one-dimensional loss analysis. |

### `pyskyfire.injector`

`injector` holds injector hardware and the flow through it. It depends on
`ducts`, whose `Orifice` is every hole, and on `thermo`; the engine network
adapts it through `network.InjectorBlock`. Each element type is its own class.

| Module | Role |
|---|---|
| `element.py` | `Element`, the base every element type implements: its holes, placed in the element's frame, and its mixing diagnostics. |
| `impinging.py` | `Jet` and `ImpingingElement`: jets aimed at a common point, covering like and unlike doublets, triplets and pentads. |
| `injector.py` | `ElementSet` (one element type at supplied face positions), `Injector`, and the per-feed hydraulics: each feed solved as one equivalent orifice, for either the flow or the plenum pressure. |

### `pyskyfire.pump` / `pyskyfire.turbine`

The `pump` and `turbine` subpackages contain tools for turbopumps. They are both under development, and currently the pyskyfire package works fine without them. `pump` sizes radial impellers and generates their blade geometry. `turbine` designs axial turbines at mean radius (`AxialTurbine`), describes them as hardware (`TurbineDefinition`, `BladeRow`) and generates their blades. `pump.PumpDuty` and `turbine.TurbineDuty` state what each machine is asked to do. Volute manifolds live in `pyskyfire.ducts`, and the meridional passage both machine types share lives in `pyskyfire.common`.

### `pyskyfire.viz`

`viz` turns simulation data into plots, reports, and geometry views. It is where the package moves from numerical results to things a user can inspect visually.

| Module | Role |
|---|---|
| `core.py` | Defines common Plotly helper functionality used by several plotting classes. It provides shared behaviour for configuring, showing, and saving figures. |
| `plot_regen.py` | Plotting tools for regenerative-cooling results, such as wall temperature, coolant temperature, pressure, heat flux, and heat-transfer coefficients. |
| `plot_film_cooling.py` | Plotting tools for film-cooling results. |
| `plot_thermo.py` | Plotting tools for aerothermodynamic and CEA-derived results. |
| `plot_common.py` | Shared plotting utilities used by several plotting modules. |
| `engine_viz.py` | 3D engine and cooling-channel visualisation tools. This includes geometry-oriented views of cooling circuits and engine structure. |
| `impeller_viz.py` | Visualisation tools for pump and impeller geometry. |
| `components.py` | Adapts engineering objects — thrust chambers, impellers, pipes, annular ducts, volutes — into scene hierarchies. |
| `report.py` | Report-generation helpers for collecting plots, tables, and simulation summaries. |
| `lookup_tables.py` | Lookup-table support used by plotting or reporting workflows. |
| `__init__.py` | Re-exports selected visualisation functions and classes for easier imports. |

## 3. How the subpackages fit together

The subpackages each have a responsibility, but in order to run a regenerative cooling analysis or a full cycle analysis, the packages work together. A typical regenerative-cooling analysis might work like this:

1. `common` defines shared, dependency-free primitives such as solids, results, and section geometry.
2. `thermo` defines fluids and provides gas and coolant state and property models.
3. `regen` defines the thrust chamber, cooling channels, wall layers, heat-transfer correlations, and regenerative-cooling solver.
4. `pump` and `turbine` provide machinery-related utilities and can be connected into an engine-cycle model.
5. `network.engine_network` and `network.blocks` allow component models to be connected into a full engine network.
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

A `ThrustChamber` object therefore acts as a container for many other objects. It may contain a `Contour`, one or more `CoolingCircuit` objects, `Wall` objects, `FluidTransport` objects, and an aerothermodynamic model for the combustion gas. The solver therefore does not need every detail written directly inside `solver.py`. Instead, it receives the thrust-chamber object and asks it for the geometry, material, gas-property, and coolant-property information it needs.

A value used in the solver may come from another object several layers away. For example, the solver may ask a cooling circuit for its hydraulic diameter, and that cooling circuit may compute it from a cross-section object, a channel-height function, and a channel-placement object.

For regenerative cooling, start with:

| Start here | Why |
|---|---|
| `examples/` or `validation/regen_vali/...` | Shows a complete setup script where the objects are actually created and connected. |
| `regen/thrust_chamber.py` | Shows how the physical thrust chamber, wall layers, and cooling circuits are represented. |
| `regen/coupled_solver.py` | Shows the regenerative and film-coupled cooling solve procedure. |
| `common/physics.py` | Shows the local equations and correlations used by the solver. |
| `thermo/aerothermodynamics.py` | Shows where the hot-gas properties come from. |
| `viz/plot_regen.py` | Shows how regenerative-cooling results are plotted. |

For engine-cycle work, start with:

| Start here | Why |
|---|---|
| `network/blocks.py` | Defines the component blocks that read and write stations/signals. |
| `network/engine_network.py` | The underlying orchestrator that connects the blocks. |

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

- `common` contains dependency-free shared primitives: solid-material helpers, results containers, and section geometry.
- `network` contains the engine-cycle blocks and the `EngineNetwork` orchestrator that connects them.
- `regen` contains thrust-chamber geometry, cooling-channel definitions, heat-transfer correlations, and the regenerative-cooling solver.
- `thermo` provides fluid definitions, state resolution, transport-property interfaces, and CEA-based aerothermodynamics.
- `pump` contains preliminary centrifugal pump and impeller design tools.
- `turbine` contains early turbine-related utilities.
- `viz` contains plotting, reporting, and visualisation tools.

Pyskyfire is a set of cooperating objects and modules. The package separates geometry, properties, physics, solvers, engine-cycle blocks, and visualisation, and orchestrates them together to run a regenerative cooling or full cycle simulation. 

### `pyskyfire.heat_exchanger`

`circuit.py` defines channel banks, prescribed gas, and explicit thermal/flow
connections. `solver.py` couples segment heat rates to enthalpy and static
pressure. `chamber.py` places assemblies in prescribed chamber flow and reports
displaced volume. Paths come from `ducts.path`, sections from `regen.cross_section`,
and scalar correlations from `common`. See [Heat exchangers](../howto/heat-exchangers.md).
