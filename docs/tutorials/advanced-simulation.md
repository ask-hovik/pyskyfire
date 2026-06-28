# Engine-Cycle Simulation

This tutorial develops the engine-cycle layer of the advanced methane/oxygen dual-expander example. It assumes that you are already comfortable with the thrust-chamber, cooling-circuit, and regenerative heat-transfer objects introduced in {doc}`Minimal Simulation <minimal-simulation>`.

The complete example consists of two scripts:

```text
examples/advanced/
├── sizer_sim.py
└── post_process.py
```

`RUN_MODE = "regen_only"` is useful while developing the thrust chamber: it runs the regenerative circuits with explicit inlet boundary conditions and does not construct the complete cycle. This tutorial concentrates on `RUN_MODE = "full_cycle"`, where the cooling circuits are embedded in a coupled pump–regenerator–turbine–injector network.

Run the simulation and report generator from the repository root:

```console
python examples/advanced/sizer_sim.py
python examples/advanced/post_process.py
```

The simulation writes `results.pkl`; post-processing reads that file and writes `methane_engine_report.html`.

## Choose the run mode

The top of the simulation script selects the calculation to perform.

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:run-mode
:end-before: tutorial:end:run-mode
```

Use `"regen_only"` when you are changing chamber geometry, walls, cooling-channel layout, or standalone cooling boundary conditions. Use `"full_cycle"` only once the regenerative configuration is plausible enough to be embedded in the cycle.


## Define the design point

`make_params()` contains the editable numerical and physical inputs. It holds four kinds of data:

- chamber design point, propellants, contour parameters, wall stack, and cooling geometry;
- tank pressures, pump efficiencies and speeds, turbine efficiencies, and duct pressure ratios;
- small recirculation and turbine-bypass fractions;
- independent inlet conditions for the standalone regeneration mode.

The combustion propellants and coolant propellants are both defined because the hot-gas model and coolant-property model have distinct responsibilities. The parameter dictionary serves as the single source of design assumptions: it is saved together with the result and displayed in the report.

## Build the thrust chamber

The full cycle still begins by constructing a thrust chamber. The block below computes the combustion-gas transport model, generates the contour, defines walls and cooling circuits, and combines them into one `ThrustChamber`. (In principle, the thrust chamber aerothermodynamic properties are dependent on the condition of the propellants at the inlet of the engine. It could therefore be a part of the simulation loop. However, the computational cost to this is great, and the added precision it gives is marginal. It is therefore more practical to guess the inlet conditions to the thrust chamber, and perhaps update them once the inlet conditions to the thrust chamber has been established by running the simulation). 

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:thrust-chamber
:end-before: tutorial:end:thrust-chamber
```

This example contains three sequential fuel cooling passes and one oxidizer pass. The fuel first cools the copper throat/chamber region, then travels through a co-current nozzle pass and returns through a counter-current nozzle pass. The oxidizer is split into identical parallel paths through the oxidizer cooling circuit.

The chamber construction is intentionally presented as one block here. For a step-by-step explanation of the combustion model, contour generation, walls, channel placement, and standalone regenerative solution, see the {doc}`Minimal Simulation <minimal-simulation>` tutorial.

## Initial station guesses

A full engine cycle has fluid states at every component interface. Pyskyfire represents each state as a `Station(p, T, mdot)`, containing pressure in Pa, temperature in K, and mass flow in kg/s.

The network needs an initial station dictionary before it can converge:

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:initial-stations
:end-before: tutorial:end:initial-stations
```

The station names are not merely labels. They are the links between blocks: a block reads its named inlet stations and writes its named outlet stations. For example, the fuel pump writes `fu_pump_out`; the recirculation splitter then reads `fu_pump_out` and writes `fu_regen_duct_in` and `fu_shaft_recirc`.

These values are **initial guesses**, not prescribed results. They should nevertheless be physically credible. A useful starting point follows the expected trend through the engine:

- pressure rises across a pump and decreases across ducts, cooling circuits, turbines, and injectors;
- temperature rises through regenerative cooling and decreases across a turbine;
- mass flow follows the intended split, merge, recirculation, and bypass fractions.

A fixed-point solver is not a global optimiser that can reliably recover from arbitrary guesses. If an early sweep sends a fluid to an invalid thermodynamic state, makes a turbine outlet temperature non-physical, or generates an enormous pressure change, property calls and the regenerative solver may fail before the network has a chance to settle. Use rough but reasonable values; they do not need to be accurate final-cycle predictions.

## Initial scalar signals

Stations represent flowing propellant. Signals represent scalar quantities shared between blocks but not carried by a single fluid stream: pump power, turbine power requirement, pressure drops, targets, or future control variables.

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:initial-signals
:end-before: tutorial:end:initial-signals
```

This example begins with estimates for the fuel and oxidizer pump powers. The transmission blocks later update the corresponding turbine power requirements. Pressure-drop signals are added separately after the blocks are created because each loss-producing block declares its own `dp_key`.

## Create the engine network

The network setup first gathers the station guesses, scalar signals, and component blocks.

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:network-setup
:end-before: tutorial:end:network-setup
```

Each block advertises the station and signal keys it consumes and emits. The `EngineNetwork` executes the blocks in the order supplied by the list, merges each block's outputs into its station and signal dictionaries, and measures the maximum relative change in the updated quantities.

At present, Pyskyfire preserves the supplied block order. It does not yet perform a dependency-aware topological sort. Arrange the blocks in a valid flow order: merge or source, pump, splitter, ducts and cooling passages, turbine, downstream ducts, injector, then transmission coupling.

### Fluid and signal blocks

The example uses the following block types.

| Block | Role in the network |
|---|---|
| `MassFlowMergerBlock` | Combines same-fluid branches by mass and enthalpy. |
| `PumpBlock` | Raises pressure to meet its required load and emits the pump shaft-power signal. |
| `MassFlowSplitterBlock` | Divides a stream into fixed fractions; the branch pressure and temperature are unchanged by the ideal split itself. |
| `SimpleDuctBlock` | Applies a fixed pressure ratio and an adiabatic, constant-enthalpy pressure loss. |
| `RegenBlock` | Runs a regenerative-cooling circuit using its inlet station, writes its outlet station, and emits its pressure drop. |
| `TurbineBlock` | Expands the inlet stream enough to deliver the shaft-power signal requested by the transmission. |
| `TransmissionBlock` | Sums shaft power demands and writes the power signal consumed by the turbine. |

### Fuel-side blocks and the pump load

The fuel-side block sequence reads much like an engine flow schematic:

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:fuel-side-blocks
:end-before: tutorial:end:fuel-side-blocks
```

`PumpBlock.overcome` deserves special attention. It lists the downstream blocks whose pressure losses the pump must overcome:

```python
overcome=[
    "duct_pump_regen_fuel",
    "regen_throat_pass",
    "duct_regen_turbine_fuel",
    "fuel_turbine",
    "duct_turbine_injector_fuel",
    "fu_injector",
]
```

A pressure-loss block named `regen_throat_pass` writes a signal named `dp_regen_throat_pass`. The pump collects `dp_<block name>` for every listed item, adds the baseline chamber-side pressure requirement, and determines the target pump outlet pressure. Its required shaft power is then emitted as `P_fuel_pump`.

This is somewhat clunky. The pressure path is not inferred automatically from the network graph, so modifying the cycle topology requires manually updating `overcome`. It remains the current solution because it keeps the block models local and the fixed-point solver simple: the pump only needs scalar pressure-drop signals, rather than a general graph traversal and algebraic loop formulation. Treat each `overcome` list as part of the cycle definition and review it whenever a component is added, removed, bypassed, or moved.

### Oxidizer-side blocks

The oxidizer side follows the same architecture, with a split into two parallel cooling branches before the turbine:

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:oxidizer-side-blocks
:end-before: tutorial:end:oxidizer-side-blocks
```

The two oxidizer `RegenBlock` instances deliberately use the same name because they represent identical parallel paths. They therefore produce the same branch pressure-drop signal, which is the relevant loss for the pump's `overcome` list.

### Seed pressure-drop signals and solve

Before the first iteration, the script creates initial pressure-drop signals from the station guesses.

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:pressure-drop-signals
:end-before: tutorial:end:pressure-drop-signals
```

This bootstrap step exists because the pumps read downstream pressure-drop signals on their first sweep, while those signals are only calculated by the loss-producing blocks during the sweep. The initial values do not need to be exact, but they should be non-negative and consistent with the station guesses.

The network is then constructed and solved:

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:network-solve
:end-before: tutorial:end:network-solve
```

During one fixed-point sweep, every block receives the current station and signal dictionaries, computes its outputs, and overwrites the corresponding dictionary entries. Pyskyfire records the largest relative update among pressure, temperature, mass flow, and scalar signals. It repeats complete sweeps until that residual falls below `tol` or the iteration limit is reached.

A converged network performs one additional **post-process sweep**. Every block has a `post_process(stations, signals)` hook. Most blocks return an empty dictionary because their ordinary station result is already sufficient. `RegenBlock` uses this hook to rerun its cooling calculation on a detailed axial grid using the final converged inlet condition. The resulting temperature, pressure, heat-flux, and wall-temperature profiles are collected in `net.block_results`, keyed by block name.

This split is important: a network iteration only needs enough information to update the coupled cycle. Detailed axial arrays are more expensive and are therefore generated after convergence, once rather than once per fixed-point sweep.

## Normalize cooling results and save the simulation

The standalone regeneration path returns named cooling data directly. The full-cycle path obtains equivalent cooling data from the regenerative blocks' post-process outputs:

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:full-cycle-cooling-data
:end-before: tutorial:end:full-cycle-cooling-data
```

The main function selects the mode, runs the requested calculation, and saves a portable `Results` object.

```{literalinclude} ../../examples/advanced/sizer_sim.py
:language: python
:start-after: tutorial:start:run-and-save
:end-before: tutorial:end:run-and-save
```

Both modes save the common result contract:

```text
mode
params
thrust_chamber
cooling_data
```

A full-cycle result also stores:

```text
net
stations
signals
residuals
block_results
```

Saving results separates the expensive calculation from visualisation. You can run the cycle once, inspect the residual history, then adjust graph selection, captions, plot order, or report layout repeatedly without rerunning the thermal solver and fixed-point cycle.

## Generate the report

`post_process.py` loads `results.pkl` and first generates report tabs that are meaningful for either run mode: input parameters, thrust-chamber geometry, cooling data, combustion properties, and through-wall thermal gradients.

```{literalinclude} ../../examples/advanced/post_process.py
:language: python
:start-after: tutorial:start:common-report
:end-before: tutorial:end:common-report
```

For a full-cycle result, the post-processing script adds a dedicated cycle tab containing the convergence history, station pressure/temperature/mass-flow plots, and fuel/oxidizer pressure-temperature paths.

```{literalinclude} ../../examples/advanced/post_process.py
:language: python
:start-after: tutorial:start:cycle-report
:end-before: tutorial:end:cycle-report
```

Finally, the script validates the saved result contract, chooses the appropriate report content from `mode`, and writes the HTML report.

```{literalinclude} ../../examples/advanced/post_process.py
:language: python
:start-after: tutorial:start:generate-report
:end-before: tutorial:end:generate-report
```

A regeneration-only result produces a thrust-chamber and cooling report without an engine-cycle tab. A full-cycle result includes both the chamber analysis and the coupled-cycle diagnostics.
