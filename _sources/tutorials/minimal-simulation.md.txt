# Build a Minimal Regeneratively Cooled Engine

This tutorial builds and analyses a small nitrous-oxide/ethanol rocket engine with a single regenerative-cooling circuit. It is intended to introduce the principal Pyskyfire objects and the normal workflow for a thrust-chamber thermal analysis, not to establish a design ready for manufacture or hot-fire testing.

The complete, runnable source is maintained in [`examples/minimal/minimal_sim.py`](https://github.com/ask-hovik/pyskyfire/blob/main/examples/minimal/minimal_sim.py). The code blocks below are included directly from that file, so the tutorial and executable example remain synchronized.

## Prerequisites

Install Pyskyfire and its simulation dependencies. From a source checkout, install the project into the active environment, then run:

```console
python examples/minimal/minimal_sim.py
```

The script writes `minimal_report.html` beside the script. The report includes an
interactive three-dimensional engine view, together
with the simulation results.

## What you will build

The example uses a 5 kN engine with a 50 bar chamber pressure, an area ratio of 10, nitrous oxide as oxidizer, and ethanol as both fuel and coolant. It uses one helical, square-channel cooling circuit from nozzle exit to chamber inlet.

The purpose is to establish the complete chain from a performance design point, through chamber geometry and cooling-channel geometry, to a steady-state thermal calculation and report. It deliberately keeps the cooling geometry simple: a constant 2 mm channel height, a uniform 45-degree helix angle, and one wall material.

## Define the engine design point

Start with the intended chamber conditions, nozzle geometry, combustion mixture ratio, coolant inlet state, and the thermal/geometry choices for the cooling circuit.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:engine-inputs
:end-before: tutorial:end:engine-inputs
:dedent: 4  
``` 

`cea_fu` and `cea_ox` define the propellants passed to the NASA CEA-backed combustion model. `coolprop_fu` defines the coolant passed to the CoolProp-backed transport-property model. The two representations must name the same physical coolant, but use the naming conventions accepted by their respective backends.

`p_c`, `F`, `eps`, `L_star`, and `MR` establish the thermodynamic design point. `AR_c` controls the generated chamber geometry. The remaining inputs define the cooling wall and flow channels: material, hot-gas-to-coolant wall thickness, channel count, rib blockage, and coolant-side roughness.

## Compute the combustion and coolant models

Pyskyfire uses the supplied design point to construct an aerothermodynamic model of the combustion gases. In this example, that model also determines the chamber volume and throat radius needed for the requested thrust, chamber pressure, nozzle area ratio, and characteristic length.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:aerothermodynamics
:end-before: tutorial:end:aerothermodynamics
:dedent: 4
```

The calculated `V_c` and `r_t` are inserted into `params` because the contour generator requires them. `CoolantTransport` supplies coolant thermodynamic and transport properties along the cooling path.

## Generate the thrust-chamber contour and wall

Next, create a Rao nozzle contour and combine it with the chamber dimensions calculated above. The resulting axial and radial coordinates are stored in a `Contour` object. Create a `Wall` to represent the material layer separating the coolant from the combustion gases.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:contour-and-wall
:end-before: tutorial:end:contour-and-wall
:dedent: 4
```

`get_contour` exposes further shaping parameters beyond those used here. The three curvature factors, `R_1f`, `R_2f`, and `R_3f`, affect the chamber-to-throat and nozzle transitions. Retain this simple contour while learning the workflow; contour optimisation is a later design problem.

## Define the regenerative-cooling circuit

A `CoolingCircuit` combines the coolant model, channel cross-section, channel placement, wall stack, roughness, and dimensions. Pyskyfire can represent multiple circuits, including circuits covering different spans of the contour or interlacing with one another. This tutorial uses one circuit running from the nozzle exit to the chamber inlet.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:cooling-circuit
:end-before: tutorial:end:cooling-circuit
:dedent: 4
```

The constant `channel_height_function` creates 2 mm-high channels everywhere on the engine. `SurfacePlacement` maps channels around the thrust-chamber surface. Here the channels follow a constant 45-degree helix. The `span=[1.0, -1.0]` sets the coolant flow direction from the bottom of the nozzle towards the chamber.

## Assemble and inspect the thrust chamber

`ThrustChamber` is the composite object that joins the contour, combustion model, and cooling circuit. It computes geometry-dependent quantities that cannot be determined by the individual objects alone.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:thrust-chamber
:end-before: tutorial:end:thrust-chamber
:dedent: 4
```

Before solving the heat-transfer problem, it may be nice to inspect the generated cooling-channel geometry.

This view is useful for checking that channels cover the intended contour span and have the expected placement before interpreting a thermal result.

## Run the steady-state cooling analysis

The cooling simulation needs a coolant inlet temperature, inlet pressure, and mass flow rate. This simple engine assumes that all fuel supplied to the chamber first passes through the cooling circuit, so the coolant mass flow is the fuel mass flow calculated by the combustion model.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:simulation
:end-before: tutorial:end:simulation
:dedent: 4
```

`steady_heating_analysis` returns `cooling_data`, which contains the axial solution for coolant state, wall temperatures, heat flux, flow velocity, and related quantities. More sophisticated engine studies may solve the cooling-loop inlet conditions as part of a complete engine-cycle iteration rather than prescribing them directly.

## Generate and inspect the report

The example uses Pyskyfire's report system to collect input data, calculated optimal values, geometry plots, cooling results, combustion transport properties, and selected through-wall temperature profiles into a portable HTML file.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:report
:end-before: tutorial:end:report
:dedent: 4
```

Open `minimal_report.html` after the script completes. The first checks should be:

- **Wall temperature:** locate the hot-side wall-temperature peak, particularly near the throat.
- **Coolant pressure:** confirm that the pressure remains above a suitable margin for the coolant state and your feed-system assumptions.
- **Coolant temperature and velocity:** check for implausible temperature rise, acceleration, or geometry discontinuities.
- **Heat flux:** compare its axial distribution with the expected throat-region thermal peak.
- **Cooling-channel area and hydraulic diameter:** verify that the generated channel geometry is physically credible along the entire contour.

This is a deliberately simple baseline. It does not size injectors, pumps, turbines, tanks, valves, structural margins, ignition hardware, or a full feed cycle. The next tutorial extends this workflow into a full engine-cycle study.

## Generated report

The documentation build runs the complete example and publishes the resulting report
with this tutorial.

<a href="../_static/reports/minimal-report.html">Open the generated minimal-engine report</a>