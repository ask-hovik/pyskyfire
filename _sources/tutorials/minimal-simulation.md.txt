# Minimal Simulation

This tutorial builds and analyses a small regeneratively cooled nitrous-oxide/ethanol rocket engine. It is intended to introduce the principal Pyskyfire objects and the normal workflow for a thrust-chamber thermal analysis. This tutorial is focused on showcasing Pyskyfire capabilities, and is not focused on engine design.

The complete, runnable source is maintained in [`examples/minimal/minimal_sim.py`](https://github.com/ask-hovik/pyskyfire/blob/main/examples/minimal/minimal_sim.py).

## Prerequisites

Install Pyskyfire into your current environemnt, then run:

```console
python examples/minimal/minimal_sim.py
```

The script builds a thrust chamber, and then runs a regenerative cooling analysis. It then postprocesses the results, giving you a few options to view the generated data: Either as standalone html graphs, or as a compiled report output as `minimal_report.html`.

## What you will build

The example uses a 5 kN engine with a 20 bar chamber pressure, an area ratio of 10, nitrous oxide as oxidizer, and ethanol as both fuel and coolant. It uses a single, helical, square-channel cooling circuit running from nozzle exit to chamber inlet.

## Define the engine design point

Start with the intended chamber conditions, thrust, nozzle geometry, coolant inlet state, and the thermal/geometry choices for the cooling circuit.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:engine-inputs
:end-before: tutorial:end:engine-inputs
:dedent: 4  
``` 

`cea_fu` and `cea_ox` define the propellants passed to the NASA CEA-backed combustion model. `coolprop_fu` defines the coolant passed to the CoolProp-backed transport-property model. Different backends are used for the coolant and the combustion gas, NASA CEA for the hot gas, CoolProp for the coolant, hence the need to define the ethanol twice. 

`p_c`, `F`, `eps`, `L_star`, and `MR` establish the thermodynamic design point. `AR_c` represents the chamber aspect ratio, defined as 

$$ AR_c = \frac{L_c^2}{S_c} $$

where $L_c$ is the chamber length, and $S_c$ is the area of the cross section of the chamber when the section plane passes through the chamber axis. 

![Thrust chamber contour](../_static/images/thrust-chamber-contour.svg)

The remaining inputs define the cooling wall and flow channels: material, hot-gas-to-coolant wall thickness, channel count, rib blockage, coolant-side roughness, helix angle and channel height.

## Compute the combustion and coolant models

Pyskyfire uses the supplied design point to construct an aerothermodynamic model of the combustion gases. Upon initialisation, the aerothermodynamics class calculates many parameters about the engine. In this example, the ideal chamber volume, throat radius and fuel mass flow the class has calculated is used further.  

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:aerothermodynamics
:end-before: tutorial:end:aerothermodynamics
:dedent: 4
```

The aerothermodynamics class has multiple constructors allowing you to choose which set of inputs is used to construct the hot gas properties. In this example `from_F_eps_Lstar` is used. See other options in {doc}`Aerothermodynamics <../autoapi/pyskyfire/skycea/aerothermodynamics/Aerothermodynamics>`. The calculated `V_c`, `r_t` and `mdot_fu` are inserted into `params` for convenience. `CoolantTransport` supplies coolant thermodynamic and transport properties.

## Generate the thrust-chamber contour

Next, we are going to create a nozzle contour. We can choose between a conical nozzle and a bell (Rao) nozzle. The resulting axial and radial coordinates are stored in a `Contour` object. 

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:contour
:end-before: tutorial:end:contour
:dedent: 4
```

`get_contour` exposes further shaping parameters beyond those used here. The three curvature factors, `R_1f`, `R_2f`, and `R_3f`, affect the chamber-to-throat and nozzle transitions.

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/tutorial-artifacts/minimal-simulation/contour.html"
    title="Interactive thrust-chamber contour"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
``` 

## Define the walls

Create a `Wall` to represent the material layer separating the coolant from the combustion gases. Multiple wall layers are possible, but we keep it simple here. 

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:walls
:end-before: tutorial:end:walls
:dedent: 4
```

## Define the regenerative-cooling circuit

A `CoolingCircuit` combines the coolant model, channel cross-section, channel placement, wall stack, roughness, and dimensions. Pyskyfire can represent multiple circuits, including circuits covering different spans of the contour or interlacing with one another. This tutorial makes it simple with one circuit running from the nozzle exit to the chamber inlet.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:cooling-circuit
:end-before: tutorial:end:cooling-circuit
:dedent: 4
```

The constant `channel_height_function` creates 2 mm-high channels everywhere on the engine. `SurfacePlacement` maps channels around the thrust-chamber surface. Here the channels follow a constant 45-degree helix. The `span=[1.0, -1.0]` sets the coolant flow direction from the bottom of the nozzle towards the chamber.

## Assemble the thrust chamber

`ThrustChamber` is the composite object that joins the contour, combustion model, and cooling circuit. It computes geometry-dependent quantities that cannot be determined by the individual objects alone.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:thrust-chamber
:end-before: tutorial:end:thrust-chamber
:dedent: 4
```

Before solving the heat-transfer problem, it may be nice to inspect the generated cooling-channel geometry. See the report section for visualisation options. Bellow is one possibility, a 3d representation of the thrust chamber created. 

```{raw} html
<div class="psf-demo-frame psf-demo-frame--large">
  <iframe
    src="../_static/tutorial-artifacts/minimal-simulation/engine-3d.html"
    title="Interactive 3D minimal-engine model"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
``` 

This view is useful for checking that channels cover the intended contour span and have the expected placement before interpreting a thermal result. Unfortunately, rib visualisation is currently not implemented, so cooling channels will always cover the entire chamber wall in the visualsation (but not in simulation). A fix for this will be implemented in the future. 

## Run the steady-state cooling analysis

The cooling simulation needs a coolant inlet temperature, inlet pressure, and mass flow rate. This simple engine assumes that all fuel supplied to the chamber first passes through the cooling circuit, so the coolant mass flow is the fuel mass flow calculated by the combustion model. In the case where there is multiple cooling circuits, the circuit index indicates which one should be used in the given simulation.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:simulation
:end-before: tutorial:end:simulation
:dedent: 4
```

`steady_heating_analysis` returns `cooling_data`, which contains the axial solution for coolant state, wall temperatures, heat flux, flow velocity, and related quantities. Below the heat flux is shown as an example. 

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/tutorial-artifacts/minimal-simulation/heat-flux.html"
    title="Interactive heat flux"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
``` 

## Generate and inspect the report

The example uses Pyskyfire's report system to collect input data, calculated optimal values, geometry plots, cooling results, combustion transport properties, and selected through-wall temperature profiles into a portable HTML file.

```{literalinclude} ../../examples/minimal/minimal_sim.py
:language: python
:start-after: tutorial:start:report
:end-before: tutorial:end:report
:dedent: 4
```

When the scripts completes, a few standalone interactive html graphs have been made. A comprehensive report is compiled into `minimal-report.html`. You can view the report here: 
<a href="../_static/tutorial-artifacts/minimal-simulation/minimal-report.html"> Minimal Simulation Report</a>
