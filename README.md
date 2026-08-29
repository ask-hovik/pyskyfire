![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/pyskyfire_header.png)

**Pyskyfire is a simulation framework for regeneratively cooled, liquid propellant rocket engines.**

Model thrust-chamber thermodynamics, cooling channels, engine cycles, and get interactive engineering reports from a unified Python interface.

------------------
[![PyPI version](https://img.shields.io/pypi/v/pyskyfire.svg)](https://pypi.org/project/pyskyfire/)
[![Python versions](https://img.shields.io/pypi/pyversions/pyskyfire.svg)](https://pypi.org/project/pyskyfire/)
[![Tests](https://img.shields.io/github/actions/workflow/status/ask-hovik/pyskyfire/tests.yml?branch=main&label=tests)](https://github.com/ask-hovik/pyskyfire/actions/workflows/tests.yml)
[![Docs](https://img.shields.io/github/actions/workflow/status/ask-hovik/pyskyfire/docs.yml?branch=main&label=docs)](https://github.com/ask-hovik/pyskyfire/actions/workflows/docs.yml)
[![License](https://img.shields.io/github/license/ask-hovik/pyskyfire.svg)](https://github.com/ask-hovik/pyskyfire/blob/main/LICENSE)

------------------
<p align="center">
  <a href="https://ask-hovik.github.io/pyskyfire/"><strong>Documentation</strong></a>
  ·
  <a href="https://ask-hovik.github.io/pyskyfire/tutorials/minimal-simulation.html"><strong>Quick Start</strong></a>
  ·
  <a href="https://pypi.org/project/pyskyfire/"><strong>PyPI</strong></a>
</p>

<p align="center">
  <a href="https://ask-hovik.github.io/pyskyfire/validation/engine-3d.html">
    <img src="https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/RL10_turntable.webp" width="820" alt="Rotating 3D view of the RL10A-3-3A thrust chamber and its regenerative cooling tubes">
  </a>
</p>
<p align="center">
  <em>RL10A-3-3A thrust chamber reconstructed in Pyskyfire. 180 long tubes (red) run the
  full chamber and throat; 180 short tubes (dark) start at the nozzle interlacing manifold.
  Downstream of it the two sets alternate to give 360 passages, while only the long tubes
  remain upstream.<br>
  <a href="https://ask-hovik.github.io/pyskyfire/validation/engine-3d.html">Open the interactive model</a> to rotate, pan, and zoom it yourself.</em>
</p>


# Description
Pyskyfire is an open-source python package, meant as an alternative to RPA, NPSS, ESPSS and other regenerative cooling and engine cycle analysis software. The solver produces regenerative cooling and engine cycle results. Compare NASA and Pratt & Whitney with Pyskyfire in this RL10 case: 

<p align="center">
  <a href="https://ask-hovik.github.io/pyskyfire/_static/readme-artifacts/rl10-wall-temperature.html">
    <img src="https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/RL10_wall_temperature.png" width="820" alt="Predicted RL10A-3-3A hot-wall temperature compared with three published reference models">
  </a>
</p>
<p align="center">
  <em><a href="https://ask-hovik.github.io/pyskyfire/_static/readme-artifacts/rl10-wall-temperature.html">Open the interactive chart</a>
  · <a href="https://ask-hovik.github.io/pyskyfire/validation/rl10a-3-3a.html">Read the full RL10A-3-3A validation report</a></em>
</p>

The first iteration of pyskyfire was written as part of the master thesis of Ask Haugerud Hovik, which can be read [here](https://drive.google.com/file/d/1sZJmt-8UWtUChprji67LmnazS3Ei_K3a/view). The motivation to start writing the software came purely from a curiosity standpoint and from an innate wish to spread the understanding of rocket engines and propel us further into the space age. Please use this software responsibly and make sure you, your team memebers and everyone else stay safe in your rocket engine endeavours.

# Program Capabilities
- Chemical-equilibrium and hot-gas property calculations
- Rao nozzle and combustion-chamber contour generation
- Multi-pass regenerative-cooling analysis
- Configurable channel geometry and wall materials
- Pump, turbine, and engine-cycle component models
- Interactive 3D visualizations and HTML engineering reports
- Validation cases based on published engine data 
- View the full overview of the [capabilities of pyskyfire](https://ask-hovik.github.io/pyskyfire/explanations/capabilities.html).

# Documentaion
[Documentation](https://ask-hovik.github.io/pyskyfire) for the project is available, and written in a pedagogic style. The package with examples and validation cases is still being developed, so the documentation does not cover everything the package is capable of. Look into the validation cases and advanced examples to see everything the package can do.

# Installation
The core package is available on PyPI and can be installed without the optional visualisation dependencies:

```
uv pip install pyskyfire
```

To include the interactive plots, 3D visualisations, and HTML reports, install the `viz` extra:

```
uv pip install "pyskyfire[viz]"
```

For the bleeding edge version, clone the repository and install the development environment. This includes the visualisation, test, and documentation dependencies.


# Contributing

Pyskyfire began as my master's thesis and is now developed in the open. I would like it to
become something students and professionals build together, and I am happy to help
contributors who are new to propulsion, to open source, or to both.


[CONTRIBUTING.md](https://github.com/ask-hovik/pyskyfire/blob/main/CONTRIBUTING.md) covers
development setup, testing, and the pull-request process. For anything substantial, please
open an [issue](https://github.com/ask-hovik/pyskyfire/issues) first so the approach can be
discussed — and if you would rather talk it through before writing any code, that is welcome
too.

# Using the package
What is it like to use the package? Pyskyfire is an object oriented library, so using it means creating a lot of objects and putting them together, for example: 

```python
import numpy as np
import pyskyfire as psf

# A wall, a channel cross section, and how the channels wrap the contour.
wall = psf.regen.Wall(material=psf.common.solids.StainlessSteel304, thickness=0.5e-3)
cross_section = psf.regen.CrossSectionSquared(blockage_ratio=0.1)
placement = psf.regen.SurfacePlacement(
    n_channel_positions=60,
    helix_angle=lambda x: np.radians(45),
)

# One regenerative cooling circuit, running from the nozzle exit to the injector.
cooling_circuit = psf.regen.CoolingCircuit(
    name="Cooling Pass",
    contour=contour,                      # from psf.regen.contour.get_contour(...)
    coolant_transport=coolant_transport,  # CoolProp-backed coolant properties
    cross_section=cross_section,
    span=[1.0, -1.0],
    placement=placement,
    walls=[wall],
    roughness=10e-6,
    channel_height=lambda x: 2e-3,
)

# The thrust chamber ties geometry, hot-gas properties, and cooling circuits together.
thrust_chamber = psf.regen.ThrustChamber(
    contour=contour,
    combustion_transport=aerothermodynamics,  # from psf.skycea.Aerothermodynamics
    cooling_circuits=[cooling_circuit],
    n_nodes=150,
)
```

Then this object can be used in a variety of ways. For the thrust chamber we created above, for example: 

```python
# Solve the coupled hot-gas, wall, and coolant problem along the circuit.
boundary_conditions = psf.regen.BoundaryConditions(
    T_coolant_in=273.15,
    p_coolant_in=23e5,
    mdot_coolant=aerothermodynamics.mdot_fu,
)

cooling_data = psf.regen.coupled_steady_heating_analysis(
    thrust_chamber,
    nodes=100,
    circuit_index=0,
    boundary_conditions=boundary_conditions,
)

# Plot a single result...
psf.viz.PlotWallTemperature(
    cooling_data, plot_hot=True, plot_coolant_wall=True
).save_html("wall-temperature.html")

# ...or collect several into one standalone, interactive HTML report.
report = psf.viz.Report("Minimal Engine")
tab = report.add_tab("Cooling Data")
tab.add_figure(psf.viz.PlotHeatFlux(cooling_data))
tab.add_figure(psf.viz.PlotCoolantTemperature(cooling_data))
report.save_html("engine-report.html")
```
