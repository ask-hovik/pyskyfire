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
Pyskyfire is an open-source python package, meant as an alternative to RPA, NPSS, ESPSS and other regenerative cooling and engine cycle analysis software. The solver produces regenerative cooling and engine cycle results: 



The first iteration of pyskyfire was written as part of the master thesis of Ask Haugerud Hovik, which can be read [here](https://drive.google.com/file/d/1sZJmt-8UWtUChprji67LmnazS3Ei_K3a/view). The motivation to start writing the software came purely from a curiosity standpoint and from an innate wish to spread the understanding of rocket engines and propel us further into the space age. Please use this software responsibly and make sure you, your team memebers and everyone else stay safe in your rocket engine endeavours. Okay. 

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

```
uv sync
```

# Contributions

The pyskyfire project started as my master thesis, but it is now out in the open. I would love for other students and professionals to contribute to pyskyfire. If you are interested in propulsion, is great at coding, and want to use simulation as a path for learning about rocket engines, please reach out. 

# Getting Started
The documentation and examples for pyskyfire is the best place to start. Validation cases can also be used to learn how the program works.
