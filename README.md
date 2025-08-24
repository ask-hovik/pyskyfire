![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/pyskyfire_header.png)

Pyskyfire is a simulation framework for regeneratively cooled, liquid propellant rocket engines.

------------------
[![License](https://img.shields.io/github/license/ask-hovik/pyskyfire.svg)](https://github.com/ask-hovik/pyskyfire/blob/main/LICENSE)

# Description
Pyskyfire is an open-source python package, meant as an alternative to RPA, NPSS, ESPSS and other regenerative cooling and engine cycle analysis software. It is however a work in progress, and no responsibility for the results of this program can be provided.  

The first iteration of pyskyfire was written as part of the master thesis of Ask Haugerud Hovik, which can be read [here](https://drive.google.com/file/d/1sZJmt-8UWtUChprji67LmnazS3Ei_K3a/view). The motivation to start writing the software came purely from a curiosity standpoint and from an innate wish to spread the understanding of rocket engines and propel us further into the space age. Please use this software responsibly and make sure you, your team memebers and everyone else stay safe in your rocket engine endeavours.  

# Program Capabilities
Features of pyskyfire include thrust chamber chemical equilibrium analysis, multi-pass regenerative cooling, thrust chamber contour generation, pump and turbine utilities, and full rocket engine cycle analysis.

## Regenerative Cooling
A cornerstone of pyskyfire is its built-in regenerative cooling solver. The solver supports multi pass regenerative cooling with multiple propellants. Here demonstrated is the one-and-a-half pass cooling arrangement of the venerable RL10 engine. 

**Coolant Pressure**
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/RL10_coolant_static_pressure.png)
**Coolant Temperature**
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/RL10_coolant_temperature.png)
**Wall Heat Flux**
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/RL10_heat_flux.png)
**Wall Temperature**
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/RL10_wall_temperature.png)

In addition, multiple propellants can be used as coolants in the same thrust chamber, and multiple wall layers can be added, simulating such things as thermal barrier coatings. 

**Multi-Pass, Thermal Barrier Coating Example**
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/methane_engine_wall_temperatures.png)
**Estimate of Temperature Profile**
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/temperature_profile_chamber.png)

## Contour Generation
Pyskyfire can generate optimal contours using either a "Rao" nozzle or a conical nozzle. 

**Bell Nozzle Example**
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/methane_engine_contour.png)

## Aerothermodynamic Properties
The aerothermodynamic properties of the hot gas as it moves through the nozzle is predicted by the program, either using Cantera or NASA CEA. 

![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/skycea_comparison_1.png)
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/skycea_comparison_2.png)

## Full Cycle Analysis
By combining built-in modules for pumps and turbines, the whole engine cycle thermodynamic state can be predicted. Here compared with the RL10 cycle. 

![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/fuel_side_cycle_vali_condensed.png)

The whole engine can be represented in an engine network, and through the use of PT-diagrams. 

**Engine Network Example**
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/engine_network_transp.png)
**Methane-Side Path**
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/methane_PT.png)
**Oxygen-Side Path**
![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/oxygen_PT.png)

## Engine Visualisation 
The engine cooling channels can be visualised in 3D, with interlacing between different cooling loops, and inner and outer wall of each cooling channel rendered. Cross sections through any plane of this render can also be made. Here, the RL10 with its semicircular cooling channels is shown

![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/RL10_visualisation_transp.png
)

# Installation
The package is available on PyPI, and is simply installed with 

```
pip install pyskyfire
```

# Contributions

The pyskyfire project started as my master thesis, but it is now out in the open. I would love for other students and professionals to contribute to pyskyfire. If you are interested in propulsion, is great at coding, and want to use simulation as a path for learning about rocket engines, please reach out. 

# Getting Started
The documentation and examples for pyskyfire is still under construction. Currently, the best way to learn how to use the package is to look in the `examples/` and `validation/` folders and study the different cases. In particular, `validation/regen_vali/RL10_simulation/RL10_simulation_cantera.py` is a good place to start. 