![](https://raw.githubusercontent.com/ask-hovik/pyskyfire/main/images/pyskyfire_header.png)

Pyskyfire is a simulation framework for regeneratively cooled, liquid propellant rocket engines.

------------------
[![License](https://img.shields.io/github/license/ask-hovik/pyskyfire.svg)](https://github.com/ask-hovik/pyskyfire/blob/main/LICENSE)

# Description
Pyskyfire is an open-source python package, meant as an alternative to RPA, NPSS, ESPSS and other regenerative cooling and engine cycle analysis software. It is however a work in progress, and no responsibility for the results of this program can be provided.  

The first iteration of pyskyfire was written as part of the master thesis of Ask Haugerud Hovik, which can be read [here](https://drive.google.com/file/d/1sZJmt-8UWtUChprji67LmnazS3Ei_K3a/view). The motivation to start writing the software came purely from a curiosity standpoint and from an innate wish to spread the understanding of rocket engines and propel us further into the space age. Please use this software responsibly and make sure you, your team memebers and everyone else stay safe in your rocket engine endeavours.  

# Program Capabilities
Features of pyskyfire include thrust chamber chemical equilibrium analysis, multi-pass regenerative cooling, thrust chamber contour generation, pump and turbine utilities, and full rocket engine cycle analysis. A full explanation of the capabilities of pyskyfire is outlined in [capabilities](https://ask-hovik.github.io/pyskyfire/explanations/capabilities.html).

# Documentaion
[Documentation](https://ask-hovik.github.io/pyskyfire) for the project is available, and written in a pedagogic style. The package with examples and validation cases is still being developed, so the documentation does not cover everything the package is capable of. Look into the validation cases and advanced examples to see everything the package can do.

# Installation
The package is available on PyPI, and is simply installed with 

```
uv pip install pyskyfire
```

For the bleeding edge version clone the repository and install in an editable environment. 

``` 
uv pip install -e . 
```  

# Contributions

The pyskyfire project started as my master thesis, but it is now out in the open. I would love for other students and professionals to contribute to pyskyfire. If you are interested in propulsion, is great at coding, and want to use simulation as a path for learning about rocket engines, please reach out. 

# Getting Started
The documentation and examples for pyskyfire is the best place to start. Validation cases can also be used to learn how the program works. 