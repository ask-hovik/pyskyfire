"""
enginepy - A Python library for rocket engine simulation and design.

Modules:
    contour         - Contour generation and nozzle shape calculations.
    engine          - Engine performance and combustion reaction classes.
    materials       - Material properties and transport classes.
    physics         - Thermodynamic and heat transfer calculations.
    plot_engine     - Visualization tools for engine geometry.
    solver          - Solvers for system equations.
    thrust_chamber  - Thrust chamber and cooling channel geometry.
"""

from .contour import *
from .physics import *
from .solver import *
from .thrust_chamber import *
from .plot import *
from .solver import *
from .cross_section import *
from .channel_height import *
from .contour_2 import *

