"""
pyskyfire - A Python library for rocket engine simulation and design.

Subpackages:
    common  - Functionality combining other subpackages. 
    regen   - Regenerative cooling analysis.
    pump    - Pump design and performance calculations.
    turbine - Turbine design and analysis.
    thermo  - Fluid properties, state resolution and chemical equilibrium
    ducts   - Duct hydraulics: friction and pressure drop
    viz     - Visualisation
"""

from importlib import import_module

try:
    from importlib.metadata import PackageNotFoundError, version

    __version__ = version("pyskyfire")
except PackageNotFoundError:  # during editable installs without metadata
    __version__ = "0+unknown"

from . import thermo
from . import ducts
from . import regen
from . import pump
from . import turbine
from . import common


__all__ = ["regen", "pump", "turbine", "common", "thermo", "ducts"]


def __getattr__(name: str):
    """Load optional subpackages only when they are accessed."""

    if name == "viz":
        module = import_module(".viz", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
