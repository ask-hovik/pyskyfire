"""
pyskyfire - A Python library for rocket engine simulation and design.

Subpackages:
    common  - Functionality combining other subpackages. 
    regen   - Regenerative cooling analysis.
    pump    - Pump design and performance calculations.
    turbine - Turbine design and analysis.
    skycea  - Chemical equilibrium wrapper
    viz     - Visualisation
"""

from importlib import import_module

try:
    from importlib.metadata import PackageNotFoundError, version

    __version__ = version("pyskyfire")
except PackageNotFoundError:  # during editable installs without metadata
    __version__ = "0+unknown"

from . import regen
from . import pump
from . import turbine
from . import common
from . import skycea


__all__ = ["regen", "pump", "turbine", "common", "skycea"]


def __getattr__(name: str):
    """Load optional subpackages only when they are accessed."""

    if name == "viz":
        module = import_module(".viz", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
