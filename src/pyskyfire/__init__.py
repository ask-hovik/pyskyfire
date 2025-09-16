"""
pyskyfire - A Python library for rocket engine simulation and design.

Subpackages:
    regen   - Regenerative cooling analysis.
    pump    - Pump design and performance calculations.
    turbine - Turbine design and analysis.
"""

try:
    from importlib.metadata import version, PackageNotFoundError
    __version__ = version("pyskyfire")
except PackageNotFoundError:  # during editable installs without metadata
    __version__ = "0+unknown"

from . import regen
from . import pump
from . import turbine
from . import common
from . import skycea
from . import viz


__all__ = ["regen", "pump", "turbine", "common", "skycea", "viz"]
