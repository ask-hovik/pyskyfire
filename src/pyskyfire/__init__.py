"""
pyskyfire - A Python library for rocket engine simulation and design.

Subpackages:
    regen   - Regenerative cooling analysis.
    pump    - Pump design and performance calculations.
    turbine - Turbine design and analysis.

Usage:
    import pyskyfire.regen as regen
    import pyskyfire.pump as pump
    import pyskyfire.turbine as turbine
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


__all__ = ["regen", "pump", "turbine", "common", "skycea"]
