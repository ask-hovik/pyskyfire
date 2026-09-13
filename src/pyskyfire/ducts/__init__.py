"""Duct hydraulics shared by cooling channels, piping and manifolds.

This subpackage contains reusable friction physics and explicit pipe geometry.
It does not import from :mod:`pyskyfire.regen` or :mod:`pyskyfire.common`, so
cooling channels and engine-network adapters can both depend on it without a
cycle. Pipe state records use the lower-level :mod:`pyskyfire.thermo` package.
"""

from .friction import *
from .pipe import *
