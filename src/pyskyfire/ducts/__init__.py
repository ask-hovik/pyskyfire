"""Duct hydraulics shared by cooling channels, piping and manifolds.

This subpackage contains shared path geometry, explicit pipe hydraulics,
annular ducts, orifices and orifice plates, and volute and collector manifolds. It depends only on :mod:`pyskyfire.common` for
planar section primitives and friction correlations and :mod:`pyskyfire.thermo` for fluid state records,
and it does not import from :mod:`pyskyfire.regen` or :mod:`pyskyfire.network`,
so cooling channels, turbopump hardware and engine-network adapters can all
depend on it without a cycle.
"""

from pyskyfire.common.friction import *
from .path import *
from .pipe import *
from .annular import *
from .orifice import *
from .volute_section import *
from .half_moon import *
from .area_law import *
from .port import *
from .exit_cone import *
from .manifold_law import *
from .manifold import *
from .volute import *
