"""Duct hydraulics shared by cooling channels, piping and manifolds.

This subpackage is a leaf: it depends only on numpy, scipy and ``fluids``, and
never imports from :mod:`pyskyfire.regen` or :mod:`pyskyfire.common`. That
keeps it importable from anywhere in the library without a cycle, which is the
whole reason it exists as a separate package rather than living in
:mod:`pyskyfire.common`.

Currently it provides the friction physics behind pressure drop. Centerline
geometry, fittings and the arc-length march will follow.
"""

from .friction import *
