"""Explicit RL10 plumbing, one module per physical component.

The cycle in ``sizer_sim.py`` currently represents every duct as a
``SimpleDuctBlock`` on a fixed pressure ratio. Those ratios are calibrated to
close Binder's station balance, so they cannot be replaced piecemeal: a leg
only becomes physical once every component along it is modelled, otherwise a
partial model displaces a whole-leg ratio and the balance moves for no gain.

This package is where those components accumulate. Each file holds the source
data for one part -- the fabrication drawing's bending table, bore, and
roughness -- and a ``build_pipe()`` returning a ``psf.ducts.Pipe``. Nothing
here is wired into the engine network. When a leg's components are all present,
they can be chained as ``psf.common.PipeBlock``s and that leg's
``eta_`` parameter retired.

Components so far
-----------------
``turbine_discharge_duct``
    Turbine housing discharge to the main fuel shut-off valve, Figure G1. Part
    of the ``duct_turbine_injector`` leg, which also needs the shut-off valve
    and the injector inlet before the ratio can go.
"""

from . import bending_table, turbine_discharge_duct
from .bending_table import (
    INCH,
    dms,
    pipe_from_bending_table,
    route_from_bending_table,
    tangent_distances,
)

__all__ = [
    "INCH",
    "bending_table",
    "dms",
    "pipe_from_bending_table",
    "route_from_bending_table",
    "tangent_distances",
    "turbine_discharge_duct",
]
