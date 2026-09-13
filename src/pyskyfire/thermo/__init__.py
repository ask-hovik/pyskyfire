"""Fluid properties, thermodynamic states and combustion-gas aerothermodynamics.

CEA-backed APIs are loaded on first access so fluid and nozzle-loss utilities
can be imported without loading CEA, including through the package root.
"""

from importlib import import_module as _import_module

from .fluid import Fluid
from .transport import FluidTransport, TransportProperties
from .state import (
    PHASE_LIQUID,
    PHASE_TWO_PHASE,
    PHASE_VAPOUR,
    PHASE_SUPERCRITICAL,
    PHASE_UNKNOWN,
    FluidCapability,
    SaturationState,
    FluidState,
    FluidStateResolver,
    probe_fluid,
)
from .nozzle_loss import (
    EfficiencyTable,
    NozzleLoss,
    DeliveredPerformance,
    RL10A_3_3A_ETA_CF,
    RL10A_3_3A_LOSS,
)

__all__ = [
    "Fluid",
    "FluidTransport",
    "TransportProperties",
    "FluidStateResolver",
    "FluidState",
    "FluidCapability",
    "SaturationState",
    "probe_fluid",
    "PHASE_LIQUID",
    "PHASE_TWO_PHASE",
    "PHASE_VAPOUR",
    "PHASE_SUPERCRITICAL",
    "PHASE_UNKNOWN",
    "Aerothermodynamics",
    "CombustionGasState",
    "verify_state_scaling",
    "EfficiencyTable",
    "NozzleLoss",
    "DeliveredPerformance",
    "RL10A_3_3A_ETA_CF",
    "RL10A_3_3A_LOSS",
]


def __getattr__(name: str):
    """Load combustion-gas APIs only when they are requested."""
    if name in {"Aerothermodynamics", "CombustionGasState", "verify_state_scaling"}:
        module = _import_module(".aerothermodynamics", __name__)
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
