"""Blade surfaces in cylindrical coordinates, shared by impellers and turbines.

A blade is two sides, suction and pressure, each an array of shape
``(span, chord, 3)`` in ``(x, r, theta)``: metres, metres, radians, about the
shaft (+x). The sides meet at the leading and trailing edges, so their first
and last chord stations coincide.
"""
from dataclasses import dataclass

import numpy as np

__all__ = ["BladeSurface", "cylindrical_to_cartesian"]


def frozen_array(value):
    result = np.array(value, dtype=float, copy=True)
    if not np.isfinite(result).all():
        raise ValueError("Geometry must contain only finite coordinates")
    result.setflags(write=False)
    return result


def cylindrical_to_cartesian(cyl):
    a = np.asarray(cyl, dtype=float)
    return np.stack((a[..., 0], a[..., 1]*np.cos(a[..., 2]),
                     a[..., 1]*np.sin(a[..., 2])), axis=-1)


@dataclass(frozen=True)
class BladeSurface:
    """One blade's closed-edge sides, (span, chord, 3), in (x,r,theta)."""
    suction: np.ndarray
    pressure: np.ndarray

    def __post_init__(self):
        for name in ("suction", "pressure"):
            p = frozen_array(getattr(self, name))
            if p.ndim != 3 or p.shape[-1] != 3 or min(p.shape[:2]) < 2 or np.any(p[..., 1] <= 0):
                raise ValueError("Blade surfaces need shape (span,chord,3) and positive radii")
            object.__setattr__(self, name, p)
        if self.suction.shape != self.pressure.shape:
            raise ValueError("Blade surface shapes must match")
        if not np.allclose(self.suction[:, [0, -1]], self.pressure[:, [0, -1]], atol=1e-9, rtol=0):
            raise ValueError("Blade leading and trailing edges must close")
