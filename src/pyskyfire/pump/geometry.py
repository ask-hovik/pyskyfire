"""Pump geometry contracts. Lengths are metres; cylindrical arrays are (x,r,theta).

The shaft is +x, and which way it turns about that axis is named by the
impeller's ``rotation``: counterclockwise means the right-handed sense about
+x, seen with the axis pointing at the viewer. A backward-swept blade has
decreasing theta from leading to trailing edge in that sense, and a clockwise
impeller is its mirror image, built by negating theta.

+x runs the way the flow does, from the inlet plane toward the back of the
wheel, so counterclockwise here means counterclockwise seen from behind the
wheel looking upstream. That anchor is a property of the impeller itself and
does not move when the machine is repackaged.

It is deliberately *not* the pump-industry convention of naming rotation as
seen from the driver end. Which end the driver sits on is a packaging choice --
a turbine can be mounted either side of a wheel -- so a driver-relative name
silently changes meaning when the drive train is rearranged, which is the
failure this vocabulary exists to prevent. Where a driver end is fixed and a
nameplate is wanted, derive that label from this sense and the drive layout
rather than the other way round.

The sense is named about the impeller's own +x rather than about a world axis,
so a placement that reverses +x also reverses the sense the impeller presents
to whatever it discharges into. A volute collecting that discharge is named the
same way about its own axis; see :mod:`pyskyfire.common.rotation`.
"""
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from pyskyfire.common.blade_surface import BladeSurface, cylindrical_to_cartesian, frozen_array
from pyskyfire.common.meridional import MeridionalPassage
from pyskyfire.common.rotation import COUNTERCLOCKWISE, parse_rotation, rotation_sign


@dataclass(frozen=True)
class MeridionalContours(MeridionalPassage):
    """A radial impeller's passage: axial inlet, radial outlet, in (x,r) metres.

    The shared :class:`~pyskyfire.common.MeridionalPassage` with the checks
    only a radial wheel satisfies: both walls rise monotonically in x and r,
    and the hub runs axially downstream of the shroud.
    """

    def __post_init__(self):
        for name in ("hub", "shroud"):
            p = np.asarray(getattr(self, name), dtype=float)
            if p.ndim != 2 or p.shape[-1] != 2 or len(p) < 4:
                raise ValueError("Meridional curves need shape (N,2), N >= 4")
            if np.any(p[:, 1] <= 0) or np.any(np.linalg.norm(np.diff(p, axis=0), axis=1) <= 0):
                raise ValueError("Meridional curves need positive radii and distinct stations")
            if np.any(np.diff(p[:, 0]) < -1e-12) or np.any(np.diff(p[:, 1]) < -1e-12):
                raise ValueError("This radial impeller requires nondecreasing axial coordinates and radii")
        super().__post_init__()
        if np.any(self.hub[:, 0] < self.shroud[:, 0]-1e-12):
            raise ValueError("Hub must be axially downstream of the shroud")

    def with_leading_edge(self, hub_station, tip_station):
        """Re-pair the stations so one span line carries a slanted leading edge.

        A blade section can only begin on a span line, so slanting the leading
        edge is a statement about which hub and shroud stations pair. The hub
        keeps its stations; the shroud is resampled so that station
        ``hub_station`` meets it where ``tip_station`` did. The monotone map
        relaxes back to the original pairing over the station distance it took
        to build up, keeping the skew near the inlet where the leading edge is.
        """
        for name, value in (("hub", hub_station), ("tip", tip_station)):
            if not np.isscalar(value) or not np.isfinite(value):
                raise ValueError(f"Leading edge {name} station must be a finite scalar")
        if hub_station == tip_station:
            return self
        if not 0 < tip_station < hub_station < 1:
            raise ValueError("Leading edge stations must satisfy 0 < tip < hub < 1")
        s = self.stations
        relaxed = hub_station+(hub_station-tip_station)
        nodes = ([0., hub_station, relaxed, 1.], [0., tip_station, relaxed, 1.]) if relaxed < 1 else (
                [0., hub_station, 1.], [0., tip_station, 1.])
        shroud = PchipInterpolator(s, self.shroud)(PchipInterpolator(*nodes)(s))
        try:
            return replace(self, shroud=shroud)
        except ValueError as error:
            raise ValueError(f"A leading edge slanted from hub station {hub_station:g} to tip "
                             f"station {tip_station:g} does not fit this passage: {error}") from None


@dataclass(frozen=True)
class BladeSettings:
    """Blade extent in normalized hub arc length and edge construction controls.

    ``start`` is the leading edge. One station puts it on an existing span line;
    a ``(hub, tip)`` pair slants it, beginning the blade further downstream on
    the hub as a blade development normally is drawn, and the passage is
    re-paired to carry that line. Both stations are always available as
    ``start_hub`` and ``start_tip``.

    Thickness is a body thickness; both edges close to a point. The leading
    edge is rounded by pyTurbo. ``trailing_edge`` selects its closure style.
    ``angle_law`` is used only to estimate total wrap.
    """
    start: float | tuple[float, float] = 0.05
    end: float = 0.98
    trailing_edge: str = "rounded"
    angle_law: str = "cosine"

    def __post_init__(self):
        start = (self.start, self.start) if np.isscalar(self.start) else tuple(self.start)
        if len(start) != 2 or not all(np.isscalar(v) and np.isfinite(v) for v in start):
            raise ValueError("Blade start must be one station or a (hub, tip) pair of stations")
        object.__setattr__(self, "start", (float(start[0]), float(start[1])))
        if not all(0 <= v < self.end <= 1 for v in self.start):
            raise ValueError("Blade extent must satisfy 0 <= start < end <= 1")
        if self.start_tip > self.start_hub:
            raise ValueError("A slanted leading edge starts no earlier on the hub than at the tip")
        if self.trailing_edge not in ("rounded", "cut"):
            raise ValueError("trailing_edge must be 'rounded' or 'cut'")
        if self.angle_law not in ("linear", "cosine"):
            raise ValueError("angle_law must be 'linear' or 'cosine'")

    @property
    def start_hub(self):
        return self.start[0]

    @property
    def start_tip(self):
        return self.start[1]


@dataclass(frozen=True)
class SplitterSettings:
    """One repeated splitter family following the main-blade camber.

    ``start`` is the fraction of the main camber parameter where it begins;
    ``pitch_fraction`` locates it between adjacent main blades.
    """
    start: float = 0.5
    pitch_fraction: float = 0.5
    thickness_scale: float = 1.0

    def __post_init__(self):
        if not 0 < self.start < 1 or not 0 < self.pitch_fraction < 1:
            raise ValueError("Splitter start and pitch_fraction must lie strictly between 0 and 1")
        if not np.isfinite(self.thickness_scale) or self.thickness_scale <= 0:
            raise ValueError("Splitter thickness_scale must be positive")


@dataclass(frozen=True)
class BladeSection:
    span_fraction: float
    beta1_deg: float
    beta2_deg: float
    wrap_deg: float


@dataclass(frozen=True)
class ImpellerBladeDefinition:
    meridional: MeridionalContours
    sections: tuple[BladeSection, ...]
    blade_count: int
    thickness: object
    settings: BladeSettings
    splitter: SplitterSettings | None
    tip_clearance: float
    rotation: str = COUNTERCLOCKWISE

    def __post_init__(self):
        object.__setattr__(self, "rotation", parse_rotation(self.rotation))


@dataclass(frozen=True)
class ImpellerGeometry:
    hub_curve: np.ndarray
    shroud_curve: np.ndarray
    main_blade: BladeSurface
    blade_count: int
    splitter: BladeSurface | None = None
    splitter_pitch_fraction: float = 0.5
    backend_version: str = ""
    thickness_max_relative_error: float = 0.0
    rotation: str = COUNTERCLOCKWISE

    def __post_init__(self):
        object.__setattr__(self, "rotation", parse_rotation(self.rotation))

    @property
    def rotation_sign(self) -> float:
        """``+1`` counterclockwise about ``+x``, ``-1`` clockwise.

        Blade surfaces are already built in this sense. It is the pitch offset
        between blades that still needs it, since stepping round the wheel has
        to step the same way the wheel is mirrored.
        """
        return rotation_sign(self.rotation)

    def export_npz(self, path):
        """Save sampled geometry without serializing callables or pyTurbo objects."""
        out = Path(path)
        data = dict(hub_curve=self.hub_curve, shroud_curve=self.shroud_curve,
                    suction=self.main_blade.suction, pressure=self.main_blade.pressure,
                    blade_count=self.blade_count, backend_version=self.backend_version,
                    splitter_pitch_fraction=self.splitter_pitch_fraction,
                    rotation=self.rotation)
        if self.splitter is not None:
            data.update(splitter_suction=self.splitter.suction, splitter_pressure=self.splitter.pressure)
        with out.open("wb") as f:
            np.savez_compressed(f, **data)
        return out
