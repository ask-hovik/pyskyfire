"""Pump geometry contracts. Lengths are metres; cylindrical arrays are (x,r,theta).

The shaft is +x. Positive shaft rotation is right-handed about +x, so a
backward-swept blade has decreasing theta from leading to trailing edge.
"""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator


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
class MeridionalContours:
    """Matched hub/shroud stations, ordered inlet to outlet, in (x,r) metres.

    Each corresponding pair defines a span line. Both curves must have the
    same number of points; independent resampling would lose that pairing.
    """
    hub: np.ndarray
    shroud: np.ndarray

    def __post_init__(self):
        for name in ("hub", "shroud"):
            p = frozen_array(getattr(self, name))
            if p.ndim != 2 or p.shape[1] != 2 or len(p) < 4:
                raise ValueError("Meridional curves need shape (N,2), N >= 4")
            if np.any(p[:, 1] <= 0) or np.any(np.linalg.norm(np.diff(p, axis=0), axis=1) <= 0):
                raise ValueError("Meridional curves need positive radii and distinct stations")
            if np.any(np.diff(p[:, 0]) < -1e-12) or np.any(np.diff(p[:, 1]) < -1e-12):
                raise ValueError("This radial impeller requires nondecreasing axial coordinates and radii")
            object.__setattr__(self, name, p)
        if self.hub.shape != self.shroud.shape:
            raise ValueError("Hub and shroud stations must be paired")
        if np.any(self.shroud[:, 1] < self.hub[:, 1]-1e-12):
            raise ValueError("Shroud radius must not be below hub radius")
        if np.any(self.hub[:, 0] < self.shroud[:, 0]-1e-12):
            raise ValueError("Hub must be axially downstream of the shroud")
        if np.any(np.linalg.norm(self.shroud-self.hub, axis=1) <= 1e-12):
            raise ValueError("Hub and shroud must remain separated, including the outlet")

    @property
    def stations(self):
        s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(self.hub, axis=0), axis=1))]
        return s/s[-1]

    def sample(self, fraction, span_fraction):
        h = PchipInterpolator(self.stations, self.hub)(fraction)
        s = PchipInterpolator(self.stations, self.shroud)(fraction)
        return h + np.asarray(span_fraction)[..., None]*(s-h)

    def wall_tangents(self, fraction):
        """Unit hub and shroud tangents at the ends of one span line."""
        tangents = [PchipInterpolator(self.stations, wall).derivative()(fraction)
                    for wall in (self.hub, self.shroud)]
        return np.array([t/np.linalg.norm(t) for t in tangents])

    def with_clearance(self, clearance):
        gap = np.linalg.norm(self.shroud-self.hub, axis=1)
        if clearance < 0 or np.any(gap <= clearance):
            raise ValueError("tip_clearance must be nonnegative and smaller than the passage height")
        return MeridionalContours(self.hub, self.shroud +
                                  clearance/gap[:, None]*(self.hub-self.shroud))

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
            return MeridionalContours(self.hub, shroud)
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

    def export_npz(self, path):
        """Save sampled geometry without serializing callables or pyTurbo objects."""
        out = Path(path)
        data = dict(hub_curve=self.hub_curve, shroud_curve=self.shroud_curve,
                    suction=self.main_blade.suction, pressure=self.main_blade.pressure,
                    blade_count=self.blade_count, backend_version=self.backend_version,
                    splitter_pitch_fraction=self.splitter_pitch_fraction)
        if self.splitter is not None:
            data.update(splitter_suction=self.splitter.suction, splitter_pressure=self.splitter.pressure)
        with out.open("wb") as f:
            np.savez_compressed(f, **data)
        return out
