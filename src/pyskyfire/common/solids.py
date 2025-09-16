from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple, Union
import numpy as np
import warnings

ArrayLike = Union[float, int, np.ndarray]

# ---------------
# Property models
# ---------------

# Base Property Model Class
class PropertyModel:
    def __call__(self, T: ArrayLike) -> np.ndarray:
        raise NotImplementedError

# Simple Constant Model
@dataclass
class ConstantModel(PropertyModel):
    value: float
    def __call__(self, T: ArrayLike) -> np.ndarray:
        T = np.asarray(T, dtype=float)
        return np.full_like(T, self.value, dtype=float)

# Ordinary Polynomial Model
@dataclass
class PolynomialModel(PropertyModel):
    """y = c0 + c1*T + c2*T^2 + ...  (no bounds unless provided)"""
    coeffs: Iterable[float]
    Tmin: Optional[float] = None
    Tmax: Optional[float] = None
    enforce_bounds: bool = False  # keep False for sub-models; Piecewise handles policy

    def __call__(self, T: ArrayLike) -> np.ndarray:
        T = np.asarray(T, dtype=float)
        if self.enforce_bounds and (self.Tmin is not None or self.Tmax is not None):
            if self.Tmin is not None and np.any(T < self.Tmin):
                raise ValueError(f"T below Tmin={self.Tmin}")
            if self.Tmax is not None and np.any(T > self.Tmax):
                raise ValueError(f"T above Tmax={self.Tmax}")
        y = np.zeros_like(T)
        for c in reversed(list(self.coeffs)):
            y = y * T + c
        return y

# Log-polynomial model (used in the NIST Cryo database)
@dataclass
class Log10PolynomialModel(PropertyModel):
    """
    y = 10 ** P(log10(T)), P(x) = sum a_i x^i
    Bounds optional; leave unenforced and let Piecewise handle range policy.
    """
    coeffs: Iterable[float]
    Tmin: Optional[float] = None
    Tmax: Optional[float] = None
    enforce_bounds: bool = False  # keep False; outer model manages clipping

    def __call__(self, T: ArrayLike) -> np.ndarray:
        T = np.asarray(T, dtype=float)
        if np.any(T <= 0):
            raise ValueError("Log10PolynomialModel requires T > 0 K.")
        if self.enforce_bounds and (self.Tmin is not None or self.Tmax is not None):
            if self.Tmin is not None and np.any(T < self.Tmin):
                raise ValueError(f"T below Tmin={self.Tmin}")
            if self.Tmax is not None and np.any(T > self.Tmax):
                raise ValueError(f"T above Tmax={self.Tmax}")
        x = np.log10(T)
        P = np.zeros_like(x)
        for a in reversed(list(self.coeffs)):
            P = P * x + a
        return np.power(10.0, P)

# Interpolate Data Model
@dataclass
class TabulatedModel(PropertyModel):
    Ts: np.ndarray
    Ys: np.ndarray
    enforce_bounds: bool = False  # let Piecewise handle warn&clip

    def __post_init__(self):
        self.Ts = np.asarray(self.Ts, dtype=float)
        self.Ys = np.asarray(self.Ys, dtype=float)
        if self.Ts.ndim != 1 or self.Ys.ndim != 1 or len(self.Ts) != len(self.Ys):
            raise ValueError("TabulatedModel: Ts and Ys must be 1D arrays of same length.")
        if not np.all(np.diff(self.Ts) > 0):
            raise ValueError("TabulatedModel: Ts must be strictly increasing.")

    def __call__(self, T: ArrayLike) -> np.ndarray:
        T = np.asarray(T, dtype=float)
        # No range policy here; Piecewise will pre-clip.
        return np.interp(T, self.Ts, self.Ys)

# Gaussian Sum Model 
@dataclass
class SumOfGaussiansModel(PropertyModel):
    params: List[Tuple[float, float, float]]  # list of (b, c, d)
    def __call__(self, T):
        T = np.asarray(T, dtype=float)
        y = np.zeros_like(T)
        for b, c, d in self.params:
            y += b * np.exp(-((T - c) / d)**2)
        return y

# ================================
# Property Model Combination Class
# ================================
@dataclass
class PiecewiseModel(PropertyModel):
    """
    Piecewise wrapper over sub-models defined on [T_lo, T_hi] segments.

    New rules:
      • GAPS: if no segment covers T but T lies between two segments, linearly
        interpolate between the *endpoint values* of those segments.
      • OVERLAPS: if multiple segments cover T, evaluate all and return the average.
      • OUT-OF-RANGE: clip to nearest edge; if 'warn_clip' in range_policy, warn.

    Notes:
      • Sub-models must be callable: y = model(T).
      • Segments are tuples: (T_lo, T_hi, model) with T_lo < T_hi (strict).
    """
    segments: List[Tuple[float, float, PropertyModel]]
    range_policy: str = "warn_clip"  # 'warn_clip' or 'clip' or 'error'
    # 'blend' retained for API compatibility (unused under the new rules)
    blend: float = 0.0

    def __post_init__(self):
        if not self.segments:
            raise ValueError("PiecewiseModel requires at least one segment.")

        # Normalize & sanity check segments (strict bounds)
        norm: List[Tuple[float, float, PropertyModel]] = []
        for (a, b, m) in self.segments:
            a = float(a); b = float(b)
            if not (b > a):
                raise ValueError(f"Segment must have T_hi > T_lo; got [{a}, {b}].")
            norm.append((a, b, m))

        # Sort by start temperature
        norm.sort(key=lambda s: (s[0], s[1]))
        self._segments = norm

        # Precompute useful arrays
        self._starts = np.array([s[0] for s in self._segments], dtype=float)
        self._ends   = np.array([s[1] for s in self._segments], dtype=float)
        self._models = [s[2] for s in self._segments]

        self._Tmin = float(self._starts.min())
        self._Tmax = float(self._ends.max())

        # Cache endpoint values for gap interpolation (evaluated lazily on first use)
        self._edge_cache = {}  # (idx, 'lo'|'hi') -> float

    # ---------- helpers ----------
    def _edge_value(self, idx: int, side: str) -> float:
        """Evaluate model at the exact segment edge (lo/hi) and cache it."""
        key = (idx, side)
        if key in self._edge_cache:
            return self._edge_cache[key]
        T = self._starts[idx] if side == 'lo' else self._ends[idx]
        val = float(np.asarray(self._models[idx](T)))
        self._edge_cache[key] = val
        return val

    def _find_neighbors(self, T: float) -> Tuple[Optional[int], Optional[int]]:
        """
        Return indices (i_left, i_right) of nearest segments such that:
          end[i_left] <= T <= start[i_right], with end[i_left] < start[i_right].
        If not found on one side, that index is None.
        """
        # Left candidate: segment with maximum end <= T
        left_candidates = np.where(self._ends <= T)[0]
        i_left = int(left_candidates.max()) if left_candidates.size else None
        # Right candidate: segment with minimum start >= T
        right_candidates = np.where(self._starts >= T)[0]
        i_right = int(right_candidates.min()) if right_candidates.size else None
        # Ensure they form a real gap around T
        if i_left is not None and i_right is not None:
            if not (self._ends[i_left] < self._starts[i_right]):
                # They overlap or touch → not a gap
                return None, None
        return i_left, i_right

    def _covering_indices(self, T: float) -> np.ndarray:
        """Return indices of all segments that cover T (inclusive bounds)."""
        return np.where((self._starts <= T) & (T <= self._ends))[0]

    def _handle_oob(self, T: float) -> float:
        """Out-of-bounds policy: clip to nearest edge, warn as configured."""
        if 'error' in self.range_policy:
            raise ValueError(f"T={T} K outside piecewise range [{self._Tmin}, {self._Tmax}] K.")
        if 'warn' in self.range_policy:
            warnings.warn(
                f"T={T:.6f} K outside piecewise range "
                f"[{self._Tmin:.6f}, {self._Tmax:.6f}] K; clipping.",
                RuntimeWarning
            )
        # Clip to nearest edge by evaluating the nearest covering endpoint model
        if T < self._Tmin:
            # nearest is the first segment's lo-edge
            return self._edge_value(0, 'lo')
        else:
            # nearest is the last segment's hi-edge
            return self._edge_value(len(self._segments) - 1, 'hi')

    # ---------- main call ----------
    def __call__(self, T: ArrayLike) -> np.ndarray:
        T = np.asarray(T, dtype=float)
        scalar = False
        if T.ndim == 0:
            T = T[None]
            scalar = True

        out = np.empty_like(T, dtype=float)

        for i, Ti in enumerate(T):
            # 1) In-range? (within global envelope)
            if (Ti < self._Tmin) or (Ti > self._Tmax):
                out[i] = self._handle_oob(Ti)
                continue

            # 2) Check overlaps: average all models that cover Ti
            cov = self._covering_indices(Ti)
            if cov.size >= 1:
                vals = []
                for idx in cov:
                    vals.append(float(np.asarray(self._models[idx](Ti))))
                out[i] = float(np.mean(vals))
                continue

            # 3) No covering segment → possible GAP between neighbors → interpolate endpoints
            i_left, i_right = self._find_neighbors(Ti)
            if (i_left is not None) and (i_right is not None):
                Tl = self._ends[i_left]
                Tr = self._starts[i_right]
                # Endpoint values
                yl = self._edge_value(i_left,  'hi')   # value at Tl
                yr = self._edge_value(i_right, 'lo')   # value at Tr

                # Linear interpolation across the gap
                if Tr == Tl:
                    out[i] = 0.5*(yl + yr)  # degenerate, but safe
                else:
                    w = (Ti - Tl) / (Tr - Tl)  # in [0,1]
                    out[i] = (1.0 - w)*yl + w*yr
                continue

            # 4) Fallback (shouldn’t happen): treat as out-of-bounds
            out[i] = self._handle_oob(Ti)

        if scalar:
            return out[0]
        return out    

# =====================
# Material Getter Class
# =====================
@dataclass
class Material:
    name: str
    k: Optional[PropertyModel] = None
    E: Optional[PropertyModel] = None
    alpha: Optional[PropertyModel] = None
    nu: Optional[PropertyModel] = None
    rho: Optional[PropertyModel] = None

    def get_k(self, T: ArrayLike) -> np.ndarray:
        if self.k is None:
            raise AttributeError(f"{self.name}: thermal conductivity not set.")
        return self.k(T)

    def get_E(self, T: ArrayLike) -> np.ndarray:
        if self.E is None:
            raise AttributeError(f"{self.name}: Young's modulus not set.")
        return self.E(T)

    def get_alpha(self, T: ArrayLike) -> np.ndarray:
        if self.alpha is None:
            raise AttributeError(f"{self.name}: thermal expansion not set.")
        return self.alpha(T)

    def get_nu(self, T: ArrayLike) -> np.ndarray:
        if self.nu is None:
            raise AttributeError(f"{self.name}: Poisson's ratio not set.")
        return self.nu(T)

    def get_rho(self, T: ArrayLike) -> np.ndarray:
        if self.rho is None:
            raise AttributeError(f"{self.name}: density not set.")
        return self.rho(T)

# ------------------------------------------
# 304 stainless, 4–1700 K
# ------------------------------------------

# Cryogenic 304 log10 polynomial (your coefficients), bounds provided but not enforced here
log10poly_304_cryo = Log10PolynomialModel(
    coeffs=[-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199],
    Tmin=4.0, Tmax=300.0, enforce_bounds=False
)

# Mills austenitic fit (300–1800 K): k = 9.2 + 0.0175 T - 2e-6 T^2
mills_304_highT = PolynomialModel(coeffs=[9.2, 0.0175, -2e-6], Tmin=300.0, Tmax=1700.0, enforce_bounds=False)

k_304_piecewise = PiecewiseModel(
    segments=[
        (4.0, 300.0, log10poly_304_cryo),
        (300.0, 1700.0, mills_304_highT),
    ],
    blend=50.0,          # smooth across 300 K
    range_policy="warn_clip"
)

StainlessSteel304 = Material(
    name="AISI 304",
    k=k_304_piecewise,
    E=ConstantModel(193e9),
    alpha=ConstantModel(16e-6),
    nu=ConstantModel(0.29),
    rho=ConstantModel(8000.0)
)

# ------------------------------------------
# Inconel 718, 4–922 K
# ------------------------------------------

k718_cryo = Log10PolynomialModel(
    coeffs=[-8.28921, 39.4470, -83.4353, 98.1690, -67.2088, 26.7082, -5.7205, 0.51115, 0.0],
    Tmin=4.0, Tmax=300.0, enforce_bounds=False
)

T_F = np.array([70,200,400,600,800,1000,1200], dtype=float)
k_BTUin = np.array([79,87,100,112,124,136,148], dtype=float)
T718_hi = (T_F - 32.0) * 5.0/9.0 + 273.15    # -> K
k718_hi = 0.1442 * k_BTUin                   # -> W/m·K

k718_table = TabulatedModel(Ts=T718_hi, Ys=k718_hi)
k718 = PiecewiseModel(segments=[(4.0, 300.0, k718_cryo),
                                (min(T718_hi), max(T718_hi), k718_table)],
                                blend=30.0, range_policy="warn_clip")

# Inconel 718 (UNS N07718)
Inconel718 = Material(
    name="Inconel 718",
    k=k718,                          # your piecewise model
    E=ConstantModel(200e9),          # ~200 GPa
    alpha=ConstantModel(13.0e-6),    # ~13.0 µm/m-K
    nu=ConstantModel(0.29),          # Poisson's ratio
    rho=ConstantModel(8190.0)        # kg/m^3 (8.19 g/cc)
)

# ------------------------------------------
# Inconel 625, 116 - 973 K
# ------------------------------------------

T625_C = np.array([-157,-129,-73,-18,0,21,38,93,204,400,500,600,700], float)
k625_W = np.array([7.3,7.4,8.3,9.2,9.2,9.9,10.0,10.7,12.6,15.3,16.9,18.3,19.8], float)
T625 = T625_C + 273.15

k625 = TabulatedModel(Ts=T625, Ys=k625_W)

# Inconel 625 (UNS N06625)
Inconel625 = Material(
    name="Inconel 625",
    k=k625,                          # your tabulated/fit model
    E=ConstantModel(205e9),          # ~205 GPa
    alpha=ConstantModel(12.8e-6),    # ~12.8 µm/m-K
    nu=ConstantModel(0.27),          # Poisson's ratio
    rho=ConstantModel(8440.0)        # kg/m^3 (8.44 g/cc)
)

# --------
# GRCop-42
# --------

# Monokrousos 2024 sum-of-Gaussians (1–1300 K)
k42_mono = SumOfGaussiansModel(params=[
    (1.404e4,  6.796,  4.891),
    (1.302e4, 12.82,   7.552),
    (4.556e3, 24.51,  15.08 ),
    (7.45e2, -3807.0, 5431.0),
])

# ---- GRCop-42: Chen 2023 measured AM correlation (≈295–973 K) ----
# Placeholder: implement with the exact equation from the paper.
# For now, treat as tabulated if you only have digitized points:
#T_chen = np.array([300, 400, 600, 800, 900, 973], float)  # example anchors
#k_chen = np.array([~340, ~350, ~360, ~330, ~310, ~300], float)  # replace with the paper’s values
#k42_chen = TabulatedModel(Ts=T_chen, Ys=k_chen)

# ---- Piecewise: favor measured AM segment where available; else use Monokrousos ----
#k_GRCop42 = PiecewiseModel(
#    segments=[
#        (1.0,     300.0,  k42_mono),
#        (300.0,   973.0,  k42_chen),  # measured correlation once you plug it in
#        (973.0,  1300.0,  k42_mono),
#    ],
#    blend=30.0,
#    range_policy="warn_clip"
#)

GRCop42 = Material(
    name="GRCop-42",
    k=k42_mono,
    E=ConstantModel(130e9),       # ~129.7 GPa typical RT; vendor/AM data
    alpha=ConstantModel(15.5e-6), # representative RT CTE; Chen et al. measure full curve
    nu=ConstantModel(0.34),       # commonly reported/used value
    rho=ConstantModel(8890.0)     # theoretical density; AM parts ~8790 kg/m^3 after HIP
)

ZirconiumOxide = Material(
    name = "Zirconium Oxide",
    k       = ConstantModel(2.2),
    rho     = ConstantModel(5600)
)

TEOS = Material(
    name = "TEOS",
    k = ConstantModel(0.5) # Somewhere in the range of 0.5 - 1?
)

