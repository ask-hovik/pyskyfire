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
    r"""Abstract base class for temperature-dependent property models.

    All subclasses must implement ``__call__(T) -> np.ndarray``.

    See Also
    --------
    ConstantModel
        Fixed value, independent of temperature.
    PolynomialModel
        Ordinary polynomial :math:`y = \sum_i c_i T^i`.
    Log10PolynomialModel
        Log-polynomial of :math:`T` as used in cryogenic data.
    TabulatedModel
        Linear interpolation from discrete data points.
    SumOfGaussiansModel
        Smooth empirical fit as a sum of Gaussian peaks.
    PiecewiseModel
        Combines multiple sub-models over temperature segments.
    """
    def __call__(self, T: ArrayLike) -> np.ndarray:
        raise NotImplementedError

# Simple Constant Model
@dataclass
class ConstantModel(PropertyModel):
    """Constant property model returning a fixed value.

    Parameters
    ----------
    value : float
        Constant property value.

    Examples
    --------
    >>> ConstantModel(10.0)([1, 2, 3])
    array([10., 10., 10.])
    """
    value: float
    def __call__(self, T: ArrayLike) -> np.ndarray:
        """Return a constant array equal to :attr:`value`.

        Parameters
        ----------
        T : ArrayLike
            Temperature array (ignored).

        Returns
        -------
        np.ndarray
            Array of same shape as ``T`` filled with :attr:`value`.
        """
        T = np.asarray(T, dtype=float)
        return np.full_like(T, self.value, dtype=float)

# Ordinary Polynomial Model
@dataclass
class PolynomialModel(PropertyModel):
    """Ordinary polynomial property model.

    Implements :math:`y = c_0 + c_1 T + c_2 T^2 + ...`.

    Parameters
    ----------
    coeffs : Iterable[float]
        Polynomial coefficients ``[c0, c1, c2, ...]``.
    Tmin, Tmax : float, optional
        Optional temperature bounds for validity.
    enforce_bounds : bool, optional
        If ``True``, out-of-range temperatures raise :class:`ValueError`.
        Default is ``False``.

    Notes
    -----
    When used inside :class:`PiecewiseModel`, range policy is typically handled
    externally.
    """
    coeffs: Iterable[float]
    Tmin: Optional[float] = None
    Tmax: Optional[float] = None
    enforce_bounds: bool = False  # keep False for sub-models; Piecewise handles policy

    def __call__(self, T: ArrayLike) -> np.ndarray:
        """Evaluate the polynomial.

        Parameters
        ----------
        T : ArrayLike
            Temperature array.

        Returns
        -------
        np.ndarray
            Evaluated property values.
        """

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
    r"""Log-polynomial model (NIST cryogenic form).

    Implements :math:`y = 10^{P(\log_{10} T)}`,
    where :math:`P(x) = \sum_i a_i x^i`.

    Parameters
    ----------
    coeffs : Iterable[float]
        Polynomial coefficients ``[a0, a1, a2, ...]``.
    Tmin, Tmax : float, optional
        Optional temperature bounds for validity.
    enforce_bounds : bool, optional
        If ``True``, out-of-range temperatures raise :class:`ValueError`.
        Default is ``False``.

    Notes
    -----
    Useful for fitting data that spans multiple orders of magnitude in
    temperature or property values.
    """

    coeffs: Iterable[float]
    Tmin: Optional[float] = None
    Tmax: Optional[float] = None
    enforce_bounds: bool = False  # keep False; outer model manages clipping

    def __call__(self, T: ArrayLike) -> np.ndarray:
        """Evaluate the log-polynomial model.

        Parameters
        ----------
        T : ArrayLike
            Temperature array (must be positive).

        Returns
        -------
        np.ndarray
            Evaluated property values.

        Raises
        ------
        ValueError
            If any ``T <= 0`` or bounds are violated.
        """

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
    """Property model based on tabulated data.

    Performs 1D linear interpolation between discrete (T, Y) points.

    Parameters
    ----------
    Ts : array_like
        Monotonic array of temperatures [K].
    Ys : array_like
        Property values corresponding to :attr:`Ts`.
    enforce_bounds : bool, optional
        If ``True``, out-of-range temperatures raise :class:`ValueError`.
        Default is ``False``.

    Raises
    ------
    ValueError
        If ``Ts`` and ``Ys`` have mismatched shapes or are not strictly increasing.
    """

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
        """Interpolate the property at temperature ``T``.

        Parameters
        ----------
        T : ArrayLike
            Temperatures [K].

        Returns
        -------
        np.ndarray
            Interpolated property values.
        """
        T = np.asarray(T, dtype=float)
        # No range policy here; Piecewise will pre-clip.
        return np.interp(T, self.Ts, self.Ys)

# Gaussian Sum Model 
@dataclass
class SumOfGaussiansModel(PropertyModel):
    r"""Empirical model as a sum of Gaussian peaks.

    Implements :math:`y(T) = \sum_i b_i \exp(-((T - c_i)/d_i)^2)`.

    Parameters
    ----------
    params : list[tuple[float, float, float]]
        List of triplets ``(b, c, d)`` representing amplitude, center, and width.
    """

    params: List[Tuple[float, float, float]]  # list of (b, c, d)
    def __call__(self, T):
        """Evaluate the Gaussian sum model.

        Parameters
        ----------
        T : ArrayLike
            Temperature array.

        Returns
        -------
        np.ndarray
            Property values.
        """
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
    """Combine multiple sub-models across temperature segments.

    Each segment is a tuple ``(T_lo, T_hi, model)`` where ``model`` is a callable
    that returns property values. The segments must have strictly increasing
    temperature bounds.

    Parameters
    ----------
    segments : list[tuple[float, float, PropertyModel]]
        Segment list defining temperature intervals and their models.
    range_policy : str, optional
        Range-handling mode: ``"warn_clip"`` (default), ``"clip"`` or ``"error"``.
    blend : float, optional
        Retained for API compatibility (unused in the new rules).
    
    Notes
    -----
    The evaluation rules are:

    - **Gaps:** linearly interpolate between the nearest segment endpoints.
    - **Overlaps:** average results from all covering segments.
    - **Out-of-range:** clip to the nearest valid edge (warns if ``"warn_clip"``).

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
        Return indices (i_left, i_right) of nearest segments such that
        ``end[i_left] <= T <= start[i_right]``, with
        ``end[i_left] < start[i_right]``. If not found on one side, that index
        is None.
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
    """Container for material name and temperature-dependent properties.

    Parameters
    ----------
    name : str
        Material name.
    k : PropertyModel, optional
        Thermal conductivity model [W/m·K].
    E : PropertyModel, optional
        Young’s modulus model [Pa].
    alpha : PropertyModel, optional
        Thermal expansion coefficient model [1/K].
    nu : PropertyModel, optional
        Poisson’s ratio model.
    rho : PropertyModel, optional
        Density model [kg/m³].

    Methods
    -------
    get_k(T)
        Evaluate thermal conductivity.
    get_E(T)
        Evaluate Young’s modulus.
    get_alpha(T)
        Evaluate thermal expansion coefficient.
    get_nu(T)
        Evaluate Poisson’s ratio.
    get_rho(T)
        Evaluate density.

    Notes
    -----
    Each getter raises :class:`AttributeError` if the property is not defined.
    """

    name: str
    k: Optional[PropertyModel] = None
    E: Optional[PropertyModel] = None
    alpha: Optional[PropertyModel] = None
    nu: Optional[PropertyModel] = None
    rho: Optional[PropertyModel] = None

    def get_k(self, T: ArrayLike) -> np.ndarray:
        """Evaluate thermal conductivity.

        Parameters
        ----------
        T : ArrayLike
            Temperature array [K].

        Returns
        -------
        np.ndarray
            Thermal conductivity [W/m·K].

        Raises
        ------
        AttributeError
            If no model is defined for :attr:`k`.
        """
        if self.k is None:
            raise AttributeError(f"{self.name}: thermal conductivity not set.")
        return self.k(T)

    def get_E(self, T: ArrayLike) -> np.ndarray:
        """Evaluate Young’s modulus."""
        if self.E is None:
            raise AttributeError(f"{self.name}: Young's modulus not set.")
        return self.E(T)

    def get_alpha(self, T: ArrayLike) -> np.ndarray:
        """Evaluate coefficient of thermal expansion."""
        if self.alpha is None:
            raise AttributeError(f"{self.name}: thermal expansion not set.")
        return self.alpha(T)

    def get_nu(self, T: ArrayLike) -> np.ndarray:
        """Evaluate Poisson’s ratio."""
        if self.nu is None:
            raise AttributeError(f"{self.name}: Poisson's ratio not set.")
        return self.nu(T)

    def get_rho(self, T: ArrayLike) -> np.ndarray:
        """Evaluate density."""
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
# 347 stainless, 22–1278 K (conductivity)
# ------------------------------------------

# Aerojet-General NERVA Materials Development, NASA CR-125889 (1970),
# Appendix A, DRM M-2, page 4 of 7 (PDF page 248):
# https://ntrs.nasa.gov/api/citations/19720012948/downloads/19720012948.pdf
#
# The source table gives temperature in degF and conductivity in
# Btu/(hr ft degF).  Its quoted uncertainty is +/-15% from -420 to -320 degF,
# +/-10% from -320 to 540 degF, and +/-5% above 540 degF.  Preserve the source
# values here and perform the unit conversion explicitly so they remain easy to
# audit against the report.
T_347_NERVA_F = np.array(
    [-420, -410, -360, -320, -260, -200, -160, -60, -10, 40, 140,
     540, 1040, 1540, 1840],
    dtype=float,
)
k_347_NERVA_BTU_HR_FT_F = np.array(
    [1.39, 1.77, 3.50, 4.50, 5.64, 6.44, 6.85, 7.62, 7.93, 8.21,
     8.73, 10.78, 12.82, 14.93, 16.11],
    dtype=float,
)
T_347_NERVA = (T_347_NERVA_F - 32.0) * 5.0 / 9.0 + 273.15
k_347_NERVA = k_347_NERVA_BTU_HR_FT_F * 1.730734666371391

k_347_table = TabulatedModel(Ts=T_347_NERVA, Ys=k_347_NERVA)
k_347 = PiecewiseModel(
    segments=[(T_347_NERVA[0], T_347_NERVA[-1], k_347_table)],
    range_policy="warn_clip",
)

# Density from the same NERVA data release, DRM M-2, page 2 of 7
# (PDF page 246), in lb/in^3.  The source quotes +/-5% uncertainty.
T_347_RHO_NERVA_F = np.array(
    [-459, -279, -99.4, 80.6, 261, 441, 621, 801, 981, 1161, 1341,
     1521, 1701, 1881, 2061, 2241],
    dtype=float,
)
rho_347_NERVA_LB_IN3 = np.array(
    [0.289, 0.288, 0.287, 0.285, 0.284, 0.283, 0.281, 0.280, 0.278,
     0.277, 0.275, 0.273, 0.271, 0.269, 0.267, 0.265],
    dtype=float,
)
T_347_RHO_NERVA = (T_347_RHO_NERVA_F - 32.0) * 5.0 / 9.0 + 273.15
rho_347_NERVA = rho_347_NERVA_LB_IN3 * 27679.904710191

rho_347_table = TabulatedModel(Ts=T_347_RHO_NERVA, Ys=rho_347_NERVA)
rho_347 = PiecewiseModel(
    segments=[(T_347_RHO_NERVA[0], T_347_RHO_NERVA[-1], rho_347_table)],
    range_policy="warn_clip",
)

StainlessSteel347 = Material(
    name="AISI 347",
    k=k_347,
    # Room-temperature modulus and 20–100 degC mean expansion coefficient:
    # https://www.sandmeyersteel.com/wp-content/uploads/Alloy347-SpecSheet.pdf
    E=ConstantModel(193e9),
    alpha=ConstantModel(16.0e-6),
    # Representative room-temperature value for austenitic stainless steel.
    nu=ConstantModel(0.29),
    rho=rho_347,
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
"""Ni-based superalloy with excellent corrosion resistance and high-temp strength.
Typical use: hot-section hardware, manifolds, and fasteners. Density ~8.44 g/cc.
"""

# --------------------------------------------------
# GRCop-42 and GRCop-84, laser powder bed fusion (AM)
# --------------------------------------------------

# Chen et al., "Thermophysical properties of additively manufactured (AM)
# GRCop-42 and GRCop-84", Materials Today Communications 36 (2023) 106665,
# doi:10.1016/j.mtcomm.2023.106665. Laser-flash conductivity of eight GRCop-42
# and two GRCop-84 L-PBF builds, from several vendors and machines, measured
# from 25 to 700 °C. The GRCop-42 builds fall within about ±4 % of their mean;
# the two GRCop-84 builds differ by about 6 %.
#
# Only 298-973 K is measured. Below it, k is held at its 298 K value: a
# Wiedemann-Franz/Matthiessen model fits the measured range to 2 % in two
# different ways that disagree by 25 % at 120 K, so the data cannot say which
# way it goes. Above it, k continues linearly with the slope at 973 K up to
# the solidus.

def _grcop_k(coeffs_per_kK, T_solidus):
    """Conductivity [W/m·K] from a cubic in T/1000 fitted over 298-973 K."""
    c0, c1, c2, c3 = coeffs_per_kK
    measured = PolynomialModel(coeffs=[c0, c1/1e3, c2/1e6, c3/1e9])
    T_lo, T_hi = 298.15, 973.15
    k_lo, k_hi = float(measured(T_lo)), float(measured(T_hi))
    slope = (c1 + 2*c2*T_hi/1e3 + 3*c3*(T_hi/1e3)**2)/1e3
    return PiecewiseModel(
        segments=[
            (10.0, T_lo, ConstantModel(k_lo)),
            (T_lo, T_hi, measured),
            (T_hi, T_solidus, PolynomialModel(coeffs=[k_hi - slope*T_hi, slope])),
        ],
        range_policy="warn_clip",
    )


def _grcop_alpha(strain_coeffs):
    """Instantaneous CTE [1/K], the derivative of a thermal-strain cubic in (T - 273).

    Chen et al. fit strain over 20-1000 °C; their printed CTE equations are
    garbled, so the CTE is taken from the strain fit instead.
    """
    derivative = np.polynomial.Polynomial(strain_coeffs).deriv()
    in_kelvin = derivative(np.polynomial.Polynomial([-273.0, 1.0]))
    return PiecewiseModel(
        segments=[(293.15, 1273.15, PolynomialModel(coeffs=list(in_kelvin.coef)))],
        range_policy="warn_clip",
    )


# GRCop-42: Chen Eq. 3, the mean of all eight builds (the paper prints its
# linear term as "0.97x"). Solidus 1040 °C.
k_GRCop42_AM = _grcop_k([329.14, 0.97, 4.17, -39.86], T_solidus=1313.15)

GRCop42_AM = Material(
    name="GRCop-42 (AM)",
    k=k_GRCop42_AM,
    E=ConstantModel(130e9),       # ~129.7 GPa typical RT; vendor/AM data
    alpha=_grcop_alpha([-5.89e-4, 1.51e-5, 6.28e-9, -1.88e-12]),  # Chen Eq. 9
    nu=ConstantModel(0.34),       # commonly reported/used value
    rho=ConstantModel(8890.0)     # theoretical density; AM parts ~8790 kg/m^3 after HIP
)
"""GRCop-42 built by laser powder bed fusion, from Chen et al. (2023).

Conductivity is measured over 298-973 K, held constant below and extended
linearly above to the solidus. About 85 % of pure copper at room temperature.
"""

# GRCop-84: Chen give no conductivity equation, only the two builds in their
# Fig. 2(b). Those were digitised at the eight test temperatures and their
# mean fitted by a cubic in T/1000 like Eq. 3 (largest residual 0.8 W/m·K).
# Solidus 1080 °C.
T_GRCop84 = np.array([25.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]) + 273.15
k_GRCop84_builds = np.array([
    [308.6, 304.7, 303.2, 300.5, 296.6, 292.3, 285.9, 275.4],   # GRCop-84-1
    [291.3, 288.4, 287.5, 285.2, 281.7, 277.1, 271.3, 262.1],   # GRCop-84-2
])
k_GRCop84_AM = _grcop_k([319.011, -106.155, 174.029, -121.201], T_solidus=1353.15)

GRCop84_AM = Material(
    name="GRCop-84 (AM)",
    k=k_GRCop84_AM,
    alpha=_grcop_alpha([-6.92e-4, 1.59e-5, 1.63e-9, 1.58e-12]),  # Chen Eq. 10
)
"""GRCop-84 built by laser powder bed fusion, from Chen et al. (2023).

Conductivity is the mean of two builds over 298-973 K, held constant below
and extended linearly above to the solidus. Modulus, Poisson's ratio and
density are not given by the source and are left unset.
"""

#: Kept for existing scripts; it is the AM model.
GRCop42 = GRCop42_AM

ZirconiumOxide = Material(
    name = "Zirconium Oxide",
    k       = ConstantModel(2.2),
    rho     = ConstantModel(5600)
)

# ---------------------------------------
# Silica deposit (dense vitreous silica)
# ---------------------------------------

# Carwile and Hoge, "Thermal Conductivity of Vitreous Silica: Selected Values",
# U.S. Army Natick Laboratories, Technical Report 67-7-PR (1966), DTIC AD0636276.
# A critical review of the published data from 30 to 2100 K.
#
# 80-660 K is their Table 1, the recommended master curve, given in
# millical/(cm s K) and converted here with 1 millical/(cm s K) = 0.4184 W/(m K).
# They put it within ±7 % between 150 and 450 K, widening to ±15 % at 660 K.
#
# 850-1960 K is the dashed "without radiation (estimated)" curve of their
# Fig. 3, digitised here. They derive it by holding k/cp at its 473-623 K value
# and scaling by Sosman's heat capacity, and call it "a rough but reasonable
# estimate". The measured data at those temperatures lie higher because heat
# radiates through clear centimetre-thick samples. That radiative part grows
# with sample thickness and does not apply to a layer tens of microns thick.
# 660-850 K is a straight line between the two.
T_vitreous_silica = np.array(
    [
        80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0,
        260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0,
        440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0,
        620.0, 640.0, 660.0,
        850.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1450.0,
        1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 1960.0,
    ],
    dtype=float,
)
k_vitreous_silica_values = 0.4184 * np.array(
    [
        1.31, 1.61, 1.88, 2.12, 2.33, 2.52, 2.69, 2.84, 2.97,
        3.09, 3.19, 3.28, 3.36, 3.43, 3.50, 3.56, 3.61, 3.66,
        3.72, 3.77, 3.83, 3.88, 3.93, 3.99, 4.04, 4.10, 4.15,
        4.20, 4.24, 4.28,
        4.55, 4.59, 4.69, 4.78, 4.83, 4.90, 4.91, 4.97,
        5.09, 5.33, 5.50, 5.67, 5.76, 5.82,
    ],
    dtype=float,
)
k_vitreous_silica = PiecewiseModel(
    segments=[
        (
            T_vitreous_silica[0],
            T_vitreous_silica[-1],
            TabulatedModel(Ts=T_vitreous_silica, Ys=k_vitreous_silica_values),
        )
    ],
    range_policy="warn_clip",
)

SilicaDeposit = Material(
    name="Silica deposit (dense vitreous SiO2)",
    k=k_vitreous_silica,
)
"""Silica laid down on the hot wall by a silicon fuel additive such as TEOS.

Conductivity is that of dense vitreous silica without radiative transfer, from
Carwile and Hoge (1966): 1.37 W/m·K at 300 K, 1.79 W/m·K at 660 K and about
2.1 W/m·K at 1500 K. Below 660 K it is measured; above, it is their estimate.

No measured conductivity of an additive-formed deposit was found. Such a
deposit is likely porous and then conducts less than this. Only the thermal
resistance ``t/k`` of the layer enters the solve, so read its thickness as the
thickness of dense silica with the same resistance, and calibrate it against a
measured drop in heat load. Use it as a ``role="liner"`` wall ahead of the
chamber wall. Density is left unset, since a porous deposit's is not known.
"""


# -------------------------
# Well-annealed pure silver
# -------------------------

# Recommended high-purity solid-silver values from Ho, Powell, and Liley,
# "Thermal Conductivity of the Elements", J. Phys. Chem. Ref. Data 1,
# 279 (1972), doi:10.1063/1.3253100. The source table is in W/(cm K),
# converted here to W/(m K). Values above 1000 K are estimated by the source.
T_pure_silver = np.array(
    [
        20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0,
        80.0, 90.0, 100.0, 123.2, 150.0, 173.2, 200.0, 223.2,
        250.0, 273.2, 298.2, 300.0, 323.2, 350.0, 373.2, 400.0,
        473.2, 500.0, 573.2, 600.0, 673.2, 700.0, 773.2, 800.0,
        873.2, 900.0, 973.2, 1000.0, 1073.2, 1100.0, 1173.2,
        1200.0, 1235.08,
    ],
    dtype=float,
)
k_pure_silver_values = 100.0 * np.array(
    [
        51.0, 29.5, 19.3, 13.7, 10.5, 8.4, 7.0, 5.30, 4.82,
        4.64, 4.51, 4.44, 4.36, 4.32, 4.31, 4.30, 4.30, 4.29,
        4.29, 4.29, 4.29, 4.28, 4.27, 4.26, 4.25, 4.20, 4.19,
        4.14, 4.12, 4.07, 4.04, 3.99, 3.96, 3.90, 3.88, 3.82,
        3.79, 3.73, 3.70, 3.63, 3.61, 3.58,
    ],
    dtype=float,
)
k_pure_silver = PiecewiseModel(
    segments=[
        (
            T_pure_silver[0],
            T_pure_silver[-1],
            TabulatedModel(Ts=T_pure_silver, Ys=k_pure_silver_values),
        )
    ],
    range_policy="warn_clip",
)

PureSilver = Material(
    name="Pure silver",
    k=k_pure_silver,
    # Room-temperature values from NIST Technical Note 1363,
    # "Cryogenic Properties of Silver".
    E=ConstantModel(78.89e9),
    alpha=ConstantModel(18.9e-6),
    nu=ConstantModel(0.37),
    rho=ConstantModel(10492.0),
)









# ---------- Yellow Brass (70/30 ~ C26000 / CuZn30) ----------
# Sources:
# - Cryogenic k(T) curve for Alloy 260 (70–30 brass): CDA Fig. 7. (units converted)  [CIT: turn5view3]
# - CuZn30 electrical: rho(20°C) ≈ 6.2e-8 Ω·m; alpha(0–100°C) ≈ 1.5e-3 1/K.          [CIT: turn7search0, turn7search14]
# - Wiedemann–Franz engineering use at low–mid T.                                       [CIT: turn9search4]

# --- 1) Cryogenic anchors from CDA Fig. 7 (approximate digitization) ---
# Temperatures chosen at common cryo references; values read from the Alloy 260 curve
# and converted to W/m·K. These reproduce the order of magnitude and trend.
T_cold = np.array([20.0, 77.0, 150.0, 200.0, 250.0, 293.15])  # K
k_cold = np.array([55.0, 70.0, 90.0, 100.0, 110.0, 120.0])     # W/m·K  (approx from CDA curve; 120 W/m·K @RT matches datasheets)

k_yb_cold = TabulatedModel(Ts=T_cold, Ys=k_cold)

# --- 2) High-T via Wiedemann–Franz using CuZn30 electrical data ---
T_ref = 293.15  # 20 °C
rho0  = 6.2e-8  # Ω·m at 20 °C (CuZn30)
alpha = 1.5e-3  # 1/K over 0–100 °C (used as engineering approx up to ~800 K)
L0    = 2.44e-8 # W·Ω·K^-2

T_hi  = np.linspace(200.0, 800.0, 121)                      # K
rho_T = rho0 * (1.0 + alpha * (T_hi - T_ref))
k_hi  = L0 * T_hi / rho_T
k_yb_hi = TabulatedModel(Ts=T_hi, Ys=k_hi)

# Piecewise: use the cryo tabulation below ~300 K; blend into WF above ~200 K.
k_YellowBrass_70_30 = PiecewiseModel(
    segments=[
        (min(T_cold), max(T_cold), k_yb_cold),
        (min(T_hi),   max(T_hi),   k_yb_hi),
    ],
    blend=30.0,
    range_policy="warn_clip"
)

YellowBrass_70_30 = Material(
    name="Yellow Brass 70/30 (≈C26000 / CuZn30)",
    k=k_YellowBrass_70_30,
    E=ConstantModel(115e9),        # representative RT modulus
    alpha=ConstantModel(19.7e-6),  # 20–300 °C mean CTE (CuZn30 datasheets)
    nu=ConstantModel(0.34),
    rho=ConstantModel(8530.0)
)

copper = Material(
    name = "copper",
    k = ConstantModel(400) # Somewhere in the range of 0.5 - 1?
)
