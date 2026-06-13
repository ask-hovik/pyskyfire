"""Preliminary centrifugal impeller geometry generation for pyskyfire.

This module is intentionally a geometry/design-preparation layer, not a pump
meanline/performance solver.  It sizes a first-pass closed radial centrifugal
impeller from a design point, generates meridional streamlines, integrates a
blade camber surface from a blade-angle law, and exposes arrays suitable for
visualisation/meshing/CAD export.

Coordinate convention
---------------------
The internal streamline convention follows the rest of pyskyfire's geometry
style: cylindrical coordinates are stored as ``(x, r, theta)`` where ``x`` is the
axial coordinate [m], ``r`` is radius [m], and ``theta`` is angle [rad].

The implemented correlations are preliminary engineering correlations based on
Gülich-style pump sizing.  They are appropriate for smooth radial, closed,
single-stage pump geometry generation at BEP-like design points; they are not a
validated efficiency/cavitation/stress model.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from math import pi
from pathlib import Path
from typing import Literal
import json

import numpy as np

G0 = 9.80665
TWO_PI = 2.0 * pi

BladeAngleLaw = Literal["linear", "cosine"]


@dataclass(slots=True)
class ImpellerInputs:
    """User-facing design inputs for a first-pass centrifugal impeller.

    Parameters
    ----------
    Q:
        Design-point volumetric flow rate [m^3/s].
    H:
        Design-point pump head [m].
    n:
        Rotational speed [rpm].
    d1:
        Optional impeller inlet / suction-eye outer diameter [m].  If omitted,
        it is estimated from the selected inlet flow coefficient and hub ratio.
    dn:
        Optional hub diameter at impeller inlet [m].  If omitted, ``hub_ratio``
        times ``d1`` is used.
    blade_count:
        Optional number of main blades.  If omitted, a conservative heuristic is
        used.  For real hardware this should normally be selected deliberately.
    beta1_deg, beta2_deg:
        Inlet and outlet blade metal angles [deg], measured relative to the
        circumferential direction in the usual pump convention.  These are used
        only to generate a camber surface here.
    phi1:
        Inlet meridional-flow coefficient used when estimating ``d1``:
        ``phi1 = c_m1/u_1``.
    fd1:
        Safety/enlargement factor applied to the calculated minimum inlet
        diameter to account for blockage and non-uniformity.
    hub_ratio:
        ``dn/d1`` used if ``dn`` is not provided.
    n_streamlines, n_points:
        Resolution of the generated meridional/blade camber surface.
    """

    Q: float
    H: float
    n: float

    d1: float | None = None
    dn: float | None = None
    blade_count: int | None = None

    beta1_deg: float = 18.0
    beta2_deg: float = 25.0
    blade_angle_law: BladeAngleLaw = "cosine"

    phi1: float = 0.20
    fd1: float = 1.05
    hub_ratio: float = 0.30

    fT: float = 1.10
    n_streamlines: int = 5
    n_points: int = 96

    # Shape tuning for the analytical meridional section.  The defaults retain
    # the spirit of the earlier implementation while exposing the knobs.
    z_in_factor: float = 0.12
    z_E_factor: float = 0.75
    outlet_straight_fraction: float = 1.0

    def validate(self) -> None:
        if self.Q <= 0.0:
            raise ValueError("Q must be positive [m^3/s]")
        if self.H <= 0.0:
            raise ValueError("H must be positive [m]")
        if self.n <= 0.0:
            raise ValueError("n must be positive [rpm]")
        if self.d1 is not None and self.d1 <= 0.0:
            raise ValueError("d1 must be positive when supplied")
        if self.dn is not None and self.dn < 0.0:
            raise ValueError("dn must be non-negative when supplied")
        if self.blade_count is not None and self.blade_count < 2:
            raise ValueError("blade_count must be >= 2")
        if not 0.05 <= self.phi1 <= 0.60:
            raise ValueError("phi1 should normally be in the range 0.05..0.60")
        if not 0.0 <= self.hub_ratio < 0.85:
            raise ValueError("hub_ratio must be in [0, 0.85)")
        if self.n_streamlines < 2:
            raise ValueError("n_streamlines must be at least 2")
        if self.n_points < 8:
            raise ValueError("n_points must be at least 8")


@dataclass(slots=True)
class ImpellerDimensions:
    """Main calculated impeller dimensions and BEP-like coefficients."""

    nq: float
    psi: float
    d2: float
    r2: float
    u2: float
    b2: float
    b2_star: float
    d1: float
    r1: float
    dn: float
    rn: float
    u1: float
    blade_count: int
    blade_count_method: str


@dataclass(slots=True)
class ImpellerGeometry:
    """Generated impeller geometry arrays.

    Arrays use cylindrical coordinates ``(x, r, theta)`` [m, m, rad].
    ``meridional`` is ``(n_streamlines, n_points, 2)`` and stores ``(x, r)``.
    ``camber_surface`` is ``(n_streamlines, n_points, 3)`` and stores
    ``(x, r, theta)`` for one blade.
    """

    meridional: np.ndarray
    camber_surface: np.ndarray
    hub_curve: np.ndarray
    shroud_curve: np.ndarray
    beta_deg: np.ndarray
    blade_surfaces: list[np.ndarray] = field(default_factory=list)


class Impeller:
    """First-pass closed radial centrifugal impeller design object.

    The class is deliberately staged internally but convenient externally: on
    initialisation it computes main dimensions and geometry, then stores them in
    explicit dataclasses.  Solvers and visualisation code should consume the
    object through ``dimensions`` and ``geometry`` rather than recomputing.
    """

    def __init__(
        self,
        Q: float,
        H: float,
        n: float,
        *,
        d1: float | None = None,
        dn: float | None = None,
        blade_count: int | None = None,
        beta1_deg: float = 18.0,
        beta2_deg: float = 25.0,
        blade_angle_law: BladeAngleLaw = "cosine",
        phi1: float = 0.20,
        fd1: float = 1.05,
        hub_ratio: float = 0.30,
        fT: float = 1.10,
        n_streamlines: int = 5,
        n_points: int = 96,
    ) -> None:
        self.inputs = ImpellerInputs(
            Q=Q,
            H=H,
            n=n,
            d1=d1,
            dn=dn,
            blade_count=blade_count,
            beta1_deg=beta1_deg,
            beta2_deg=beta2_deg,
            blade_angle_law=blade_angle_law,
            phi1=phi1,
            fd1=fd1,
            hub_ratio=hub_ratio,
            fT=fT,
            n_streamlines=n_streamlines,
            n_points=n_points,
        )
        self.inputs.validate()

        self.dimensions = self._compute_dimensions()
        self.geometry = self._compute_geometry()

        # Backwards-compatible attribute aliases for old scripts.
        self.Q = Q
        self.H = H
        self.n = n
        self.n_q = self.dimensions.nq
        self.psi = self.dimensions.psi
        self.d_2 = self.dimensions.d2
        self.u_2 = self.dimensions.u2
        self.b_2 = self.dimensions.b2
        self.d_1 = self.dimensions.d1
        self.d_n = self.dimensions.dn
        self.z = self.dimensions.blade_count

    @classmethod
    def from_pressure_rise(
        cls,
        *,
        mdot: float,
        rho: float,
        dp: float,
        n: float,
        **kwargs,
    ) -> "Impeller":
        """Construct from mass flow, density and pressure rise.

        This is useful for rocket-engine cycle models where the pump load is
        usually specified as ``dp`` rather than as head.
        """
        if rho <= 0.0:
            raise ValueError("rho must be positive [kg/m^3]")
        if mdot <= 0.0:
            raise ValueError("mdot must be positive [kg/s]")
        if dp <= 0.0:
            raise ValueError("dp must be positive [Pa]")
        Q = mdot / rho
        H = dp / (rho * G0)
        return cls(Q=Q, H=H, n=n, **kwargs)

    def _compute_dimensions(self) -> ImpellerDimensions:
        i = self.inputs
        nq = specific_speed(i.n, i.Q, i.H)
        psi = pressure_coefficient(nq, fT=i.fT)
        d2 = outlet_diameter(psi, i.n, i.H)
        r2 = 0.5 * d2
        u2 = tangential_velocity(i.n, d2)
        b2_star = outlet_width_ratio(nq)
        b2 = b2_star * d2

        d1 = i.d1 if i.d1 is not None else inlet_diameter_from_phi(
            Q=i.Q,
            n=i.n,
            phi1=i.phi1,
            hub_ratio=i.hub_ratio,
            fd1=i.fd1,
        )
        d1 = min(d1, 0.90 * d2)  # keep preliminary radial geometry sane
        r1 = 0.5 * d1

        dn = i.dn if i.dn is not None else i.hub_ratio * d1
        dn = min(dn, 0.95 * d1)
        rn = 0.5 * dn
        u1 = tangential_velocity(i.n, d1)

        if i.blade_count is None:
            z, method = estimate_blade_count(nq)
        else:
            z, method = int(i.blade_count), "user"

        return ImpellerDimensions(
            nq=nq,
            psi=psi,
            d2=d2,
            r2=r2,
            u2=u2,
            b2=b2,
            b2_star=b2_star,
            d1=d1,
            r1=r1,
            dn=dn,
            rn=rn,
            u1=u1,
            blade_count=z,
            blade_count_method=method,
        )

    def _compute_geometry(self) -> ImpellerGeometry:
        mer = meridional_streamlines(
            nq=self.dimensions.nq,
            b2=self.dimensions.b2,
            d1=self.dimensions.d1,
            dn=self.dimensions.dn,
            d2=self.dimensions.d2,
            n_streamlines=self.inputs.n_streamlines,
            n_points=self.inputs.n_points,
        )
        camber, beta_deg = add_blade_theta(
            mer,
            beta1_deg=self.inputs.beta1_deg,
            beta2_deg=self.inputs.beta2_deg,
            law=self.inputs.blade_angle_law,
        )
        surfaces = replicate_blades(camber, self.dimensions.blade_count)
        return ImpellerGeometry(
            meridional=mer,
            camber_surface=camber,
            hub_curve=mer[0],
            shroud_curve=mer[-1],
            beta_deg=beta_deg,
            blade_surfaces=surfaces,
        )

    def to_dict(self) -> dict:
        """Return scalar design data plus geometry arrays as JSON-ready lists."""
        return {
            "inputs": asdict(self.inputs),
            "dimensions": asdict(self.dimensions),
            "geometry": {
                "meridional": self.geometry.meridional.tolist(),
                "camber_surface": self.geometry.camber_surface.tolist(),
                "hub_curve": self.geometry.hub_curve.tolist(),
                "shroud_curve": self.geometry.shroud_curve.tolist(),
                "beta_deg": self.geometry.beta_deg.tolist(),
            },
        }

    def export_json(self, path: str | Path) -> Path:
        """Export dimensions and construction curves for CAD/visualisation."""
        out = Path(path)
        out.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")
        return out

    def __repr__(self) -> str:
        return (
            "Impeller("
            f"Q={self.inputs.Q:.6g} m^3/s, H={self.inputs.H:.6g} m, "
            f"n={self.inputs.n:.6g} rpm, d2={self.dimensions.d2:.6g} m, "
            f"d1={self.dimensions.d1:.6g} m, b2={self.dimensions.b2:.6g} m, "
            f"z={self.dimensions.blade_count})"
        )

    def __str__(self) -> str:
        d = self.dimensions
        return (
            "\nImpeller preliminary geometry:\n"
            f"  Q                         = {self.inputs.Q:.6g} m^3/s\n"
            f"  H                         = {self.inputs.H:.6g} m\n"
            f"  n                         = {self.inputs.n:.6g} rpm\n"
            f"  specific speed n_q        = {d.nq:.6g}\n"
            f"  pressure coefficient psi  = {d.psi:.6g}\n"
            f"  outlet diameter d2        = {d.d2:.6g} m\n"
            f"  outlet width b2           = {d.b2:.6g} m\n"
            f"  inlet diameter d1         = {d.d1:.6g} m\n"
            f"  hub diameter dn           = {d.dn:.6g} m\n"
            f"  tip speed u2              = {d.u2:.6g} m/s\n"
            f"  blade count z             = {d.blade_count} ({d.blade_count_method})\n"
            f"  beta1 -> beta2            = {self.inputs.beta1_deg:.3g} -> {self.inputs.beta2_deg:.3g} deg\n"
        )


# ---------------------------------------------------------------------------
# Main sizing correlations
# ---------------------------------------------------------------------------


def specific_speed(n: float, Q: float, H: float) -> float:
    """Gülich-style pump specific speed ``n_q`` [rpm, m^3/s, m]."""
    return float(n * np.sqrt(Q) / (H ** 0.75))


def pressure_coefficient(nq: float, *, fT: float = 1.10) -> float:
    """Approximate optimum pressure/head coefficient.

    The coefficient convention used here is ``psi = 2 g H / u2^2``, matching
    the outlet-diameter expression below.  ``fT`` is retained as an explicit
    tuning factor rather than buried as a magic number.
    """
    nq_ref = 100.0
    return float(1.21 * fT * np.exp(-0.77 * nq / nq_ref))


def outlet_diameter(psi: float, n: float, H: float) -> float:
    """Outlet diameter from ``psi = 2 g H / u2^2`` and ``u2 = pi d2 n/60``."""
    if psi <= 0.0:
        raise ValueError("psi must be positive")
    return float(60.0 / (pi * n) * np.sqrt(2.0 * G0 * H / psi))


def tangential_velocity(n: float, d: float) -> float:
    """Circumferential speed [m/s] for rpm and diameter [m]."""
    return float(pi * d * n / 60.0)


def outlet_width_ratio(nq: float) -> float:
    """Approximate impeller outlet width ratio ``b2/d2``.

    This preserves the polynomial already used in the prototype code, with the
    reference value made explicit.
    """
    x = nq / 100.0
    return float(0.017 + 0.262 * x - 0.080 * x**2 + 0.0093 * x**3)


def inlet_diameter_from_phi(
    *,
    Q: float,
    n: float,
    phi1: float,
    hub_ratio: float,
    fd1: float = 1.05,
) -> float:
    """Estimate inlet diameter from continuity and an inlet flow coefficient.

    ``Q = c_m1 A1`` with ``c_m1 = phi1 u1`` and
    ``A1 = pi/4 * d1^2 * (1 - hub_ratio^2)``.
    """
    denom = phi1 * pi**2 * n * (1.0 - hub_ratio**2)
    if denom <= 0.0:
        raise ValueError("invalid phi1/hub_ratio/n combination")
    d1_min = (240.0 * Q / denom) ** (1.0 / 3.0)
    return float(fd1 * d1_min)


def estimate_blade_count(nq: float) -> tuple[int, str]:
    """Return a conservative preliminary main-blade count.

    Blade count is usually a design decision, not a reliable scalar correlation.
    The heuristic below exists to keep geometry generation automatic; serious
    designs should pass ``blade_count`` explicitly.
    """
    if nq < 15.0:
        return 7, "heuristic-low-nq"
    if nq < 45.0:
        return 6, "heuristic-radial"
    if nq < 80.0:
        return 5, "heuristic-mixed"
    return 4, "heuristic-high-nq"


# ---------------------------------------------------------------------------
# Meridional section construction
# ---------------------------------------------------------------------------

# Normalized outer/shroud and inner/hub streamline coordinates from the existing
# prototype.  Stored as arrays for vectorized transformations.  These are used to
# produce smooth construction curves; they are not a substitute for validation.
_ZA_STAR = np.array([
    1.0000, 0.9986, 0.9945, 0.9878, 0.9784, 0.9664, 0.9519, 0.9349,
    0.9155, 0.8938, 0.8698, 0.8437, 0.8156, 0.7855, 0.7537, 0.7201,
    0.6850, 0.6484, 0.6106, 0.5716, 0.5317, 0.4910, 0.4496, 0.4079,
    0.3658, 0.3237, 0.2818, 0.2402, 0.1992, 0.1723, 0.1458, 0.1199,
    0.0944, 0.0820, 0.0697, 0.0576, 0.0456, 0.0335, 0.0224, 0.0111,
    0.0000,
])
_RA_STAR = np.array([
    1.0000, 0.9335, 0.8692, 0.8072, 0.7475, 0.6901, 0.6351, 0.5825,
    0.5325, 0.4849, 0.4401, 0.3971, 0.3569, 0.3196, 0.2839, 0.2506,
    0.2206, 0.1925, 0.1667, 0.1431, 0.1215, 0.1025, 0.0852, 0.0700,
    0.0565, 0.0449, 0.0348, 0.0264, 0.0193, 0.0153, 0.0119, 0.0089,
    0.0064, 0.0053, 0.0043, 0.0034, 0.0026, 0.0018, 0.0012, 0.0005,
    0.0000,
])
_ZI_STAR = np.array([
    1.0000, 0.9911, 0.9735, 0.9526, 0.9302, 0.9050, 0.8834, 0.8603,
    0.8378, 0.8108, 0.7863, 0.7614, 0.7362, 0.7129, 0.6872, 0.6581,
    0.6310, 0.6033, 0.5749, 0.5458, 0.5158, 0.4849, 0.4531, 0.4210,
    0.3861, 0.3508, 0.3143, 0.2763, 0.2370, 0.2099, 0.1821, 0.1540,
    0.1244, 0.1095, 0.0942, 0.0792, 0.0637, 0.0481, 0.0323, 0.0162,
    0.0000,
])
_RI_STAR = np.array([
    1.0000, 0.8068, 0.6969, 0.6195, 0.5610, 0.5134, 0.4729, 0.4374,
    0.4056, 0.3767, 0.3503, 0.3253, 0.3021, 0.2803, 0.2597, 0.2404,
    0.2214, 0.2036, 0.1865, 0.1701, 0.1543, 0.1392, 0.1246, 0.1106,
    0.0972, 0.0843, 0.0720, 0.0602, 0.0498, 0.0418, 0.0349, 0.0283,
    0.0220, 0.0190, 0.0159, 0.0131, 0.0107, 0.0076, 0.0050, 0.0024,
    0.0000,
])


def meridional_streamlines(
    *,
    nq: float,
    b2: float,
    d1: float,
    dn: float,
    d2: float,
    n_streamlines: int = 5,
    n_points: int = 96,
) -> np.ndarray:
    """Generate hub-to-shroud meridional streamlines ``(x, r)``.

    Returns
    -------
    np.ndarray
        Shape ``(n_streamlines, n_points, 2)``.  The first streamline is the hub
        side, the last streamline is the shroud side.  Points run from inlet to
        outlet.
    """
    r1 = 0.5 * d1
    rn = 0.5 * dn
    r2 = 0.5 * d2

    # The existing prototype used these Gülich-style dependencies.  Keep them
    # explicit and bounded for sane first-pass geometry over common radial nq.
    nq_eff = max(float(nq), 1.0)
    r_ga = 2.75 * (1.0 / nq_eff) ** 0.16 * r1
    r_ga = float(np.clip(r_ga, 1.05 * r1, 0.95 * r2))
    r_gi = 0.83 * (nq_eff) ** 0.021 * r_ga
    r_gi = float(np.clip(r_gi, max(1.05 * rn, 0.50 * r_ga), 0.98 * r_ga))

    x_in = 0.12 * max(d1 - dn, 1e-9)
    x_E = 0.75 * 0.5 * max(d1 - dn, 1e-9) * nq_eff ** -0.05
    x_iSL = (x_E + b2) * (0.20 + 0.002 * nq_eff)

    x_ga = x_E
    x_gi = max(x_E + b2 - x_iSL, 0.10 * b2)

    # Curved parts.  Prototype arrays run outlet-ish to inlet-ish; we reverse at
    # the end so all streamlines run from inlet to outlet.
    x_shroud = x_in + x_ga * _ZA_STAR
    r_shroud = r1 + (r_ga - r1) * _RA_STAR
    x_hub = x_in + x_iSL + x_gi * _ZI_STAR
    r_hub = rn + (r_gi - rn) * _RI_STAR

    # Add outlet point where both hub and shroud reach r2 with separated axial
    # locations by b2.  This is a construction surface, not a finished casing.
    x_shroud = np.r_[x_in + x_ga + b2, x_shroud]
    r_shroud = np.r_[r2, r_shroud]
    x_hub = np.r_[x_in + x_gi + x_iSL, x_hub]
    r_hub = np.r_[r2, r_hub]

    shroud = np.column_stack([x_shroud, r_shroud])[::-1]
    hub = np.column_stack([x_hub, r_hub])[::-1]

    shroud_i = interpolate_curve_by_arclength(shroud, n_points)
    hub_i = interpolate_curve_by_arclength(hub, n_points)

    curves = []
    for alpha in np.linspace(0.0, 1.0, n_streamlines):
        curves.append((1.0 - alpha) * hub_i + alpha * shroud_i)
    return np.asarray(curves, dtype=float)


def interpolate_curve_by_arclength(points: np.ndarray, n_points: int) -> np.ndarray:
    """Resample a 2D/3D polyline to approximately uniform arc length."""
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[0] < 2:
        raise ValueError("points must have shape (N, D) with N >= 2")
    ds = np.linalg.norm(np.diff(p, axis=0), axis=1)
    s = np.r_[0.0, np.cumsum(ds)]
    if s[-1] <= 0.0:
        return np.repeat(p[:1], n_points, axis=0)
    s_new = np.linspace(0.0, s[-1], n_points)
    out = np.column_stack([np.interp(s_new, s, p[:, j]) for j in range(p.shape[1])])
    return out


# ---------------------------------------------------------------------------
# Blade camber surface construction
# ---------------------------------------------------------------------------


def blade_angle_distribution(
    n_points: int,
    *,
    beta1_deg: float,
    beta2_deg: float,
    law: BladeAngleLaw = "cosine",
) -> np.ndarray:
    """Blade metal angle distribution from inlet to outlet [deg]."""
    y = np.linspace(0.0, 1.0, n_points)
    if law == "linear":
        f = y
    elif law == "cosine":
        # Smooth zero-slope loading at inlet/outlet.
        f = 0.5 * (1.0 - np.cos(pi * y))
    else:
        raise ValueError(f"unknown blade angle law: {law!r}")
    return beta1_deg + (beta2_deg - beta1_deg) * f


def add_blade_theta(
    meridional: np.ndarray,
    *,
    beta1_deg: float,
    beta2_deg: float,
    law: BladeAngleLaw = "cosine",
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate blade wrap angle on each meridional streamline.

    Uses ``dtheta = dm / (r tan(beta))``.  This is a geometric camber-surface
    construction law, not a solved velocity-triangle model.
    """
    mer = np.asarray(meridional, dtype=float)
    if mer.ndim != 3 or mer.shape[2] != 2:
        raise ValueError("meridional must have shape (n_streamlines, n_points, 2)")

    n_streamlines, n_points, _ = mer.shape
    beta_deg = blade_angle_distribution(
        n_points,
        beta1_deg=beta1_deg,
        beta2_deg=beta2_deg,
        law=law,
    )
    beta_rad = np.deg2rad(np.clip(beta_deg, 1.0, 89.0))

    camber = np.zeros((n_streamlines, n_points, 3), dtype=float)
    camber[:, :, :2] = mer

    for i in range(n_streamlines):
        x = mer[i, :, 0]
        r = np.maximum(mer[i, :, 1], 1e-9)
        dm = np.sqrt(np.diff(x) ** 2 + np.diff(r) ** 2)
        r_mid = 0.5 * (r[:-1] + r[1:])
        beta_mid = 0.5 * (beta_rad[:-1] + beta_rad[1:])
        dtheta = dm / (r_mid * np.tan(beta_mid))
        theta = np.r_[0.0, np.cumsum(dtheta)]
        camber[i, :, 2] = theta

    return camber, beta_deg


def replicate_blades(camber_surface: np.ndarray, blade_count: int) -> list[np.ndarray]:
    """Replicate one blade camber surface around the pump axis."""
    blades: list[np.ndarray] = []
    for k in range(int(blade_count)):
        blade = np.array(camber_surface, dtype=float, copy=True)
        blade[:, :, 2] += k * TWO_PI / blade_count
        blades.append(blade)
    return blades


def cylindrical_to_cartesian(cyl: np.ndarray) -> np.ndarray:
    """Convert ``(..., 3)`` array from ``(x, r, theta)`` to ``(x, y, z)``."""
    a = np.asarray(cyl, dtype=float)
    x = a[..., 0]
    r = a[..., 1]
    th = a[..., 2]
    return np.stack([x, r * np.cos(th), r * np.sin(th)], axis=-1)
