from __future__ import annotations
from dataclasses import dataclass
from abc import ABC, abstractmethod
import numpy as np

@dataclass(frozen=True)
class SectionProfiles:
    """All geometry profiles along the cooling-channel centerline.

    Each array must have consistent length ``N`` (number of stations).

    Attributes
    ----------
    h : np.ndarray
        Measurable channel height array [m, shape (N,)]. For rounded channels,
        this spans from the bottom of the inner arc to the top of the outer arc.
    theta : np.ndarray
        Included wedge angle at each station [rad, shape (N,)].
    t_wall : np.ndarray
        Wall thickness between coolant and hot-gas side [m, shape (N,)].
    centerline : np.ndarray
        Centerline coordinates ``(x, r, θ)`` or equivalent system [shape (N, 3)].
    local_coords : np.ndarray
        Local orthonormal frames ``(t, n, b)`` at each station [shape (N, 3, 3)].
    blockage_ratio : np.ndarray
        Fraction [0–1] describing coolant blockage (1 = fully blocked).

    Raises
    ------
    ValueError
        If any array length or shape is inconsistent.
    """

    h: np.ndarray            # channel height [N]
    theta: np.ndarray        # included angle [rad, N]
    t_wall: np.ndarray       # wall thickness [N]
    centerline: np.ndarray   # (N,3) [x, r, z] or whatever convention you use

    def __post_init__(self):
        N = self.centerline.shape[0]
        for name in ("h", "theta", "t_wall"):
            v = getattr(self, name)
            if v.shape[0] != N:
                raise ValueError(f"{name} length {v.shape[0]} != N {N}")
        if self.centerline.shape != (N, 3):
            raise ValueError("centerline must be (N,3)")

class ChannelSection(ABC):
    """Abstract base class for cooling-channel cross-section definitions.

    Derived classes must provide analytical expressions for coolant area,
    hydraulic diameter, wetted perimeters, and optionally geometry export.

    Parameters
    ----------
    n_points : int, optional
        Resolution hint for discretized representations. Default is 16.
    """

    def __init__(self, n_points: int = 16):
        self.n_points = int(n_points)

    # ---- Required API (vectorized over stations) -----------------------------
    @abstractmethod
    def A_coolant(self, prof: SectionProfiles) -> np.ndarray: ...
    @abstractmethod
    def Dh_coolant(self, prof: SectionProfiles) -> np.ndarray: ...
    @abstractmethod
    def P_thermal(self, prof: SectionProfiles) -> np.ndarray: ...
    @abstractmethod
    def P_coolant(self, prof: SectionProfiles) -> np.ndarray: ...
    def effective_hot_gas_perimeter(
        self,
        prof: SectionProfiles,
        surface_area_multiplier: np.ndarray,
    ) -> np.ndarray:
        """Return hot-gas perimeter including a circuit-supplied area factor."""
        return self.P_thermal(prof) * np.asarray(
            surface_area_multiplier,
            dtype=float,
        )

    def R_coolant_per_len(
        self,
        prof: SectionProfiles,
        h_c: np.ndarray,
        k_wall: np.ndarray | float,
        ) -> np.ndarray:
        """
        Effective coolant-side thermal resistance per unit axial length [K m / W].

        This is where ribs/fin efficiency/spreading conduction can be embedded.

        Default fallback (no rib model):
            R' = 1 / (h_c * A'_cool)
        where A'_cool is the geometric coolant thermal area per unit axial length.
        """
        raise NotImplementedError
    def R_wall_per_len(
        self,
        prof: SectionProfiles,
        walls: list,                 # or thickness/k arrays
        T_rep: np.ndarray,
        A_hot_per_len: np.ndarray,   # fallback reference
        ) -> np.ndarray:
        """
        Effective wall-stack conduction resistance per unit axial length [K m / W].
        """
        raise NotImplementedError


# ===================== Squared implementation ===============================
class CrossSectionSquared(ChannelSection):
    """Simplified rectangular channel section (wedge-sector approximation)."""

    def __init__(self, blockage_ratio: float, n_points: int = 8):
        super().__init__(n_points=n_points)
        self._blockage_ratio = blockage_ratio
        assert type(blockage_ratio) in (float, int) or callable(blockage_ratio), \
            "'blockage_ratio' input must be a float, int or callable"

    def blockage_ratio(self, x):
        # x can be float or ndarray
        if callable(self._blockage_ratio):
            return np.asarray(self._blockage_ratio(x), dtype=float)
        x_arr = np.asarray(x)
        if x_arr.ndim:
            return np.full_like(x_arr, float(self._blockage_ratio), dtype=float)
        return float(self._blockage_ratio)

    def _theta_real(self, prof: SectionProfiles) -> np.ndarray:
        """Apply blockage ratio to effective included angle."""
        x = prof.centerline[:, 0]
        br = self.blockage_ratio(x)
        return prof.theta * (1.0 - br)

    def A_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute coolant cross-sectional area [m²]."""
        r = prof.centerline[:, 1]
        r_inner = r + prof.t_wall
        r_outer = r + prof.t_wall + prof.h
        th = self._theta_real(prof)
        A_sector = (th/2) * (r_outer**2 - r_inner**2)
        return A_sector

    def Dh_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute hydraulic diameter [m]."""
        A = self.A_coolant(prof)
        r = prof.centerline[:, 1]
        r_inner = r + prof.t_wall
        r_outer = r + prof.t_wall + prof.h
        th = self._theta_real(prof)
        P = r_inner*th + r_outer*th + 2*prof.h
        return 4.0 * A / P
    
    def P_thermal(self, prof: SectionProfiles) -> np.ndarray:
        """Return geometric hot-gas thermal-contact perimeter [m]."""
        r = prof.centerline[:, 1]
        th = self._theta_real(prof)
        return r * th

    
    def P_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Return coolant-wetted perimeter [m]."""
        r = prof.centerline[:, 1]
        th = self._theta_real(prof)
        r_inner = r + prof.t_wall
        return r_inner * th
    
    def R_coolant_per_len(
        self,
        prof: SectionProfiles,
        h_c: np.ndarray,
        k_wall: np.ndarray | float,
    ) -> np.ndarray:
        """
        Coolant-side thermal resistance per unit *channel length* [K m / W],
        including rib sidewalls as fins (first-order model).

        Model:
        - Base perimeter = P_coolant(prof)  (what you already count)
        - Fin perimeter  = 2*h             (two side walls per channel)
        - Effective perimeter = P_base + eta_f * P_fin
        - R_s = 1 / (h_c * P_eff)
        """
        # Base wetted perimeter per channel length
        P_base = self.P_coolant(prof)

        # Two sidewalls of height h contribute as "fin area": per unit channel length, area = 2*h,
        # and we can treat this as an equivalent perimeter term on the coolant side
        h = prof.h
        P_fin = 2.0 * h

        # Rib thickness estimate at the hot wall (circumferential ligament width)
        # Use your existing theta-real logic:
        theta_real = self._theta_real(prof)                 # open sector [rad]
        theta_tot = prof.theta                              # total local sector [rad]
        r_inner = prof.centerline[:, 1]

        # blocked angle -> ligament width
        theta_block = np.maximum(theta_tot - theta_real, 0.0)
        t_rib = np.maximum(r_inner * theta_block, 1e-9)

        # representative conductivity (allow scalar or array)
        k = np.asarray(k_wall, dtype=float)
        if k.ndim == 0:
            k = np.full_like(h, float(k))

        # Fin efficiency (adiabatic tip) for a rectangular fin:
        # m = sqrt(2*h_c/(k*t_rib)), eta = tanh(m*h)/(m*h)
        m = np.sqrt(2.0 * h_c / (k * t_rib))
        mh = np.maximum(m * h, 1e-12)
        eta = np.tanh(mh) / mh

        P_eff = P_base + eta * P_fin
        P_eff = np.maximum(P_eff, 1e-30)

        return 1.0 / (h_c * P_eff)


# ===================== Rounded implementation ===============================
class _RoundedGeometryMixin:
    """Shared conversion from measurable to base height for rounded sections."""

    def _beta_alpha(self, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return the inner (β) and outer (α) complementary angles."""
        beta = (np.pi - theta) * 0.5
        alpha = np.pi - beta
        return beta, alpha

    def minimum_channel_height(self, prof: SectionProfiles) -> np.ndarray:
        """Return the smallest measurable height accepted at each station [m].

        At zero base height the inner and outer lobes form a circle whose
        radius is ``r*sin(theta/2) - t_wall``. The measurable height is then
        that circle's diameter.
        """
        r = prof.centerline[:, 1]
        inner_arc_radius = r * np.sin(prof.theta / 2.0) - prof.t_wall
        return 2.0 * inner_arc_radius

    def base_height(self, prof: SectionProfiles) -> np.ndarray:
        """Convert measurable full height to the rounded model's base height.

        If ``h_b`` is the old/base height and ``H`` is the measurable height,
        the arc radii and straight-side length are::

            q_inner = r sin(theta/2) - t_wall
            q_outer = q_inner + h_b sin(theta/2)
            L_side  = h_b cos(theta/2)

        Hence ``H = q_inner + L_side + q_outer`` and this method evaluates
        the exact inverse. A base height of zero is retained as the valid,
        completely round limiting case.
        """
        minimum = self.minimum_channel_height(prof)
        invalid_radius = minimum <= 0.0
        if np.any(invalid_radius):
            i = int(np.flatnonzero(invalid_radius)[0])
            raise ValueError(
                "rounded channel inner arc radius must be positive; "
                f"station {i} gives a minimum measurable height of {minimum[i]:.9g} m"
            )

        measured = np.asarray(prof.h, dtype=float)
        tolerance = 16.0 * np.finfo(float).eps * np.maximum(minimum, 1.0)
        invalid_height = measured < minimum - tolerance
        if np.any(invalid_height):
            i = int(np.flatnonzero(invalid_height)[0])
            raise ValueError(
                "rounded channel height is below the geometric minimum; "
                f"station {i} has {measured[i]:.9g} m but requires at least "
                f"{minimum[i]:.9g} m"
            )

        half_theta = prof.theta / 2.0
        denominator = np.sin(half_theta) + np.cos(half_theta)
        return np.maximum((measured - minimum) / denominator, 0.0)


class CrossSectionRounded(_RoundedGeometryMixin, ChannelSection):
    """Rounded cooling channel whose input height is its full radial extent."""
    def __init__(self, n_points: int = 16): # TODO: Does n-points have any funcion left? 
        super().__init__(n_points=n_points)

    def A_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute effective coolant cross-sectional area [m²]."""
        r = prof.centerline[:, 1]
        h_base = self.base_height(prof)
        r1 = r * np.sin(prof.theta/2)
        r2 = (r + h_base) * np.sin(prof.theta/2)
        beta, alpha = self._beta_alpha(prof.theta)
        L_side = h_base * np.cos(prof.theta/2)

        A1 = 0.5 * (r1 - prof.t_wall)**2 * 2 * beta
        A2 = 0.5 * (r2 - prof.t_wall)**2 * 2 * alpha
        S1 = (r1 - prof.t_wall) * L_side
        S2 = (r2 - prof.t_wall) * L_side
        return A1 + A2 + S1 + S2

    def Dh_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute hydraulic diameter [m]."""
        A = self.A_coolant(prof)
        r = prof.centerline[:, 1]
        h_base = self.base_height(prof)
        r1 = r * np.sin(prof.theta/2)
        r2 = (r + h_base) * np.sin(prof.theta/2)
        beta, alpha = self._beta_alpha(prof.theta)
        L_side = h_base * np.cos(prof.theta/2)
        arc1 = (r1 - prof.t_wall) * 2 * beta
        arc2 = (r2 - prof.t_wall) * 2 * alpha
        P = arc1 + arc2 + 2 * L_side
        return 4.0 * A / P

    def P_thermal(self, prof: SectionProfiles) -> np.ndarray:
        """Compute thermal-contact perimeter [m]."""
        r = prof.centerline[:, 1]
        r1 = r * np.sin(prof.theta/2)
        angle = np.radians(112.0)     # your chosen effective contact angle
        return r1 * angle

    def P_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute coolant-wetted perimeter [m]."""
        r = prof.centerline[:, 1]
        r1 = r * np.sin(prof.theta/2)
        angle = np.radians(112.0)
        return (r1 - prof.t_wall) * angle
    
    def R_coolant_per_len(
        self,
        prof: SectionProfiles,
        h_c: np.ndarray,
        k_wall: np.ndarray | float,
    ) -> np.ndarray:
        """
        Coolant-side thermal resistance per unit channel length [K m / W]. 
        In this function I included a rib on the backside of the cooling channel, 
        basically the same way as standard rib calculations with rib efficiency. 
        Unsure if this is completely appropriate. But it had only a minor effect
        on the RL10 validation case, which use low conductivity stainless. For 
        thin copper walls adding the rib has a large effect. 
        """
        P_base = self.P_coolant(prof)
        h = self.base_height(prof)

        # Fin "perimeter contribution" from two side walls
        P_fin = 2.0 * h * np.cos(prof.theta/2)

        # Rib thickness is 2x wall thickness
        t_rib = 2.0 * prof.t_wall

        # Conductivity array
        k = np.asarray(k_wall, dtype=float)
        if k.ndim == 0:
            k = np.full_like(h, float(k))

        # Fin efficiency (adiabatic tip rectangular fin)
        m = np.sqrt(2.0 * h_c / (k * t_rib))
        mh = m * h
        eta = np.ones_like(mh)
        np.divide(np.tanh(mh), mh, out=eta, where=mh != 0.0)

        P_eff = P_base + eta * P_fin

        return 1.0 / (h_c * P_eff)

    

# ===================== Rounded internal implementation ===============================
# TODO: not changed from normal rounded!!!!
class CrossSectionRoundedInternal(_RoundedGeometryMixin, ChannelSection):
    """Internal rounded channel whose input height is its full radial extent."""
    def __init__(self, n_points: int = 16): # TODO: Does n-points have any funcion left? 
        super().__init__(n_points=n_points)

    def A_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute effective coolant cross-sectional area [m²]."""
        r = prof.centerline[:, 1]
        h_base = self.base_height(prof)
        r1 = r * np.sin(prof.theta/2)
        r2 = (r + h_base) * np.sin(prof.theta/2)
        beta, alpha = self._beta_alpha(prof.theta)
        L_side = h_base * np.cos(prof.theta/2)

        A1 = 0.5 * (r1 - prof.t_wall)**2 * 2 * beta
        A2 = 0.5 * (r2 - prof.t_wall)**2 * 2 * alpha
        S1 = (r1 - prof.t_wall) * L_side
        S2 = (r2 - prof.t_wall) * L_side
        return A1 + A2 + S1 + S2

    def Dh_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute hydraulic diameter [m]."""
        A = self.A_coolant(prof)
        r = prof.centerline[:, 1]
        h_base = self.base_height(prof)
        r1 = r * np.sin(prof.theta/2)
        r2 = (r + h_base) * np.sin(prof.theta/2)
        beta, alpha = self._beta_alpha(prof.theta)
        L_side = h_base * np.cos(prof.theta/2)
        arc1 = (r1 - prof.t_wall) * 2 * beta
        arc2 = (r2 - prof.t_wall) * 2 * alpha
        P = arc1 + arc2 + 2 * L_side
        return 4.0 * A / P

    def P_thermal(self, prof: SectionProfiles) -> np.ndarray:
        """Compute thermal-contact perimeter [m]."""
        r = prof.centerline[:, 1]
        r1 = r * np.sin(prof.theta/2)
        angle = np.radians(112.0)     # your chosen effective contact angle
        return (r1 * angle)*2

    def P_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute coolant-wetted perimeter [m]."""
        r = prof.centerline[:, 1]
        r1 = r * np.sin(prof.theta/2)
        angle = np.radians(112.0)
        return ((r1 - prof.t_wall) * angle)*2
