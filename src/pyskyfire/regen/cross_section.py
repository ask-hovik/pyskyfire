from __future__ import annotations
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
from typing import Any, Optional
import numpy as np

from pyskyfire.regen.section_geometry import (
    arc_boundary_metrics, rounded_boundary, squared_boundary,
)

from pyskyfire.regen.design_inputs import (
    DesignContext,
    contour_providers,
    layout_providers,
    make_evaluator,
    section_providers,
)

@dataclass(frozen=True)
class SectionProfiles:
    """All geometry profiles along the cooling-channel centerline.

    Each required array must have consistent length ``N`` (number of stations).

    Attributes
    ----------
    h : np.ndarray
        Input channel height array [m, shape (N,)]. Rounded sections interpret
        this using their ``height_input`` setting: full height by default,
        or the legacy base-height parameter when explicitly selected.
    theta : np.ndarray
        Included wedge angle at each station [rad, shape (N,)].
    t_wall : np.ndarray
        Combined main-wall thickness, excluding hot-side liners
        [m, shape (N,)].
    cross_section_base_point_path : np.ndarray
        Physical section anchors in cylindrical ``(x, r, angle)`` coordinates
        [shape (N, 3)]. Rounded: inner-lobe center. Squared: coolant-floor midpoint.
    sector_radius : np.ndarray, optional
        Rounded sections only: the solved distance from the engine axis to the
        inner lobe center, measured in the section plane [m, shape (N,)]. The
        lobes are tangent to the sector rays at this distance, so it, and not
        the base point's own radius, sets the arc radii. Absent for a section
        that has not been placed, which then falls back to the base radius.
    projection : np.ndarray, optional
        Squared sections only: the 2x2 linear map from section-plane
        ``(u, v)`` to transverse ``(radial, azimuthal)`` coordinates
        [shape (N, 2, 2)]. Section boundaries are built in the transverse frame
        and measured through its inverse.
    azimuthal_pitch : np.ndarray, optional
        Channel-to-channel angular pitch in the transverse plane [rad, shape
        (N,)], before the ``cos(gamma)`` reduction that gives ``theta``.
    winding_rate : np.ndarray, optional
        How fast the channel turns about the engine axis, ``dtheta/dx``
        [rad/m, shape (N,)]. A radius times this rate is the tangent of the
        helix angle seen at that radius, which is what makes a rib's faces
        close in on each other above the channel floor.

    Other Parameters
    ----------------
    contour, gamma, n_channel_positions, n_channels, n_channels_all, \
    t_wall_total, t_hot_layer : optional
        Context needed by design-input functions that are evaluated against a
        resolved section, such as ``rib_width``. They are optional so that
        sections can still be built for pure geometry queries, but a design
        input that requests a quantity backed by a missing field raises
        :class:`~pyskyfire.regen.design_inputs.DesignInputError`.

    Raises
    ------
    ValueError
        If any array length or shape is inconsistent.
    """

    h: np.ndarray            # channel height [N]
    theta: np.ndarray        # included angle [rad, N]
    t_wall: np.ndarray       # wall thickness [N]
    cross_section_base_point_path: np.ndarray
    sector_radius: Optional[np.ndarray] = None
    projection: Optional[np.ndarray] = None
    azimuthal_pitch: Optional[np.ndarray] = None
    winding_rate: Optional[np.ndarray] = None

    # --- optional design-input context -------------------------------------
    contour: Optional[Any] = None
    gamma: Optional[np.ndarray] = None
    n_channel_positions: Optional[int] = None
    n_channels: Optional[int] = None
    n_channels_all: Optional[np.ndarray] = None
    t_wall_total: Optional[np.ndarray] = None
    t_hot_layer: Optional[np.ndarray] = None

    def __post_init__(self):
        N = self.cross_section_base_point_path.shape[0]
        for name in ("h", "theta", "t_wall"):
            v = getattr(self, name)
            if v.shape[0] != N:
                raise ValueError(f"{name} length {v.shape[0]} != N {N}")
        if self.cross_section_base_point_path.shape != (N, 3):
            raise ValueError("centerline must be (N,3)")
        # Not a field, so ``dataclasses.replace`` starts a derived profile off
        # with an empty memo of its own. See :meth:`cached`.
        object.__setattr__(self, "_derived", {})

    def cached(self, key, compute):
        """Return ``compute()`` for ``key``, evaluating it once per profile.

        A profile is immutable, so anything derived from one stays valid for as
        long as it lives. The solver depends on that: it holds a station fixed
        while it iterates that station's temperatures, so it asks the same
        one-station profile for the same section geometry many times over, and
        both the boundary integral and the design-input evaluation behind it
        are far too expensive to repeat.

        ``key`` identifies the computation, and must include whatever object is
        doing it, since one profile can be measured by more than one section.
        """
        cache = self._derived
        try:
            return cache[key]
        except KeyError:
            value = compute()
            cache[key] = value
            return value

    def __getstate__(self):
        """Serialise the fields only; the derived memo is rebuilt on demand."""
        state = self.__dict__.copy()
        state.pop("_derived", None)
        return state

    def __setstate__(self, state):
        # Frozen dataclasses block ``setattr``, and ``__post_init__`` does not
        # run on unpickling, so seed the memo here.
        self.__dict__.update(state)
        self.__dict__.setdefault("_derived", {})

    def design_contexts(self) -> list[DesignContext]:
        """Return one :class:`DesignContext` per station.

        The contexts expose tier A, B and C quantities, so they are what
        section-level design inputs such as ``rib_width`` are evaluated
        against.
        """
        radius = self.cross_section_base_point_path[:, 1]
        contexts = []
        for index, x in enumerate(self.cross_section_base_point_path[:, 0]):
            providers = contour_providers(self.contour, float(x))
            providers.update(
                layout_providers(
                    n_channel_positions=self.n_channel_positions,
                    n_channels=self.n_channels,
                    n_channels_all_circuits=(
                        None
                        if self.n_channels_all is None
                        else float(self.n_channels_all[index])
                    ),
                    t_wall=float(self.t_wall[index]),
                    t_wall_total=(
                        None
                        if self.t_wall_total is None
                        else float(self.t_wall_total[index])
                    ),
                    t_hot_layer=(
                        None
                        if self.t_hot_layer is None
                        else float(self.t_hot_layer[index])
                    ),
                )
            )
            providers.update(
                section_providers(
                    r=(
                        float(radius[index])
                        if self.contour is None
                        else float(self.contour.r(float(self.cross_section_base_point_path[index, 0])))
                    ),
                    r_base=float(radius[index]),
                    h=float(self.h[index]),
                    gamma=(
                        None if self.gamma is None else float(self.gamma[index])
                    ),
                    theta=float(self.theta[index]),
                    t_wall=float(self.t_wall[index]),
                )
            )
            contexts.append(DesignContext(providers))
        return contexts

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

    def boundary_metrics(self, prof: SectionProfiles):
        """Return ``(area, centroid, edge_lengths)`` for this section on ``prof``.

        One numerical contour integral of :meth:`boundary` supplies the coolant
        area, the centroid and every individual wetted edge, so it is computed
        once per profile and shared between the accessors that read it.
        """
        return prof.cached(
            (self, "boundary_metrics"),
            lambda: self._compute_boundary_metrics(prof),
        )

    def _compute_boundary_metrics(self, prof: SectionProfiles):
        """Integrate this section's boundary; override to apply a transform."""
        return arc_boundary_metrics(*self.boundary(prof))

    def full_height(self, prof: SectionProfiles) -> np.ndarray:
        """Return the full internal coolant height at each station [m]."""
        return np.asarray(prof.h, dtype=float).copy()

    def base_height(self, prof: SectionProfiles) -> np.ndarray:
        """Return the base height [m]; ordinary sections use their full height."""
        return self.full_height(prof)

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
    """Rectangular channel separated from its neighbours by constant-width ribs.

    Each channel occupies an angular pitch ``theta`` of the wall. The rib
    between two neighbouring channels has the same thickness ``rib_width`` at
    every radius, from the channel floor to the channel roof, so the channel
    itself widens with radius while the rib does not. Consequently the coolant
    area is the full-pitch sector minus the rib rectangle::

        A = theta * (r_top**2 - r_base**2) / 2 - rib_width * h

    Parameters
    ----------
    rib_width : float or callable
        Rib thickness [m], measured circumferentially in the plane normal to
        the channel axis. May be a constant or a design-input function; see
        :mod:`pyskyfire.regen.design_inputs` for the quantities such a
        function may request. Tier A, B and C quantities are all available
        here.
    n_points : int, optional
        Resolution hint for discretised representations.

    Examples
    --------
    A constant one-millimetre rib::

        CrossSectionSquared(rib_width=1.0e-3)

    A rib that always takes up a tenth of the local channel pitch, which is
    the constant-blockage behaviour::

        CrossSectionSquared(rib_width=lambda pitch_mid: 0.1 * pitch_mid)
    """

    def __init__(self, rib_width, n_points: int = 8):
        super().__init__(n_points=n_points)
        self._rib_width = rib_width
        self._rib_width_evaluator = make_evaluator(
            rib_width,
            slot="rib_width",
            label="CrossSectionSquared.rib_width",
        )

    # ---- rib geometry ------------------------------------------------------
    def rib_width(self, prof: SectionProfiles) -> np.ndarray:
        """Return the rib thickness at every station of ``prof`` [m]."""
        return prof.cached(
            (self, "rib_width"),
            lambda: self._compute_rib_width(prof),
        ).copy()

    def _compute_rib_width(self, prof: SectionProfiles) -> np.ndarray:
        widths = np.array(
            [self._rib_width_evaluator(ctx) for ctx in prof.design_contexts()],
            dtype=float,
        )
        if np.any(~np.isfinite(widths)) or np.any(widths < 0.0):
            index = int(np.flatnonzero(~np.isfinite(widths) | (widths < 0.0))[0])
            raise ValueError(
                "rib_width must be finite and non-negative; station "
                f"{index} gives {widths[index]:.9g} m"
            )

        pitch_base = self._pitch_base(prof)
        too_wide = widths >= pitch_base
        if np.any(too_wide):
            index = int(np.flatnonzero(too_wide)[0])
            raise ValueError(
                "rib_width leaves no room for the channel; station "
                f"{index} has a rib of {widths[index]:.9g} m but a "
                f"channel-floor pitch of only {pitch_base[index]:.9g} m"
            )
        return widths

    def _pitch_base(self, prof: SectionProfiles) -> np.ndarray:
        """Return the channel-to-channel pitch at the channel floor [m]."""
        r_base = prof.cross_section_base_point_path[:, 1]
        return r_base * prof.theta

    def _radii(self, prof: SectionProfiles) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r_base = prof.cross_section_base_point_path[:, 1]
        return r_base - prof.t_wall, r_base, r_base + prof.h

    def equivalent_angular_width(self, prof: SectionProfiles) -> np.ndarray:
        """Return the wedge angle enclosing the same coolant area [rad].

        A constant-width rib does not produce a wedge-shaped channel, so
        consumers that need a scalar angular opening can use this
        area-equivalent angle. Physical construction retains the rib faces.
        """
        _, r_base, r_top = self._radii(prof)
        span = r_top**2 - r_base**2
        equivalent = np.zeros_like(span, dtype=float)
        np.divide(2.0 * self.A_coolant(prof), span, out=equivalent, where=span > 0.0)
        return equivalent

    def open_fraction(self, prof: SectionProfiles) -> np.ndarray:
        """Return the fraction of the angular pitch left open by the ribs [-]."""
        fraction = np.zeros_like(prof.theta, dtype=float)
        np.divide(
            self.equivalent_angular_width(prof),
            prof.theta,
            out=fraction,
            where=prof.theta > 0.0,
        )
        return fraction

    # ---- section analytics -------------------------------------------------
    def boundary(self, prof):
        return squared_boundary(prof, self.rib_width(prof))

    def _compute_boundary_metrics(self, prof):
        transform = None if prof.projection is None else np.linalg.inv(prof.projection)
        return arc_boundary_metrics(*self.boundary(prof), transform=transform)

    def centroid(self, prof):
        """Coolant-area centroid in the normal section plane [u, v]."""
        return self.boundary_metrics(prof)[1].copy()

    def A_coolant(self, prof: SectionProfiles) -> np.ndarray:
        return self.boundary_metrics(prof)[0].copy()

    def P_wetted(self, prof: SectionProfiles) -> np.ndarray:
        return self.boundary_metrics(prof)[2].sum(axis=1)

    def Dh_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute hydraulic diameter [m]."""
        return 4.0 * self.A_coolant(prof) / self.P_wetted(prof)

    def P_thermal(self, prof: SectionProfiles) -> np.ndarray:
        """Return geometric hot-gas thermal-contact perimeter [m].

        The rib footprint is excluded because the heat entering it is carried
        by the fin term in :meth:`R_coolant_per_len` rather than by the
        channel floor.
        """
        r_hot, _, _ = self._radii(prof)
        return np.maximum(r_hot * prof.theta - self.rib_width(prof), 0.0)

    def P_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Return the channel-floor perimeter wetted by coolant [m]."""
        return self.boundary_metrics(prof)[2][:, 0].copy()

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
        - Base perimeter = P_coolant(prof)  (the channel floor)
        - Fin perimeter  = 2*h             (two side walls per channel)
        - Effective perimeter = P_base + eta_f * P_fin
        - R_s = 1 / (h_c * P_eff)

        Because the rib has a constant thickness, its fin-root thickness is
        exactly ``rib_width`` rather than an estimate.
        """
        P_base = self.P_coolant(prof)

        # Two sidewalls of height h contribute as "fin area": per unit channel
        # length, area = 2*h, treated as an equivalent coolant-side perimeter.
        h = prof.h
        P_fin = 2.0 * h

        t_rib = np.maximum(self.rib_width(prof), 1e-9)

        # representative conductivity (allow scalar or array)
        k = np.asarray(k_wall, dtype=float)
        if k.ndim == 0:
            k = np.full_like(h, float(k))

        # Fin efficiency (adiabatic tip) for a rectangular fin:
        # m = sqrt(2*h_c/(k*t_rib)), eta = tanh(m*h)/(m*h)
        m = np.sqrt(2.0 * h_c / (k * t_rib))
        mh = np.maximum(m * h, 1e-12)
        eta = np.tanh(mh) / mh

        P_eff = np.maximum(P_base + eta * P_fin, 1e-30)

        return 1.0 / (h_c * P_eff)


# ===================== Rounded implementation ===============================
class _RoundedGeometryMixin:
    """Shared interpretation of full or legacy base-height inputs."""

    def __init__(self, n_points: int = 16, *, height_input: str = "full"):
        if height_input not in ("full", "base"):
            raise ValueError("height_input must be 'full' or 'base'")
        super().__init__(n_points=n_points)
        self._height_input = height_input

    @property
    def height_input(self) -> str:
        """Interpretation of ``SectionProfiles.h`` (defaults for older pickles)."""
        return getattr(self, "_height_input", "full")

    def _beta_alpha(self, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return the inner (β) and outer (α) complementary angles."""
        beta = (np.pi - theta) * 0.5
        alpha = np.pi - beta
        return beta, alpha

    def minimum_channel_height(self, prof: SectionProfiles) -> np.ndarray:
        """Return the smallest full height at each station [m].

        At zero base height the inner and outer lobes form a circle whose
        radius is ``r*sin(theta/2) - t_wall``. The measurable height is then
        that circle's diameter.
        """
        r = prof.sector_radius
        if r is None:
            r = prof.cross_section_base_point_path[:, 1]
        inner_arc_radius = r * np.sin(prof.theta / 2.0) - prof.t_wall
        return 2.0 * inner_arc_radius

    def base_height(self, prof: SectionProfiles) -> np.ndarray:
        """Return the rounded model's base height for the selected input mode.

        If ``h_b`` is the old/base height and ``H`` is the measurable height,
        the arc radii and straight-side length are::

            q_inner = r sin(theta/2) - t_wall
            q_outer = q_inner + h_b sin(theta/2)
            L_side  = h_b cos(theta/2)

        Hence ``H = q_inner + h_b + q_outer`` and this method evaluates
        the exact inverse. A base height of zero is retained as the valid,
        completely round limiting case.
        """
        minimum = self.minimum_channel_height(prof)
        invalid_radius = ~np.isfinite(minimum) | (minimum <= 0.0)
        if np.any(invalid_radius):
            i = int(np.flatnonzero(invalid_radius)[0])
            raise ValueError(
                "rounded channel inner arc radius must be positive; "
                f"station {i} gives a minimum measurable height of {minimum[i]:.9g} m"
            )

        measured = np.asarray(prof.h, dtype=float)
        if np.any(~np.isfinite(measured)):
            raise ValueError("rounded channel height must be finite")
        if self.height_input == "base":
            if np.any(measured < 0.0):
                raise ValueError("rounded channel base height must be non-negative")
            return measured.copy()

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
        denominator = 1.0 + np.sin(half_theta)
        return np.maximum((measured - minimum) / denominator, 0.0)

    def full_height(self, prof: SectionProfiles) -> np.ndarray:
        """Return full internal height, independent of the input convention.

        Uses the rounded model's existing relation
        ``H = 2*q_inner + h_base*(1 + sin(theta/2))``.
        The base parameter is distinct from the straight-side length
        ``h_base*cos(theta/2)`` and excludes the rounded end lobes.
        """
        base = self.base_height(prof)
        if self.height_input == "full":
            return np.asarray(prof.h, dtype=float).copy()
        half_theta = prof.theta / 2.0
        return self.minimum_channel_height(prof) + base * (
            1.0 + np.sin(half_theta)
        )


class CrossSectionRounded(_RoundedGeometryMixin, ChannelSection):
    """Rounded cooling channel with a selectable input-height convention.

    Parameters
    ----------
    n_points : int, optional
        Resolution hint for discretized geometry. Default is 16.
    height_input : {"full", "base"}, optional
        Interpret the supplied height as full internal coolant height
        (default), or as the legacy base-height parameter, excluding the
        rounded lobes. Reports still default to :meth:`full_height`.

    Examples
    --------
    Use the legacy base-height input::

        CrossSectionRounded(height_input="base")
    """

    def boundary(self, prof):
        return rounded_boundary(prof, self.base_height(prof))

    def centroid(self, prof):
        """Coolant-area centroid relative to the inner-lobe center [u, v]."""
        return self.boundary_metrics(prof)[1].copy()

    def A_coolant(self, prof: SectionProfiles) -> np.ndarray:
        return self.boundary_metrics(prof)[0].copy()

    def Dh_coolant(self, prof: SectionProfiles) -> np.ndarray:
        area, _, lengths = self.boundary_metrics(prof)
        return 4 * area / lengths.sum(axis=1)

    def P_thermal(self, prof: SectionProfiles) -> np.ndarray:
        """Compute thermal-contact perimeter [m]."""
        r = prof.sector_radius
        if r is None:
            r = prof.cross_section_base_point_path[:, 1]
        r1 = r * np.sin(prof.theta/2)
        angle = np.radians(112.0)     # your chosen effective contact angle
        return r1 * angle

    def P_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute coolant-wetted perimeter [m]."""
        r = prof.sector_radius
        if r is None:
            r = prof.cross_section_base_point_path[:, 1]
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
    """Internal rounded channel with the same height options as
    :class:`CrossSectionRounded`.

    ``height_input="full"`` is the default; ``height_input="base"`` accepts
    the legacy base-height parameter. ``n_points`` defaults to 16.
    """

    def boundary(self, prof):
        return rounded_boundary(prof, self.base_height(prof))

    def centroid(self, prof):
        """Coolant-area centroid relative to the inner-lobe center [u, v]."""
        return self.boundary_metrics(prof)[1].copy()

    def A_coolant(self, prof: SectionProfiles) -> np.ndarray:
        return self.boundary_metrics(prof)[0].copy()

    def Dh_coolant(self, prof: SectionProfiles) -> np.ndarray:
        area, _, lengths = self.boundary_metrics(prof)
        return 4 * area / lengths.sum(axis=1)

    def P_thermal(self, prof: SectionProfiles) -> np.ndarray:
        """Compute thermal-contact perimeter [m]."""
        r = prof.sector_radius
        if r is None:
            r = prof.cross_section_base_point_path[:, 1]
        r1 = r * np.sin(prof.theta/2)
        angle = np.radians(112.0)     # your chosen effective contact angle
        return (r1 * angle)*2

    def P_coolant(self, prof: SectionProfiles) -> np.ndarray:
        """Compute coolant-wetted perimeter [m]."""
        r = prof.sector_radius
        if r is None:
            r = prof.cross_section_base_point_path[:, 1]
        r1 = r * np.sin(prof.theta/2)
        angle = np.radians(112.0)
        return ((r1 - prof.t_wall) * angle)*2
