from __future__ import annotations
from dataclasses import dataclass, replace
from typing import Any, Optional
import numpy as np

from pyskyfire.common.section_geometry import (
    SectionArc, SectionLine, arc_boundary_metrics, rounded_boundary, squared_boundary,
    rounded_width_boundary, boundary_surfaces, surface_lengths,
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
        resolved section, such as ``sidewall_half_thickness``. They are optional so that
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
    floor_width: Optional[np.ndarray] = None
    roof_width: Optional[np.ndarray] = None
    section_radius: Optional[np.ndarray] = None

    @classmethod
    def from_dimensions(cls, *, height, wall_thickness=0.0, floor_width=None,
                        roof_width=None, radius=None, angle=None, stations=0.0):
        """Build local section profiles without a chamber or placed path.

        Supply both internal endpoint widths, or a finite radius and full
        occupied sector angle (radians). With widths, radius only controls
        squared-arc curvature; omit it for straight ends. Dimensions broadcast
        over stations. ``stations`` is the local path coordinate in metres,
        exposed as ``x`` to design-input callables; it is not a placement.
        """
        direct = floor_width is not None or roof_width is not None
        if direct:
            if floor_width is None or roof_width is None or angle is not None:
                raise ValueError("Supply both floor_width and roof_width, without angle")
        elif radius is None or angle is None:
            raise ValueError("Supply both widths or radius and angle")
        values = [height, wall_thickness, stations,
                  np.inf if radius is None else radius,
                  0.0 if angle is None else angle]
        if direct:
            values += [floor_width, roof_width]
        arrays = [np.atleast_1d(a).astype(float) for a in np.broadcast_arrays(*values)]
        if any(a.ndim != 1 for a in arrays):
            raise ValueError("Section dimensions must be scalars or one-dimensional arrays")
        h, wall, s, r, theta = arrays[:5]
        if np.any(~np.isfinite(h + wall + s)) or np.any(h < 0) or np.any(wall < 0):
            raise ValueError("Height and wall thickness must be finite and nonnegative; stations finite")
        if np.any(np.isnan(r)) or np.any(r <= 0):
            raise ValueError("Radius must be positive or infinity")
        if not direct and (np.any(~np.isfinite(r + theta)) or np.any((theta <= 0) | (theta >= np.pi))):
            raise ValueError("Sector inputs require finite radius and 0 < angle < pi; use widths for parallel walls")
        if direct and (np.any(~np.isfinite(arrays[5] + arrays[6])) or
                       np.any(arrays[5] <= 0) or np.any(arrays[6] <= 0)):
            raise ValueError("Floor and roof widths must be finite and positive")
        return cls(h=h, t_wall=wall, theta=theta,
                   cross_section_base_point_path=np.column_stack((s, np.where(np.isfinite(r), r, 0), np.zeros_like(s))),
                   section_radius=r, floor_width=arrays[5] if direct else None,
                   roof_width=arrays[6] if direct else None)

    def __post_init__(self):
        N = self.cross_section_base_point_path.shape[0]
        for name in ("h", "theta", "t_wall"):
            v = getattr(self, name)
            if v.shape[0] != N:
                raise ValueError(f"{name} length {v.shape[0]} != N {N}")
        if self.cross_section_base_point_path.shape != (N, 3):
            raise ValueError("centerline must be (N,3)")
        if np.any(~np.isfinite(self.t_wall)) or np.any(self.t_wall < 0):
            raise ValueError("Wall thickness must be finite and nonnegative")
        if (self.floor_width is None) != (self.roof_width is None):
            raise ValueError("Supply both floor_width and roof_width")
        for name in ("floor_width", "roof_width", "section_radius"):
            value = getattr(self, name)
            if value is not None and np.shape(value) != (N,):
                raise ValueError(f"{name} must have shape (N,)")
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
        section-level design inputs such as ``sidewall_half_thickness`` are evaluated
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

class ChannelSection:
    """A section with floor, roof, and two sidewalls.

    exposure maps surface names to centred fractions in [0, 1], optionally
    design-input callables. Unspecified surfaces are inactive; the default is
    the complete floor. exposure_angles selects rounded-lobe spans in radians
    instead. A surface cannot specify both forms.

    enable_fin applies one symmetric rib model to both sidewalls, incompatible
    with direct sidewall exposure. CoolingCircuit supplies its own fin switch;
    independent sections use this switch. All surfaces count toward Dh.
    """

    surface_names = ("floor", "side_right", "roof", "side_left")

    def __init__(self, n_points=16, *, exposure=None, exposure_angles=None, enable_fin=False):
        self.n_points = int(n_points)
        if self.n_points < 2:
            raise ValueError("n_points must be at least two")
        self.enable_fin = bool(enable_fin)
        angles = dict(exposure_angles or {})
        fractions = dict(exposure) if exposure is not None else ({} if angles else {"floor": 1.0})
        unknown = (fractions.keys() | angles.keys()) - set(self.surface_names)
        if unknown:
            raise ValueError(f"Unknown thermal surfaces: {sorted(unknown)}")
        if fractions.keys() & angles.keys():
            raise ValueError("A surface cannot specify both exposure fraction and angle")
        if angles and not isinstance(self, CrossSectionRounded):
            raise ValueError("Exposure angles are only supported for rounded floor and roof arcs")
        if set(angles) - {"floor", "roof"}:
            raise ValueError("Exposure angles are only supported for floor and roof")
        self._exposure = fractions
        self._exposure_angles = angles
        self._surface_evaluators = {
            name: make_evaluator(value, slot="surface_exposure", label=f"{type(self).__name__}.{name}.exposure")
            for name, value in (fractions | angles).items()
        }
        for name, value in fractions.items():
            if not callable(value) and (not np.isfinite(value) or not 0 <= value <= 1):
                raise ValueError(f"{name} exposure must be a fraction in [0, 1]")
        for name, value in angles.items():
            if not callable(value) and (not np.isfinite(value) or not 0 <= value <= 2*np.pi):
                raise ValueError(f"{name} exposure angle must be in [0, 2*pi] radians")
        if self.enable_fin:
            self.validate_fin()

    def validate_fin(self, prof=None):
        """Reject direct sidewall transfer combined with the symmetric fin model."""
        for name in ("side_left", "side_right"):
            value = self._exposure.get(name, 0.0)
            if callable(value):
                if prof is None:
                    continue
                value = self.exposure_fractions(prof)[name]
            if np.any(np.asarray(value) > 0):
                raise ValueError("Sidewall exposure and the symmetric fin/rib model are mutually exclusive")

    def __getstate__(self):
        return {key: value for key, value in self.__dict__.items()
                if key != "_surface_evaluators" and not key.endswith("_evaluator")}

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._exposure = state.get("_exposure", {"floor": 1.0})
        self._exposure_angles = state.get("_exposure_angles", {})
        self._surface_evaluators = {
            name: make_evaluator(value, slot="surface_exposure", label=f"{type(self).__name__}.{name}.exposure")
            for name, value in (self._exposure | self._exposure_angles).items()
        }
        if "_sidewall_value" in state:
            self._sidewall_evaluator = make_evaluator(
                self._sidewall_value, slot="sidewall_half_thickness", label="CrossSectionSquared.sidewall_half_thickness")
            self._roof_evaluator = None if self._roof_value is None else make_evaluator(
                self._roof_value, slot="section_thickness", label="CrossSectionSquared.roof_thickness")

    def _evaluate(self, evaluator, prof):
        constant = evaluator.constant_value
        if constant is not None:
            return np.full_like(prof.h, constant, dtype=float)
        return np.array([evaluator(ctx) for ctx in prof.design_contexts()], dtype=float)

    def _transform(self, prof):
        return None

    def surfaces(self, prof, *, exterior=False):
        """Named directed inner or outer curves, including both straight sides."""
        ends = self.outer_boundary(prof) if exterior else self.boundary(prof)
        return dict(zip(self.surface_names, boundary_surfaces(*ends)))

    def boundary_metrics(self, prof):
        return prof.cached((self, "boundary_metrics"), lambda: self._compute_boundary_metrics(prof))

    def _compute_boundary_metrics(self, prof):
        return arc_boundary_metrics(*self.boundary(prof), transform=self._transform(prof))

    def exposure_fractions(self, prof):
        """Resolved centred fractions, with angular requests converted per lobe."""
        def compute():
            surfaces = self.surfaces(prof)
            result = {}
            for name in self.surface_names:
                if name not in self._surface_evaluators:
                    result[name] = np.zeros_like(prof.h)
                    continue
                values = self._evaluate(self._surface_evaluators[name], prof)
                if name in self._exposure_angles:
                    curve = surfaces[name]
                    span = np.abs(curve.end_angles-curve.start_angles)
                    if np.any(values > span + 1e-12):
                        raise ValueError(f"{name} exposure angle exceeds its available arc")
                    values = values/span
                if np.any(~np.isfinite(values)) or np.any((values < 0) | (values > 1+1e-12)):
                    raise ValueError(f"{name} exposure must be a fraction in [0, 1]")
                result[name] = np.clip(values, 0, 1)
            return result
        return {name: value.copy() for name, value in
                prof.cached((self, "exposure"), compute).items()}

    def thermal_perimeters(self, prof, *, exterior=False):
        """Active length of each named surface [m], excluding fin enhancement."""
        key = (self, "thermal_perimeters", exterior)
        def compute():
            fractions = self.exposure_fractions(prof)
            if not exterior:
                lengths = self.boundary_metrics(prof)[2]
                if all(np.all((f == 0) | (f == 1)) for f in fractions.values()):
                    return {name: lengths[:, i]*fractions[name]
                            for i, name in enumerate(self.surface_names)}
            return {name: surface_lengths(curve, fractions[name], self._transform(prof))
                    for name, curve in self.surfaces(prof, exterior=exterior).items()}
        return {name: value.copy() for name, value in prof.cached(key, compute).items()}

    def widths(self, prof):
        """Internal endpoint-to-endpoint floor and roof widths [m]."""
        result = []
        for curve in self.boundary(prof):
            vector = curve.evaluate([1])[:, 0]-curve.evaluate([0])[:, 0]
            transform = self._transform(prof)
            if transform is not None:
                vector = np.einsum("nij,nj->ni", transform, vector)
            result.append(np.linalg.norm(vector, axis=1))
        return tuple(result)

    def full_height(self, prof):
        return np.asarray(prof.h, float).copy()

    def base_height(self, prof):
        return self.full_height(prof)

    def centroid(self, prof):
        if type(self).boundary is ChannelSection.boundary:
            return np.column_stack((prof.h/2, np.zeros_like(prof.h)))
        return self.boundary_metrics(prof)[1].copy()

    def A_coolant(self, prof):
        return self.boundary_metrics(prof)[0].copy()

    def P_wetted(self, prof):
        return self.boundary_metrics(prof)[2].sum(axis=1)

    def Dh_coolant(self, prof):
        return 4*self.A_coolant(prof)/self.P_wetted(prof)

    def A_exterior(self, prof):
        """Area enclosed by the outside of the tube, including the passage [m²]."""
        return arc_boundary_metrics(*self.outer_boundary(prof), transform=self._transform(prof))[0]

    def P_coolant(self, prof):
        return sum(self.thermal_perimeters(prof).values())

    def P_thermal(self, prof):
        return sum(self.thermal_perimeters(prof, exterior=True).values())

    def effective_hot_gas_perimeter(self, prof, surface_area_multiplier):
        return self.P_thermal(prof)*np.asarray(surface_area_multiplier, float)

    def R_coolant_per_len(self, prof, h_c, k_wall, *, enable_fin=None):
        """Convection resistance per path length [K m/W], optionally with ribs."""
        fin = self.enable_fin if enable_fin is None else bool(enable_fin)
        h_c = np.broadcast_to(np.asarray(h_c, float), prof.h.shape)
        if np.any(~np.isfinite(h_c)) or np.any(h_c <= 0):
            raise ValueError("Heat-transfer coefficient must be finite and positive")
        perimeter = self.P_coolant(prof)
        if fin:
            self.validate_fin(prof)
            length, thickness, fin_perimeter = self._fin_dimensions(prof)
            k = np.broadcast_to(np.asarray(k_wall, float), prof.h.shape)
            if np.any(~np.isfinite(k)) or np.any(k <= 0):
                raise ValueError("Fin conductivity must be finite and positive")
            eta = np.zeros_like(length)
            active = thickness > 0
            mh = np.sqrt(2*h_c[active]/(k[active]*thickness[active]))*length[active]
            efficiency = np.ones_like(mh)
            np.divide(np.tanh(mh), mh, out=efficiency, where=mh != 0)
            eta[active] = efficiency
            perimeter = perimeter + eta*fin_perimeter
        resistance = np.full_like(perimeter, np.inf)
        np.divide(1.0, h_c*perimeter, out=resistance, where=perimeter > 0)
        return resistance

    def boundary(self, prof):
        """Default construction for analytical-only custom channel sections."""
        return squared_boundary(prof, np.zeros_like(prof.h))

    def outer_boundary(self, prof):
        raise NotImplementedError("This section does not define an exterior boundary")


class CrossSectionSquared(ChannelSection):
    """Squared passages with independently straight or curved floor and roof.

    sidewall_half_thickness [m] is ONE tube sidewall's normal thickness,
    equivalent to HALF the shared rib between adjacent chamber channels. A
    constant or design-input callable is accepted. A former 1 mm full rib is
    therefore specified as 0.5 mm here.

    Floor thickness comes from SectionProfiles.t_wall. roof_thickness defaults
    to that thickness and may instead be a constant or callable. floor and roof
    independently accept "arc" or "line". Direct width profiles use the
    distance between internal sidewall endpoints. With sector profiles, the
    existing constant-material-width rib construction determines those
    endpoints; choosing a line joins the same endpoints with a chord.
    """

    def __init__(self, sidewall_half_thickness, n_points=8, *, roof_thickness=None,
                 floor="arc", roof="arc", exposure=None, enable_fin=False):
        super().__init__(n_points, exposure=exposure, enable_fin=enable_fin)
        if floor not in ("arc", "line") or roof not in ("arc", "line"):
            raise ValueError("floor and roof must be 'arc' or 'line'")
        self.floor, self.roof = floor, roof
        self._sidewall_value = sidewall_half_thickness
        self._roof_value = roof_thickness
        self._sidewall_evaluator = make_evaluator(
            sidewall_half_thickness, slot="sidewall_half_thickness",
            label="CrossSectionSquared.sidewall_half_thickness")
        self._roof_evaluator = None if roof_thickness is None else make_evaluator(
            roof_thickness, slot="section_thickness", label="CrossSectionSquared.roof_thickness")

    def sidewall_half_thickness(self, prof):
        def compute():
            values = self._evaluate(self._sidewall_evaluator, prof)
            if np.any(~np.isfinite(values)) or np.any(values < 0):
                raise ValueError("sidewall_half_thickness must be finite and non-negative")
            if prof.floor_width is None and np.any(2*values >= self._pitch_base(prof)):
                raise ValueError("sidewall_half_thickness leaves no room for the channel")
            return values
        return prof.cached((self, "sidewall_half_thickness"), compute).copy()

    def roof_thickness(self, prof):
        values = prof.t_wall.copy() if self._roof_evaluator is None else self._evaluate(self._roof_evaluator, prof)
        if np.any(~np.isfinite(values)) or np.any(values < 0):
            raise ValueError("roof_thickness must be finite and nonnegative")
        return values

    def surface_thicknesses(self, prof):
        side = self.sidewall_half_thickness(prof)
        return {"floor": prof.t_wall.copy(), "roof": self.roof_thickness(prof),
                "side_left": side.copy(), "side_right": side.copy()}

    def _pitch_base(self, prof):
        return prof.cross_section_base_point_path[:, 1]*prof.theta

    def _radii(self, prof):
        base = prof.cross_section_base_point_path[:, 1]
        return base-prof.t_wall, base, base+prof.h

    def _transform(self, prof):
        return None if prof.projection is None else np.linalg.inv(prof.projection)

    def full_height(self, prof):
        if self.floor == self.roof == "arc":
            return prof.h.copy()
        floor, roof = self.boundary(prof)
        vector = roof.evaluate([0.5])[:, 0]-floor.evaluate([0.5])[:, 0]
        transform = self._transform(prof)
        if transform is not None:
            vector = np.einsum("nij,nj->ni", transform, vector)
        return vector[:, 0]

    def boundary(self, prof):
        def compute():
            if prof.floor_width is not None:
                from pyskyfire.common.section_geometry import squared_width_boundary
                if prof.projection is not None:
                    raise ValueError("Direct widths are measured in the local section plane; omit projection")
                radius = np.inf if prof.section_radius is None else prof.section_radius
                return squared_width_boundary(prof.floor_width, prof.roof_width, prof.h,
                                              radius, self.floor, self.roof)
            ends = squared_boundary(prof, self.sidewall_half_thickness(prof))
            return tuple(SectionLine(curve.evaluate([0])[:, 0], curve.evaluate([1])[:, 0])
                         if kind == "line" else curve
                         for kind, curve in zip((self.floor, self.roof), ends))
        return prof.cached((self, "boundary"), compute)

    def outer_boundary(self, prof):
        from pyskyfire.common.section_geometry import offset_squared_boundary
        return prof.cached((self, "outer_boundary"), lambda: offset_squared_boundary(
            *self.boundary(prof), prof.t_wall, self.roof_thickness(prof),
            self.sidewall_half_thickness(prof)))

    def thermal_perimeters(self, prof, *, exterior=False):
        # Keep the chamber's established floor-area approximation. Direct
        # dimensions use the actual offset boundary on every surface.
        fractions = self.exposure_fractions(prof)
        if exterior and prof.floor_width is None and self.floor == "arc":
            hot, _, _ = self._radii(prof)
            floor = np.maximum(hot*prof.theta-2*self.sidewall_half_thickness(prof), 0)
            if all(np.all(fractions[name] == 0) for name in ("roof", "side_left", "side_right")):
                result = {name: np.zeros_like(floor) for name in self.surface_names}
            else:
                result = super().thermal_perimeters(prof, exterior=True)
            result["floor"] = floor*fractions["floor"]
            return result
        return super().thermal_perimeters(prof, exterior=exterior)

    def equivalent_angular_width(self, prof):
        _, base, top = self._radii(prof)
        span = top**2-base**2
        result = np.zeros_like(span)
        np.divide(2*self.A_coolant(prof), span, out=result, where=span > 0)
        return result

    def open_fraction(self, prof):
        result = np.zeros_like(prof.theta)
        np.divide(self.equivalent_angular_width(prof), prof.theta, out=result, where=prof.theta > 0)
        return result

    def _fin_dimensions(self, prof):
        # Retain the chamber's first-order axial-height rib model.
        length = prof.h
        return length, 2*self.sidewall_half_thickness(prof), 2*length


class CrossSectionRounded(ChannelSection):
    """Two circular lobes connected by straight tangent sidewalls.

    Both input routes build the same boundary: sector radius/angle/wall
    thickness, or internal floor and roof endpoint widths and height. Equal
    widths give parallel sides; equal widths with full height equal to width
    give a circle. A single wall thickness offsets all four surfaces.

    height_input="full" measures bottom-to-top passage height; "base" measures
    separation of lobe centres. For unequal endpoint widths and full height,
    the less steep of the two possible tangent constructions is selected.
    Invalid combinations raise rather than distorting the lobes. Angular
    input retains its exact prescribed sidewall angle.
    """

    def __init__(self, n_points=16, *, height_input="full", exposure=None,
                 exposure_angles=None, enable_fin=False):
        if height_input not in ("full", "base"):
            raise ValueError("height_input must be 'full' or 'base'")
        super().__init__(n_points, exposure=exposure, exposure_angles=exposure_angles,
                         enable_fin=enable_fin)
        self._height_input = height_input

    @property
    def height_input(self):
        return getattr(self, "_height_input", "full")

    def _sector_inner_radius(self, prof):
        radius = prof.sector_radius
        if radius is None:
            radius = prof.cross_section_base_point_path[:, 1]
        q = radius*np.sin(prof.theta/2)-prof.t_wall
        if np.any(~np.isfinite(q)) or np.any(q <= 0):
            raise ValueError("Rounded inner arc radius must be finite and positive")
        return q

    def minimum_channel_height(self, prof):
        if prof.floor_width is not None:
            from scipy.optimize import brentq
            lower, upper = prof.floor_width, prof.roof_width
            mean, delta = (lower+upper)/2, np.abs(upper-lower)/2
            result = mean.copy()
            for i in range(len(result)):
                if delta[i] == 0:
                    continue
                ratio = delta[i]/mean[i]
                sine = brentq(lambda z: z**3+2*ratio*z*z-ratio, 0, 1/np.sqrt(2), xtol=1e-15)
                result[i] = (mean[i]+delta[i]/sine)/np.sqrt(1-sine*sine)
            return result
        return 2*self._sector_inner_radius(prof)

    def base_height(self, prof):
        if prof.floor_width is not None:
            lower, upper = self.boundary(prof)
            return (upper.centers[:, 0]-lower.centers[:, 0]).copy()
        minimum = self.minimum_channel_height(prof)
        measured = np.asarray(prof.h, float)
        if np.any(~np.isfinite(measured)):
            raise ValueError("Rounded channel height must be finite")
        if self.height_input == "base":
            if np.any(measured < 0):
                raise ValueError("Rounded channel base height must be non-negative")
            return measured.copy()
        tolerance = 16*np.finfo(float).eps*np.maximum(minimum, 1)
        if np.any(measured < minimum-tolerance):
            raise ValueError("Rounded channel height is below the geometric minimum")
        return np.maximum((measured-minimum)/(1+np.sin(prof.theta/2)), 0)

    def full_height(self, prof):
        if prof.floor_width is not None:
            lower, upper = self.boundary(prof)
            return upper.centers[:, 0]-lower.centers[:, 0]+lower.radii+upper.radii
        base = self.base_height(prof)
        return self.minimum_channel_height(prof)+base*(1+np.sin(prof.theta/2))

    def boundary(self, prof):
        def compute():
            if prof.floor_width is not None:
                return rounded_width_boundary(prof.floor_width, prof.roof_width, prof.h,
                                              height_input=self.height_input)
            return rounded_boundary(prof, self.base_height(prof))
        return prof.cached((self, "boundary"), compute)

    def outer_boundary(self, prof):
        return prof.cached((self, "outer_boundary"), lambda: tuple(
            replace(curve, radii=curve.radii+prof.t_wall) for curve in self.boundary(prof)))

    def surface_thicknesses(self, prof):
        return {name: prof.t_wall.copy() for name in self.surface_names}

    def _fin_dimensions(self, prof):
        lower, _ = self.boundary(prof)
        angle = lower.start_angles-np.pi/2
        length = self.base_height(prof)
        return length, 2*prof.t_wall, 2*length*np.cos(angle)
