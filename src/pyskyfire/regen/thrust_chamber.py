from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from pyskyfire.regen.cross_section import SectionProfiles
from pyskyfire.regen.design_inputs import (
    DesignContext,
    DesignInputError,
    contour_providers,
    layout_providers,
    make_evaluator,
    section_providers,
)
from pyskyfire.thermo import transport as _transport


class Wall:
    """Single structural wall or coating layer in a thrust chamber.

    Encapsulates material and thickness information for heat-conduction
    calculations.

    Parameters
    ----------
    material : Material
        Material object providing thermal conductivity and other data.
    thickness : float or callable
        Constant thickness [m] or function `t(x)` returning thickness.
    name : str, optional
        Descriptive label of the layer (e.g. “Copper liner”).
    role : {"main", "liner"}, keyword-only, optional
        Defaults to ``"main"``. Main layers form the channel walls and all
        contribute to channel thickness. Liners lie on the combustion side
        and contribute to the hot-side offset and conduction, but not to
        channel-wall thickness. Circuit layers are ordered hot to cold,
        with liners before main layers.

    Attributes
    ----------
    material : Material
        Thermal material definition.
    _thickness : float or callable
        Underlying storage for the thickness definition.
    name : str or None
        Optional descriptive name.
    role : {"main", "liner"}
        Whether the layer forms the channel wall or lies only on its hot side.

    """

    def __init__(self, material, thickness, name=None, *, role="main"):

        if role not in ("main", "liner"):
            raise ValueError("Wall.role must be 'main' or 'liner'")
        self.role = role
        self.name = name
        self.material = material
        self._thickness = thickness
        self._contour = None
        self._evaluator = make_evaluator(
            thickness,
            slot="wall_thickness",
            label=f"Wall({name!r}).thickness" if name else "Wall.thickness",
        )

    def __setstate__(self, state):
        """Restore legacy walls with a main role and a thickness evaluator."""
        state.setdefault("role", "main")
        state.setdefault("_contour", None)
        if "_evaluator" not in state:
            state["_evaluator"] = make_evaluator(
                state["_thickness"],
                slot="wall_thickness",
                label="Wall.thickness",
            )
        self.__dict__.update(state)

    def bind_contour(self, contour):
        """Attach the contour that tier A design inputs are evaluated against.

        Called by :class:`CoolingCircuit`. A wall whose thickness only depends
        on ``x`` never needs this.
        """
        self._contour = contour

    def thickness(self, x):
        """Return the layer thickness [m] at ``x``, which may be an array."""
        constant = self._evaluator.constant_value
        values = np.asarray(x, dtype=float)
        if constant is not None:
            if values.ndim:
                return np.full(values.shape, constant, dtype=float)
            return constant

        result = np.array(
            [
                self._evaluator(
                    DesignContext(contour_providers(self._contour, float(value)))
                )
                for value in values.ravel()
            ],
            dtype=float,
        ).reshape(values.shape)
        return float(result) if result.ndim == 0 else result


@dataclass
class FilmCooling:
    """
    User-facing film cooling inputs attached to a ThrustChamber.

    Parameters
    ----------
    x_fraction : float
        Signed axial fraction in [-1, 1], where:
            -1 -> chamber start
             0 -> throat
            +1 -> nozzle exit
        This is resolved to ``x`` by ``ThrustChamber``.
    coolant_mass_flow_rate : float
        Total coolant mass flow injected into the film [kg/s].
    film_injection_perimeter : float
        Wetted injection perimeter [m].
    liquid_absorptivity : float
        Liquid-film absorptivity for radiation [0..1].
    mole_fraction_H2O : float
        User-supplied mole fraction of H2O used by the current Grisson
        gas-radiation submodel.
    mole_fraction_CO2 : float
        User-supplied mole fraction of CO2 used by the current Grisson
        gas-radiation submodel.
    turbulence_intensity : float, optional
        Free-stream turbulence intensity used by Grisson's correlations.
    x : float | None, optional
        Absolute axial injection location [m]. This should be populated by
        ``ThrustChamber`` after it resolves ``x_fraction``.
    """

    coolant_transport: _transport.FluidTransport
    x_fraction: float
    coolant_mass_flow_rate: float
    film_injection_perimeter: float
    liquid_absorptivity: float
    mole_fraction_H2O: float
    mole_fraction_CO2: float
    turbulence_intensity: float = 0.0
    x: Optional[float] = None


class CoolingCircuit:
    """
    Simulation-only representation of a cooling circuit.

    ``walls`` are ordered from hot gas to coolant. Any ``role="liner"``
    layers precede the main layers; at least one main layer is required.
    Wall-like objects without a role are treated as main layers.

    ``enable_fin=True`` enables this circuit's rib/fin heat-transfer model.
    ``rib_material`` optionally supplies its material (with ``get_k(T)``).
    When omitted, the coolant-facing main-wall material is used. Disabling
    fins removes their thermal enhancement, not their channel geometry.
    """

    #: Derived section geometry; see :meth:`_prof` and
    #: :meth:`section_profiles_at`. Defined on the class so that a partially
    #: populated circuit, built through ``__new__`` for a plain geometry query,
    #: still has somewhere to look.
    _section_cache = None

    def __init__(
        self,
        name,
        contour,
        cross_section,  # ChannelSection
        span: list[float],
        placement,  # ChannelPlacement-like (sim uses only counts/lanes if needed)
        channel_height,  # callable x -> h
        walls,  # list of wall objects with .thickness(x)
        coolant_transport,
        roughness,
        hot_gas_surface_area_multiplier=1.0,
        *,
        enable_fin: bool = True,
        rib_material=None,
    ):
        self.name = name
        self.contour = contour
        self.cross_section = cross_section
        self.placement = placement
        self.coolant_transport = coolant_transport
        self.walls = walls
        has_main = False
        for wall in self.walls:
            role = getattr(wall, "role", "main")
            if role not in ("main", "liner"):
                raise ValueError("Wall.role must be 'main' or 'liner'")
            if role == "liner" and has_main:
                raise ValueError("Liner walls must precede main walls in hot-to-cold order")
            has_main = has_main or role == "main"
            # Duck-typed wall stubs only have to provide ``thickness(x)``; they
            # simply forgo contour-derived design inputs.
            if hasattr(wall, "bind_contour"):
                wall.bind_contour(contour)
        if not has_main:
            raise ValueError("A cooling circuit requires at least one main wall")

        self._chamber = None

        self._channel_height = channel_height
        self._channel_height_evaluator = make_evaluator(
            channel_height,
            slot="channel_height",
            label=f"{name}: channel_height",
        )
        self._roughness = roughness
        self._roughness_evaluator = make_evaluator(
            roughness,
            slot="roughness",
            label=f"{name}: roughness",
        )
        self._hot_gas_surface_area_multiplier = hot_gas_surface_area_multiplier
        self._hot_gas_surface_area_multiplier_evaluator = make_evaluator(
            hot_gas_surface_area_multiplier,
            slot="hot_gas_surface_area_multiplier",
            label=f"{name}: hot_gas_surface_area_multiplier",
        )
        constant = self._hot_gas_surface_area_multiplier_evaluator.constant_value
        if constant is not None and (not np.isfinite(constant) or constant <= 0.0):
            raise ValueError(
                "hot_gas_surface_area_multiplier must be finite and positive"
            )


        self.enable_fin = bool(enable_fin)
        self.rib_material = rib_material

        if span[0] > span[1]:
            self.span = (span[1], span[0])
            self.direction = -1
        else:
            self.span = span
            self.direction = 1

        # to be set later by set_centerline / set_channel_* / set_x_domain
        self.reference_path = None
        self.cross_section_base_point_path = None
        self.cooling_channel_centroid_path = None
        self.x_domain = None
        self.channel_heights = None
        self.channel_angular_width = None
        self.helix_angle_vals = None
        self.is_helical = False
        self.reference_angle_origin_rad = None

    def __getstate__(self):
        """Leave the derived section cache out; it rebuilds on demand."""
        state = self.__dict__.copy()
        state["_section_cache"] = None
        return state

    def __setstate__(self, state):
        """Load circuits serialized before the angular-width rename."""
        state.setdefault("_section_cache", None)
        if "channel_angular_width" not in state and "channel_width" in state:
            state["channel_angular_width"] = state.pop("channel_width")
        if "centerlines" in state or "coolant_centerlines" in state:
            raise ValueError("Saved channel geometry uses obsolete path definitions; regenerate the simulation")
        state.setdefault("reference_angle_origin_rad", None)
        state.setdefault("enable_fin", True)
        state.setdefault("rib_material", None)
        self.__dict__.update(state)

    # ---------------- basic getters ----------------
    @property
    def n_channels(self) -> int:
        """Total number of physical channels carrying this circuit's flow.

        Leaves times channels per leaf. ``placement.channel_count()`` counts
        leaves only, because it feeds the circumferential packing; a radial
        stack shares a leaf but not a passage, so anything that divides mass
        flow or accumulates volume needs this instead.
        """
        placement = self.placement
        return int(placement.n_channel_positions) * int(
            getattr(placement, "n_channels_per_leaf", 1)
        )

    # ---------------- design-input contexts ----------------
    def bind_chamber(self, chamber):
        """Attach the owning chamber, which supplies cross-circuit counts."""
        self._chamber = chamber
        # The chamber is what makes 'n_channels_all_circuits' resolvable.
        self._invalidate_section_cache()

    def _channels_all_circuits_provider(self, x):
        def provider():
            if self._chamber is None:
                raise DesignInputError(
                    f"{self.name}: 'n_channels_all_circuits' counts channels "
                    "across every circuit, so it is only available once the "
                    "circuit belongs to a ThrustChamber"
                )
            return float(self._chamber.number_of_channels(x, occluding_only=True))

        return provider

    def _section_provider_values(self, x):
        """Interpolate whichever section quantities have been built at ``x``.

        Anything not built yet comes back as ``None`` and is then omitted from
        the context, which is what makes the later part of tier C unavailable
        to slots that are evaluated earlier in the geometry build.
        """
        if self.x_domain is None:
            return None
        heights = self.channel_heights
        local_sector = getattr(self, "channel_local_sector", None)
        gamma = self.helix_angle_vals
        return {
            "h": (
                None
                if heights is None
                else float(np.interp(x, self.x_domain, heights))
            ),
            "theta": (
                None
                if local_sector is None
                else float(np.interp(x, self.x_domain, local_sector))
            ),
            "gamma": (
                None if gamma is None else float(np.interp(x, self.x_domain, gamma))
            ),
        }

    def design_context(self, x, *, include_section: bool) -> DesignContext:
        """Build the design-input context for this circuit at ``x``."""
        x = float(x)
        providers = contour_providers(self.contour, x)
        providers.update(
            layout_providers(
                n_channel_positions=int(self.placement.n_channel_positions),
                n_channels=int(self.n_channels),
                n_channels_all_circuits=self._channels_all_circuits_provider(x),
                t_wall=lambda: float(self.channel_wall_thickness(x)),
                t_wall_total=lambda: float(self.total_thickness(x)),
                t_hot_layer=lambda: float(self.hot_side_layer_thickness(x)),
            )
        )
        if include_section:
            section = self._section_provider_values(x)
            if section is not None:
                providers.update(
                    section_providers(
                        r=float(self.contour.r(x)),
                        r_base=(
                            float(self.contour.r(x))
                            if self.cross_section_base_point_path is None
                            else float(np.interp(
                                x, self.x_domain, self.cross_section_base_point_path[:, 1]))
                        ),
                        h=section["h"],
                        gamma=section["gamma"],
                        theta=section["theta"],
                        t_wall=float(self.channel_wall_thickness(x)),
                    )
                )
        return DesignContext(providers)

    def _evaluate(self, evaluator, x, *, include_section: bool):
        """Evaluate a design input over a scalar or array of ``x``."""
        constant = evaluator.constant_value
        values = np.asarray(x, dtype=float)
        if constant is not None:
            if values.ndim:
                return np.full(values.shape, constant, dtype=float)
            return constant
        result = np.array(
            [
                evaluator(
                    self.design_context(value, include_section=include_section)
                )
                for value in values.ravel()
            ],
            dtype=float,
        ).reshape(values.shape)
        return float(result) if result.ndim == 0 else result

    # ---------------- design inputs ----------------
    def channel_height(self, x):
        """Return the supplied height [m], using the section's input convention.

        Use ``channel_full_heights`` or ``channel_base_heights`` to report
        resolved physical heights on ``x_domain`` after chamber assembly.
        """
        return self._evaluate(
            self._channel_height_evaluator, x, include_section=False
        )

    def roughness(self, x):
        """Return the absolute channel roughness [m] at ``x``."""
        return self._evaluate(self._roughness_evaluator, x, include_section=True)

    def hot_gas_surface_area_multiplier(self, x):
        """Return the circuit-specific effective hot-side area multiplier."""
        values = np.asarray(
            self._evaluate(
                self._hot_gas_surface_area_multiplier_evaluator,
                x,
                include_section=True,
            ),
            dtype=float,
        )
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError(
                f"{self.name}: hot_gas_surface_area_multiplier must return "
                "finite, positive values"
            )
        return float(values) if values.ndim == 0 else values

    def _wall_thickness(self, x, role=None):
        thickness = np.zeros_like(np.asarray(x, dtype=float))
        for wall in self.walls:
            if role is None or getattr(wall, "role", "main") == role:
                thickness = thickness + np.asarray(wall.thickness(x), dtype=float)
        return float(thickness) if thickness.ndim == 0 else thickness

    def total_thickness(self, x):
        """Return the full hot-to-coolant stack thickness, including liners."""
        return self._wall_thickness(x)

    def channel_wall_thickness(self, x):
        """Return the combined main-wall thickness used for channel geometry."""
        return self._wall_thickness(x, role="main")

    def hot_side_layer_thickness(self, x):
        """Return the combined thickness of explicitly marked hot-side liners."""
        return self._wall_thickness(x, role="liner")

    # ---------------- profiles bundle ----------------
    def _prof(self, centerline: np.ndarray) -> SectionProfiles:
        """
        Build a SectionProfiles with trivial frames (analytics don't use them).

        The bundle is pure geometry, so it is cached against the path it was
        built from and reused until :meth:`_invalidate_section_cache` drops it.
        Assembling one walks every station of every wall, evaluating any
        design-input function a station at a time, and the solver reaches a
        single station through :meth:`section_profiles_at` on every residual
        it forms, so rebuilding it there dominated the whole simulation.

        Staleness is bounded the same way the ``*_vals`` arrays that
        :meth:`precompute_thermal_properties` writes are: change the geometry
        and you have to rebuild, which for the profile setters happens
        automatically and for a wholesale path swap falls out of the path no
        longer being the one the bundle was built from.
        """
        return self._section_cache_for(centerline)[1]

    def _section_cache_for(self, centerline: np.ndarray):
        """Return ``(path, bundle, stations)`` for ``centerline``, building it once.

        The single-station profiles live in the same entry as the bundle they
        were interpolated from, so swapping the path cannot leave stations
        behind that describe the old geometry.
        """
        cache = self._section_cache
        if cache is None or cache[0] is not centerline:
            cache = (centerline, self._build_prof(centerline), {})
            self._section_cache = cache
        return cache

    def _invalidate_section_cache(self):
        """Drop the cached section geometry after a change that redefines it."""
        self._section_cache = None

    def _build_prof(self, centerline: np.ndarray) -> SectionProfiles:
        x = centerline[:, 0]

        # All main layers contribute to arc/ligament thickness. Hot-side
        # liners offset the reference but do not form the tube walls.
        t_channel_wall = self.channel_wall_thickness(x)

        return SectionProfiles(
            h=np.asarray(self.channel_heights, float),
            theta=np.asarray(self.channel_local_sector, float),
            # theta=np.asarray(self.channel_angular_width, float),
            t_wall=np.asarray(t_channel_wall, float),
            cross_section_base_point_path=np.asarray(centerline, float),
            sector_radius=getattr(self, "section_sector_radius", None),
            projection=getattr(self, "section_projection", None),
            azimuthal_pitch=getattr(self, "channel_angular_width", None),
            winding_rate=getattr(self, "section_winding_rate", None),
            **self._design_input_context(x),
        )

    def _design_input_context(self, x):
        """Return the optional ``SectionProfiles`` fields, best effort.

        Everything here only feeds design-input functions. A quantity that
        cannot be determined is left out rather than raising, so that a
        partially built circuit can still be asked for plain geometry; a
        design input that needs the missing quantity is what raises.
        """
        helix = getattr(self, "helix_angle_vals", None)
        placement = getattr(self, "placement", None)
        context = {
            "contour": getattr(self, "contour", None),
            "gamma": None if helix is None else np.asarray(helix, float),
            "n_channels_all": self._occluding_channel_counts(x),
            "t_wall_total": np.asarray(self.total_thickness(x), float),
            "t_hot_layer": np.asarray(self.hot_side_layer_thickness(x), float),
        }
        if placement is not None:
            context["n_channel_positions"] = int(placement.n_channel_positions)
            context["n_channels"] = int(self.n_channels)
        return context

    def _occluding_channel_counts(self, x):
        """Return the occluding channel count at each station, if knowable."""
        if getattr(self, "_chamber", None) is None:
            return None
        return np.array(
            [
                self._chamber.number_of_channels(float(value), occluding_only=True)
                for value in np.asarray(x, dtype=float)
            ],
            dtype=float,
        )

    # ---------------- preprocessing ----------------
    def precompute_thermal_properties(self):
        """
        Precompute passage geometry along a representative coolant path.

        ``cross_section_base_point_path`` anchors section geometry;
        ``cooling_channel_centroid_path`` controls passage length and curvature.
        """
        if self.cross_section_base_point_path is None or self.cooling_channel_centroid_path is None:
            raise RuntimeError("Channel paths have not been constructed")
        section_centerline = self.cross_section_base_point_path
        flow_centerline = self.cooling_channel_centroid_path
        x_vals = self.x_domain
        # Normal-plane centroids can be axially displaced from their nominal
        # station. Differentiate against the solver stations, not their own X.
        derivatives = np.gradient(flow_centerline, x_vals, axis=0, edge_order=min(2, len(x_vals)-1))
        self.ds_dx_vals = np.sqrt(derivatives[:, 0]**2 + derivatives[:, 1]**2
                                 + (flow_centerline[:, 1]*derivatives[:, 2])**2)
        # Sections are measured in the plane normal to the nominal winding, so
        # that winding -- not the centroid it encloses -- is the path a section
        # area has to be swept along to recover a volume.
        self.ds_dx_reference_vals = cylindrical_curve_arc_length_scale(self.reference_path)

        prof = self._prof(section_centerline)

        hot_perimeter = self.cross_section.effective_hot_gas_perimeter(
            prof,
            self.hot_gas_surface_area_multiplier(x_vals),
        )
        cold_perimeter = self.cross_section.P_coolant(prof)

        self.dA_dx_thermal_exhaust_vals = hot_perimeter * self.ds_dx_vals
        self.dA_dx_thermal_coolant_vals = cold_perimeter * self.ds_dx_vals

        A_vals = self.cross_section.A_coolant(prof)
        self.A_coolant_vals = A_vals
        self.dA_dx_coolant_vals = np.gradient(A_vals, x_vals)

        self.Dh_coolant_vals = self.cross_section.Dh_coolant(prof)

        self.signed_curvature_vals = signed_curvature(flow_centerline)
        self.radius_of_curvature_vals = _radius_from_curvature(
            self.signed_curvature_vals
        )

    def section_profiles_at(self, x: float) -> SectionProfiles:
        """Interpolate the resolved physical section on the solver grid.

        Cached per station, because the solver holds ``x`` fixed while it
        iterates one wall node's temperatures and so asks for the same station
        repeatedly. Sharing the profile also shares the geometry each section
        memoises on it; see :meth:`SectionProfiles.cached`.
        """
        station = float(x)
        stations = self._section_cache_for(self.cross_section_base_point_path)[2]
        cached = stations.get(station)
        if cached is None:
            # Solver grids are small and fixed, so this only grows to a few
            # hundred entries; the bound just stops a sweep over many stations
            # from retaining all of them.
            if len(stations) >= 4096:
                stations.clear()
            cached = self._build_section_profiles_at(station)
            stations[station] = cached
        return cached

    def _build_section_profiles_at(self, x: float) -> SectionProfiles:
        from dataclasses import replace
        prof = self._prof(self.cross_section_base_point_path)
        values = {}
        for name in ("h", "theta", "t_wall", "cross_section_base_point_path",
                     "sector_radius", "projection", "azimuthal_pitch",
                     "winding_rate", "gamma", "t_wall_total", "t_hot_layer",
                     "n_channels_all"):
            value = getattr(prof, name)
            if value is None:
                continue
            array = np.asarray(value)
            flat = array.reshape(len(self.x_domain), -1)
            values[name] = np.array([np.interp(x, self.x_domain, flat[:, i])
                                     for i in range(flat.shape[1])]).reshape((1,)+array.shape[1:])
        return replace(prof, **values)

    def R_coolant_per_len(self, x: float, h_c: float, T_wall_rep: float) -> float:
        """
        Effective coolant-side thermal resistance per unit *axial length* [K m / W].

        Cross-section returns resistance per unit channel length (s). Convert to per-x using:
            R_x = R_s / (ds/dx)
        """

        Aprime_cool = float(self.dA_dx_thermal_coolant(x))

        # If fin/rib enhancement is disabled, use baseline (old) behavior.
        if not getattr(self, "enable_fin", False):
            var = 1.0 / (h_c * Aprime_cool)
            # print(f"disabled: {var}")
            return var

        # Get 1-point profiles at x (needed by cross-section closure)
        prof_x = self.section_profiles_at(x)

        # A material override applies only to the fin model. Through-wall
        # conduction continues to use each wall layer's own material.
        material = self.rib_material
        if material is None:
            material = self.walls[-1].material
        k_rib = float(material.get_k(T_wall_rep))

        # Cross-section provides R per unit channel length
        R_s = float(
            self.cross_section.R_coolant_per_len(prof_x, np.array([h_c]), k_rib)[0]
        )

        # Convert to per unit axial length using ds/dx from your already-defined geometry
        ds_dx = float(self.ds_dx(x))
        ds_dx = max(ds_dx, 1e-30)

        var = R_s / ds_dx
        # print(f"enabled: {var}")
        return var

    def compute_volume(self):
        """Set ``self.volume``: the wetted coolant volume of the circuit [m³].

        One channel encloses :math:`\\int A\\,ds` along its own centreline, so
        the circuit holds :math:`n_{chan}\\int A\\,(ds/dx)\\,dx`.

        ``A_coolant`` is the flow area in the plane normal to the *reference
        path*, not an axial slice, so it is that path's ``ds/dx`` the area is
        swept along here. :meth:`ThrustChamber.build_channel_angular_widths`
        narrows the sector by :math:`\\cos\\gamma` for a helix angle
        :math:`\\gamma`, while ``ds/dx`` lengthens the path by the reciprocal
        on a cylinder. A helical circuit therefore encloses the same volume as
        an axial one of the same azimuthal pitch, which is the physically
        correct answer -- only the path is longer and the section
        correspondingly narrower. ``ds_dx_vals`` follows the coolant centroid
        instead, which sits further out and so runs slightly further still;
        pairing it with the section area would make the swept volume depend on
        the helix angle.

        Trapezoidal on the circuit's own ``x_domain``, so the result converges
        with axial resolution rather than being exact at any given ``n_nodes``.
        """
        x_vals = np.asarray(self.x_domain, dtype=float)
        volume_per_channel = np.trapezoid(
            self.A_coolant_vals * self.ds_dx_reference_vals,
            x_vals,
        )
        # Volume is signed by the integration direction; the grid is built
        # ascending today, but do not depend on that.
        self.volume = abs(float(volume_per_channel)) * self.n_channels

    # ---------------- setters ----------------
    def section_planes(self, reference_path=None, *, axis_normal=False):
        """Return physical section planes, optionally at a resolved layout seed."""
        from .channel_construction import build_channel_construction
        from .channel_helpers import build_axis_normal_section_planes
        construction = build_channel_construction(self, reference_path)
        if axis_normal:
            return build_axis_normal_section_planes(construction.profiles.cross_section_base_point_path)
        return construction.flow_planes

    def set_channel_angular_width(self, angular_widths_rad):
        self.channel_angular_width = np.asarray(angular_widths_rad, float)
        self._invalidate_section_cache()

    def set_channel_local_sector(self, local_sectors):
        self.channel_local_sector = np.asarray(local_sectors, float)
        self._invalidate_section_cache()

    def set_channel_height(self, heights):
        self.channel_heights = np.asarray(heights, float)
        self._invalidate_section_cache()

    @property
    def channel_full_heights(self):
        """Full internal coolant heights on ``x_domain`` [m], for reporting."""
        converter = getattr(self.cross_section, "full_height", None)
        if converter is None:
            return np.asarray(self.channel_heights, float).copy()
        return converter(self._prof(self.cross_section_base_point_path))

    @property
    def channel_base_heights(self):
        """Rounded base heights on ``x_domain`` [m]; full heights for squared."""
        converter = getattr(self.cross_section, "base_height", None)
        if converter is None:
            return np.asarray(self.channel_heights, float).copy()
        return converter(self._prof(self.cross_section_base_point_path))

    def set_x_domain(self, x_domain):
        self.x_domain = np.asarray(x_domain, float)
        self._invalidate_section_cache()

    def finalize(self):
        self.precompute_thermal_properties()
        self.compute_volume()

    # ---------------- interpolants ----------------
    def dA_dx_thermal_exhaust(self, x):
        return np.interp(x, self.x_domain, self.dA_dx_thermal_exhaust_vals)

    def dA_dx_thermal_coolant(self, x):
        return np.interp(x, self.x_domain, self.dA_dx_thermal_coolant_vals)

    def A_coolant(self, x):
        return np.interp(x, self.x_domain, self.A_coolant_vals)

    def dA_dx_coolant(self, x):
        return np.interp(x, self.x_domain, self.dA_dx_coolant_vals)

    def Dh_coolant(self, x):
        return np.interp(x, self.x_domain, self.Dh_coolant_vals)

    def radius_of_curvature(self, x):
        # Interpolate curvature, not radius: interpolation involving infinity
        # is ill-conditioned at a straight-to-curved transition.
        curvature = np.interp(x, self.x_domain, self.signed_curvature_vals)
        return float(_radius_from_curvature(curvature))

    def wedge_angle(self, x):
        return np.interp(x, self.x_domain, self.channel_angular_width)

    def local_sector_angle(self, x):
        return np.interp(x, self.x_domain, self.channel_local_sector)

    def ds_dx(self, x):
        return np.interp(x, self.x_domain, self.ds_dx_vals)  # added for helical support


class ThrustChamber:
    """
    Simulation-only thrust-chamber assembly.

    Notes
    -----
    • Define x-domains for circuits from their fractional spans
    • Build nominal reference, physical base-point and coolant centroid paths
    • Compute per-circuit wedge angles (theta) and heights (h)
    • Trigger each circuit's precompute/finalize (A, Dh, perimeters, volume)
    • (Optionally) trigger combustion aerothermodynamics
    """

    def __init__(
        self,
        contour,
        cooling_circuits,
        combustion_transport,
        optimal_values=None,
        K_factor: float = 0.3,
        n_nodes: int = 50,
        h_gas_corr: float = 1.0,
        h_cold_corr: float = 1.0,
        compute_gas: bool = True,
        film_cooling: FilmCooling | None = None,
    ):
        self.contour = contour
        self.cooling_circuits = cooling_circuits
        self.combustion_transport = combustion_transport
        self.optimal_values = optimal_values

        self.n_nodes = int(n_nodes)
        self.K_factor = float(K_factor)
        self.h_gas_corr = float(h_gas_corr)
        self.h_cold_corr = float(h_cold_corr)

        # --- geometry inputs for circuits (sim only) ---
        # Circuits need a back-reference before any design input is evaluated:
        # 'n_channels_all_circuits' is a cross-circuit quantity.
        for circuit in self.cooling_circuits:
            circuit.bind_chamber(self)

        self.build_circuit_x_domain()
        self.build_channel_heights()
        self.build_reference_paths()
        self.build_channel_angular_widths()
        self.build_physical_channel_paths()

        # --- per-circuit preprocessing for solver ---
        for circuit in self.cooling_circuits:
            circuit.finalize()

        # --- hot-gas side (optional/on-demand) ---
        if compute_gas and hasattr(
            self.combustion_transport, "compute_aerothermodynamics"
        ):
            self.combustion_transport.compute_aerothermodynamics(self.contour)

        self.film_cooling = film_cooling
        if self.film_cooling is not None:
            self.film_cooling.x = self._resolve_signed_fraction_to_x(
                self.film_cooling.x_fraction
            )

    # ------------------------------------------------------------------
    # Domain → nominal winding → physical sections → coolant centroid
    # ------------------------------------------------------------------
    def build_circuit_x_domain(self):
        """
        Convert each circuit's fractional span into an x-grid whose size
        scales with the fraction of the overall contour length covered.
        Uses at least 3 nodes per circuit.
        """
        x_min = float(self.contour.xs[0])
        x_max = float(self.contour.xs[-1])
        L_tot = max(1e-30, x_max - x_min)  # guard zero-length contours

        MIN_NODES = 3

        for circuit in self.cooling_circuits:
            f0, f1 = circuit.span

            # Map signed fractional spans to absolute x, keeping your original convention
            x0 = f0 * (x_max if f0 >= 0.0 else -x_min)
            x1 = f1 * (x_max if f1 >= 0.0 else -x_min)

            L_seg = abs(x1 - x0)
            frac = L_seg / L_tot
            n_here = max(MIN_NODES, int(round(self.n_nodes * frac)))

            if x0 <= x1:
                x_domain = np.linspace(x0, x1, n_here)
            else:
                # reversed span: keep coolant marching direction by reversing
                x_domain = np.linspace(x1, x0, n_here)[::-1]

            circuit.set_x_domain(x_domain)

    def build_reference_paths(self):
        """Build the ideal-contour winding on a shared integration grid."""
        if not self.cooling_circuits:
            return
        shared_x = np.unique(np.concatenate([c.x_domain for c in self.cooling_circuits]))
        for circuit in self.cooling_circuits:
            x = circuit.x_domain
            rates = np.array([circuit.placement.dtheta_dx(float(v), self.contour, circuit)
                              for v in shared_x])
            angle = np.zeros_like(shared_x)
            angle[1:] = np.cumsum((rates[1:]+rates[:-1])*np.diff(shared_x)/2)
            sampled = np.interp(x, shared_x, angle)
            circuit.reference_angle_origin_rad = float(sampled[0])
            sampled -= sampled[0]
            radius = np.array([float(self.contour.r(float(v))) for v in x])
            circuit.reference_path = np.column_stack((x, radius, sampled))
            circuit.helix_angle_vals = np.array([
                circuit.placement.helix_angle(float(v), self.contour, circuit) for v in x])
            circuit.is_helical = bool(np.any(np.abs(circuit.helix_angle_vals) > 1e-12))
            # The helix angle reaches design inputs as 'gamma'.
            circuit._invalidate_section_cache()

    def build_physical_channel_paths(self):
        """Resolve section anchors, lobe sizes and centroid paths together."""
        from .channel_construction import build_channel_construction
        for circuit in self.cooling_circuits:
            construction = build_channel_construction(circuit)
            prof = construction.profiles
            circuit.cross_section_base_point_path = prof.cross_section_base_point_path
            circuit.section_sector_radius = prof.sector_radius
            circuit.section_projection = prof.projection
            circuit.section_winding_rate = prof.winding_rate
            circuit.channel_local_sector = prof.theta
            circuit.cooling_channel_centroid_path = construction.centroid_path(circuit.cross_section)
            circuit.gas_sector_area_vals = construction.gas_area
            circuit._invalidate_section_cache()

    def build_channel_angular_widths(self):
        """
        Compute wedge angle (theta) arrays used by cross-section analytics.

        Priority:
        1) If placement provides channel_angular_width(x), use it.
        2) If placement is 'internal' (non-occluding), use height as an angular-width proxy.
        3) Otherwise, default to an even packing among surface-occluding channels:
           theta(x) = 2π / n_occ(x), where n_occ is provided by the circuit group.
        """
        for circuit in self.cooling_circuits:
            xs = circuit.x_domain
            p = circuit.placement

            if p.has_channel_angular_width:
                angular_widths = np.array(
                    [p.channel_angular_width(x, self.contour, circuit) for x in xs],
                    dtype=float,
                )
                local_sector = angular_widths.copy()

            elif getattr(p, "occludes", True) is False:
                # internal / non-occluding: use height as an angular proxy (sim only)
                angular_widths = np.array(
                    [circuit.channel_height(x) for x in xs],
                    dtype=float,
                )
                local_sector = angular_widths.copy()

            else:
                # uniform angular share among all occluding channels present at x
                n_occ = np.array(
                    [self.number_of_channels(x, occluding_only=True) for x in xs],
                    dtype=float,
                )
                gamma = np.array(
                    [p.helix_angle(x, self.contour, circuit) for x in xs],
                    dtype=float,
                )
                # print
                angular_widths = np.where(
                    n_occ > 0.0,
                    2.0 * np.pi / n_occ,
                    0.0,
                )
                local_sector = np.where(
                    n_occ > 0.0, 2.0 * np.pi * np.cos(gamma) / n_occ, 0.0
                )  # multiplying with cos gamma for helical implementation

            circuit.set_channel_angular_width(angular_widths)
            circuit.set_channel_local_sector(local_sector)

    def build_channel_heights(self):
        """Evaluate the per-circuit height law h(x) on each circuit domain."""
        for circuit in self.cooling_circuits:
            xs = circuit.x_domain
            h = np.array([circuit.channel_height(x) for x in xs], dtype=float)
            circuit.set_channel_height(h)

    @property
    def coolant_volume(self) -> float:
        """Wetted coolant volume of every cooling circuit combined [m³].

        See :meth:`CoolingCircuit.compute_volume` for what one circuit's
        ``volume`` covers. Circuits that overlap axially, such as an RL10-style
        half pass and full pass, each contribute their own channels.
        """
        return float(sum(circuit.volume for circuit in self.cooling_circuits))

    def number_of_channels(self, x, *, occluding_only=False):
        """
        Return the total number of active channels at position `x`.

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        occluding_only : bool, optional
            If True, count only circuits that occlude the wall surface.

        Returns
        -------
        int
            Total number of channels currently active at `x`.
        """
        total_channels = 0
        for circuit in self.cooling_circuits:
            # Assume circuit.x_domain is a numpy array of x values for this circuit.
            x_start = min(circuit.x_domain[0], circuit.x_domain[-1])
            x_end = max(circuit.x_domain[0], circuit.x_domain[-1])
            if x_start <= x <= x_end:
                if occluding_only and not circuit.placement.occludes:
                    continue
                # total_channels += circuit.placement.n_channel_positions
                total_channels += circuit.placement.channel_count()
        return total_channels

    def _resolve_signed_fraction_to_x(self, f: float) -> float:
        if not (-1.0 <= f <= 1.0):
            raise ValueError("film_cooling.x_fraction must be between -1 and 1")

        x_start = float(self.contour.xs[0])
        x_throat = float(self.contour.x_t)
        x_end = float(self.contour.xs[-1])

        if f < 0.0:
            return x_throat + f * (x_throat - x_start)
        return x_throat + f * (x_end - x_throat)


def _validate_cylindrical_curve(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (N, 3): [x, r, theta]")
    if points.shape[0] < 3:
        raise ValueError("at least three centerline points are required")
    if np.any(np.diff(points[:, 0]) == 0.0):
        raise ValueError("centerline x coordinates must be distinct")
    return points


def cylindrical_curve_arc_length_scale(points: np.ndarray) -> np.ndarray:
    """Return ``ds/dx`` for a cylindrical curve ``[x, r, theta]``."""
    points = _validate_cylindrical_curve(points)
    x, radius, theta = points.T
    dr_dx = np.gradient(radius, x, edge_order=2)
    dtheta_dx = np.gradient(theta, x, edge_order=2)
    return np.sqrt(1.0 + dr_dx**2 + (radius * dtheta_dx) ** 2)


def signed_curvature(points: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """Return signed 3-D curvature for ``[x, r, theta]`` centerline points.

    Cylindrical coordinates are first converted to Cartesian coordinates, so
    azimuthal/helical bending contributes to curvature.  Positive curvature
    bends radially away from the chamber axis (the RL10 concave convention),
    while negative curvature bends toward it (convex).  Locally straight
    portions return exactly zero.

    Curvature is estimated from the circumcircle through each point and its
    two neighbours (Menger curvature).  This avoids taking two numerical
    derivatives of the sampled centerline, which strongly amplifies small
    interpolation errors and makes nearly straight portions oscillate between
    very large positive and negative radii.
    """
    points = _validate_cylindrical_curve(points)
    x, radius, theta = points.T
    xyz = np.column_stack([x, radius * np.cos(theta), radius * np.sin(theta)])

    # Use centered triples in the interior and the nearest one-sided triple at
    # each endpoint.  Thus curvature at a straight-to-curved transition is
    # confined to the triple that actually spans the transition.
    middle = np.clip(np.arange(len(points)), 1, len(points) - 2)
    before = xyz[middle - 1]
    center = xyz[middle]
    after = xyz[middle + 1]

    chord_before = center - before
    chord_after = after - center
    chord_across = after - before
    length_before = np.linalg.norm(chord_before, axis=1)
    length_after = np.linalg.norm(chord_after, axis=1)
    length_across = np.linalg.norm(chord_across, axis=1)
    denominator = length_before * length_after * length_across
    if np.any(denominator <= np.finfo(float).tiny):
        raise ValueError("centerline contains coincident points")

    # 4 * triangle area / (a*b*c), with twice the triangle area given
    # by the cross-product magnitude.
    magnitude = (
        2.0 * np.linalg.norm(np.cross(chord_before, chord_after), axis=1) / denominator
    )

    # The change in unit chord direction points toward the local centre of
    # curvature.  Only its radial projection is needed for the RL10 sign.
    tangent_change = (
        chord_after / length_after[:, None] - chord_before / length_before[:, None]
    )
    turning = np.linalg.norm(tangent_change, axis=1)
    principal_normal = np.divide(
        tangent_change,
        turning[:, None],
        out=np.zeros_like(tangent_change),
        where=turning[:, None] > 0.0,
    )
    radial_unit = np.column_stack(
        [
            np.zeros_like(middle, dtype=float),
            np.cos(theta[middle]),
            np.sin(theta[middle]),
        ]
    )
    radial_component = magnitude * np.einsum("ij,ij->i", principal_normal, radial_unit)
    direction = np.sign(radial_component)

    # A locally axial curvature vector can have a vanishing radial projection.
    # Retain the established meridional sign convention as a fallback.
    slope_before = (radius[middle] - radius[middle - 1]) / (x[middle] - x[middle - 1])
    slope_after = (radius[middle + 1] - radius[middle]) / (x[middle + 1] - x[middle])
    d2r_dx2 = 2.0 * (slope_after - slope_before) / (x[middle + 1] - x[middle - 1])
    ambiguous = np.abs(radial_component) <= eps
    direction[ambiguous] = np.sign(d2r_dx2[ambiguous])
    direction[direction == 0.0] = 1.0

    curvature = direction * magnitude
    curvature[magnitude < eps] = 0.0
    return curvature


def _radius_from_curvature(curvature, eps: float = 1e-9):
    curvature = np.asarray(curvature, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        radius = np.where(np.abs(curvature) < eps, np.inf, 1.0 / curvature)
    return float(radius) if radius.ndim == 0 else radius


def radius_of_curvature(
    points: np.ndarray,
    axis: str = "x",
    eps: float = 1e-9,
) -> np.ndarray:
    """Return signed 3-D radius of curvature for ``[x, r, theta]`` points.

    Positive radii are concave (bending away from the symmetry axis), negative
    radii are convex, and straight portions are represented by ``np.inf``.
    """
    if axis != "x":
        raise NotImplementedError("Only a cylindrical x-axis is supported for now.")
    return _radius_from_curvature(signed_curvature(points, eps=eps), eps=eps)
