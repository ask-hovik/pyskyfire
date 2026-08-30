from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from pyskyfire.regen.cross_section import SectionProfiles
from pyskyfire.skycea.coolant_transport import CoolantTransport


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

    Attributes
    ----------
    material : Material
        Thermal material definition.
    _thickness : float or callable
        Underlying storage for the thickness definition.
    name : str or None
        Optional descriptive name.

    """
    def __init__(self, material, thickness, name=None):

        self.name = name
        self.material = material
        self._thickness = thickness

        assert type(thickness) is float or type(thickness) is int or callable(thickness), "'thickness' input must be a float, int or callable"

    def thickness(self, x):
        # x can be float or ndarray
        if callable(self._thickness):
            return np.asarray(self._thickness(x), dtype=float)
        # if x is an array, return a same-length array filled with the constant
        x_arr = np.asarray(x)
        if x_arr.ndim:
            return np.full_like(x_arr, float(self._thickness), dtype=float)
        return float(self._thickness)


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
    coolant_transport: CoolantTransport
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
    """

    def __init__(
        self,
        name,
        contour,
        cross_section,             # ChannelSection
        span: list[float],
        placement,                 # ChannelPlacement-like (sim uses only counts/lanes if needed)
        channel_height,            # callable x -> h
        walls,                     # list of wall objects with .thickness(x)
        coolant_transport,
        roughness,
        hot_gas_surface_area_multiplier=1.0,
    ):
        self.name = name
        self.contour = contour
        self.cross_section = cross_section
        self.placement = placement
        self.channel_height = channel_height
        self.coolant_transport = coolant_transport
        self.walls = walls
        self._roughness = roughness
        self._hot_gas_surface_area_multiplier = hot_gas_surface_area_multiplier
        if not callable(hot_gas_surface_area_multiplier):
            multiplier = float(hot_gas_surface_area_multiplier)
            if not np.isfinite(multiplier) or multiplier <= 0.0:
                raise ValueError(
                    "hot_gas_surface_area_multiplier must be finite and positive"
                )
            self._hot_gas_surface_area_multiplier = multiplier

        # fin/rib thermal enhancement toggle (default: disabled)
        self.enable_fin = getattr(self, "enable_fin", True)


        if span[0] > span[1]:
            self.span = (span[1], span[0])
            self.direction = -1
        else:
            self.span = span
            self.direction = 1

        # to be set later by set_centerline / set_channel_* / set_x_domain
        self.centerlines = None
        self.coolant_centerlines = None
        self.x_domain = None
        self.channel_heights = None
        self.channel_width = None
        self.helix_angle_vals = None
        self.is_helical = False


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

    def roughness(self, x):
        return self._roughness(x) if callable(self._roughness) else float(self._roughness)

    def hot_gas_surface_area_multiplier(self, x):
        """Return the circuit-specific effective hot-side area multiplier."""
        coordinates = np.asarray(x, dtype=float)
        modifier = self._hot_gas_surface_area_multiplier
        if callable(modifier):
            values = np.asarray(
                [modifier(float(value)) for value in coordinates.flat],
                dtype=float,
            ).reshape(coordinates.shape)
        else:
            values = np.full_like(coordinates, modifier, dtype=float)

        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError(
                f"{self.name}: hot_gas_surface_area_multiplier must return "
                "finite, positive values"
            )
        return float(values) if values.ndim == 0 else values

    def total_thickness(self, x):
        parts = [wall.thickness(x) for wall in self.walls]
        return np.sum(parts, axis=0, dtype=float)

    def channel_wall_thickness(self, x):
        """Return the coolant-adjacent layer used for channel geometry."""
        return self.walls[-1].thickness(x)

    def hot_side_layer_thickness(self, x):
        """Return liner/coating thickness between hot contour and tube wall."""
        if len(self.walls) == 1:
            values = np.asarray(x, dtype=float)
            result = np.zeros_like(values, dtype=float)
            return float(result) if result.ndim == 0 else result
        parts = [wall.thickness(x) for wall in self.walls[:-1]]
        return np.sum(parts, axis=0, dtype=float)

    # ---------------- profiles bundle ----------------
    def _prof(self, centerline: np.ndarray) -> SectionProfiles:
        """
        Build a SectionProfiles with trivial frames (analytics don't use them).
        """
        N = centerline.shape[0]
        x = centerline[:, 0]

        # blockage handling
        #br = getattr(self, "blockage_ratio", None)
        #if br is None:
        #    br_arr = np.full(N, 0.5, dtype=float)   # your previous default
        #else:
        #    br = np.asarray(br, dtype=float)
        #    br_arr = np.full(N, float(br), dtype=float) if br.ndim == 0 else br
        #    if br_arr.shape[0] != N:
         #       raise ValueError(f"blockage_ratio length {br_arr.shape[0]} != N {N}")

        # Cross-section arc/ligament geometry belongs to the layer touching
        # the coolant. Hot-side liners move the tube-facing reference surface
        # outward but are not themselves tube-wall ligament material.
        t_channel_wall = self.channel_wall_thickness(x)

        return SectionProfiles(
            h=np.asarray(self.channel_heights, float),
            theta=np.asarray(self.channel_local_sector, float),
            #theta=np.asarray(self.channel_width, float),
            t_wall=np.asarray(t_channel_wall, float),
            centerline=np.asarray(centerline, float),
            #blockage_ratio=br_arr,
        )

    # ---------------- preprocessing ----------------
    def precompute_thermal_properties(self):
        """
        Precompute passage geometry along a representative coolant path.

        ``centerlines`` remains the section-reference curve expected by the
        cross-section models.  ``coolant_centerlines`` is the actual modeled
        flow-centre curve and controls passage length and curvature.
        """
        if self.centerlines is None:
            raise RuntimeError("centerlines not set")
        section_centerline = self.centerlines[0]
        flow_centerline = (
            self.coolant_centerlines[0]
            if self.coolant_centerlines is not None
            else section_centerline
        )
        x_vals = section_centerline[:, 0]

        self.ds_dx_vals = cylindrical_curve_arc_length_scale(flow_centerline)

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
        """
        Single-station SectionProfiles for use in local closures (rib/fin, LUTs, etc.).
        Returns arrays of length 1.
        """
        # representative centerline
        cl = self.centerlines[0]
        xs = cl[:, 0]

        # interpolate centerline position
        r = float(np.interp(x, xs, cl[:, 1]))
        th = float(np.interp(x, xs, cl[:, 2]))
        centerline = np.array([[float(x), r, th]], dtype=float)

        # interpolate section scalars
        h = float(np.interp(x, self.x_domain, self.channel_heights))
        theta = float(np.interp(x, self.x_domain, self.channel_local_sector))
        t_wall = float(
            np.interp(
                x,
                self.x_domain,
                self.channel_wall_thickness(self.x_domain),
            )
        )

        return SectionProfiles(
            h=np.array([h], dtype=float),
            theta=np.array([theta], dtype=float),
            t_wall=np.array([t_wall], dtype=float),
            centerline=centerline,
        )


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
            #print(f"disabled: {var}")
            return var


        # Get 1-point profiles at x (needed by cross-section closure)
        prof_x = self.section_profiles_at(x)

        # Representative wall conductivity consistent with your wall stack usage)
        # Uses material next to coolant as material. Beware if added nickel coating is used TODO: implement a callout for which meterial is used as rib. 
        k_wall = float(self.walls[-1].material.get_k(T_wall_rep))

        # Cross-section provides R per unit channel length
        R_s = float(self.cross_section.R_coolant_per_len(prof_x, np.array([h_c]), k_wall)[0])

        # Convert to per unit axial length using ds/dx from your already-defined geometry
        ds_dx = float(self.ds_dx(x))
        ds_dx = max(ds_dx, 1e-30)

        var = R_s / ds_dx
        #print(f"enabled: {var}")
        return var


    def compute_volume(self):
        """Set ``self.volume``: the wetted coolant volume of the circuit [m³].

        One channel encloses :math:`\\int A\\,ds` along its own centreline, so
        the circuit holds :math:`n_{chan}\\int A\\,(ds/dx)\\,dx`.

        ``A_coolant`` is the flow area in the plane normal to the channel, not
        an axial slice: :meth:`ThrustChamber.build_channel_widths` narrows the
        sector by :math:`\\cos\\gamma` for a helix angle :math:`\\gamma`, while
        ``ds/dx`` lengthens the path by the reciprocal on a cylinder. A helical
        circuit therefore encloses the same volume as an axial one of the same
        azimuthal pitch, which is the physically correct answer -- only the
        path is longer and the section correspondingly narrower.

        Trapezoidal on the circuit's own ``x_domain``, so the result converges
        with axial resolution rather than being exact at any given ``n_nodes``.
        """
        x_vals = np.asarray(self.x_domain, dtype=float)
        volume_per_channel = np.trapezoid(
            self.A_coolant_vals * self.ds_dx_vals,
            x_vals,
        )
        # Volume is signed by the integration direction; the grid is built
        # ascending today, but do not depend on that.
        self.volume = abs(float(volume_per_channel)) * self.n_channels

    # ---------------- setters ----------------
    def set_centerline(self, centerline_list):
        """
        Provide one or more centerlines as arrays of shape (N, 3): [x, r, theta].
        For simulation we use the first as the representative path.
        """
        self.centerlines = [np.asarray(cl, float) for cl in centerline_list]

    def set_coolant_centerline(self, centerline_list):
        """Set modeled coolant-flow-centre curves ``[x, r, theta]``."""
        self.coolant_centerlines = [
            np.asarray(cl, float) for cl in centerline_list
        ]

    def set_channel_width(self, widths_rad):
        self.channel_width = np.asarray(widths_rad, float)

    def set_channel_local_sector(self, local_sectors):
        self.channel_local_sector = np.asarray(local_sectors, float)

    def set_channel_height(self, heights):
        self.channel_heights = np.asarray(heights, float)

    def set_x_domain(self, x_domain):
        self.x_domain = np.asarray(x_domain, float)

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
        return np.interp(x, self.x_domain, self.channel_width)
    
    def local_sector_angle(self, x): 
        return np.interp(x, self.x_domain, self.channel_local_sector)
    
    def ds_dx(self, x):
        return np.interp(x, self.x_domain, self.ds_dx_vals) # added for helical support
    

class ThrustChamber:
    """
    Simulation-only thrust-chamber assembly.

    Notes
    -----
    • Define x-domains for circuits from their fractional spans
    • Build section-reference and coolant-flow centerlines per circuit
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
        enable_fin: bool = True,
        film_cooling: FilmCooling | None = None,
    ):
        self.contour = contour
        self.cooling_circuits = cooling_circuits
        self.combustion_transport = combustion_transport
        self.optimal_values = optimal_values

        self.n_nodes   = int(n_nodes)
        self.K_factor  = float(K_factor)
        self.h_gas_corr  = float(h_gas_corr)
        self.h_cold_corr = float(h_cold_corr)
        

        self.enable_fin = bool(enable_fin)

        # propagate to circuits
        for c in self.cooling_circuits:
            c.enable_fin = self.enable_fin

        # --- geometry inputs for circuits (sim only) ---
        self.build_circuit_x_domain()
        self.build_channel_heights()
        self.build_channel_centerlines()
        self.build_channel_widths()

        # --- per-circuit preprocessing for solver ---
        for circuit in self.cooling_circuits:
            circuit.finalize()

        # --- hot-gas side (optional/on-demand) ---
        if compute_gas and hasattr(self.combustion_transport, "compute_aerothermodynamics"):
            self.combustion_transport.compute_aerothermodynamics(self.contour)


        self.film_cooling = film_cooling
        if self.film_cooling is not None:
            self.film_cooling.x = self._resolve_signed_fraction_to_x(self.film_cooling.x_fraction)

    # ------------------------------------------------------------------
    # Domain → centerline → section scalars (sim-only)
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
            frac  = L_seg / L_tot
            n_here = max(MIN_NODES, int(round(self.n_nodes * frac)))

            if x0 <= x1:
                x_domain = np.linspace(x0, x1, n_here)
            else:
                # reversed span: keep coolant marching direction by reversing
                x_domain = np.linspace(x1, x0, n_here)[::-1]

            circuit.set_x_domain(x_domain)

    def build_channel_centerlines(self):
        for circuit in self.cooling_circuits:
            xs = circuit.x_domain

            # The cross-section models use a wall/surface reference, whereas
            # length and curvature must follow the coolant centroid.
            section_radii = np.array(
                [
                    circuit.placement.compute_centerline_radius(x, self.contour)
                    for x in xs
                ],
                dtype=float,
            )
            section_radii += np.asarray(
                circuit.hot_side_layer_thickness(xs),
                dtype=float,
            )
            wall_thicknesses = np.asarray(
                circuit.total_thickness(xs),
                dtype=float,
            )
            flow_radii = np.array(
                [
                    circuit.placement.compute_flow_centerline_radius(
                        x,
                        self.contour,
                        wall_thickness,
                        channel_height,
                    )
                    for x, wall_thickness, channel_height in zip(
                        xs,
                        wall_thicknesses,
                        circuit.channel_heights,
                    )
                ],
                dtype=float,
            )
            flow_dr_dx = np.gradient(flow_radii, xs, edge_order=2)
            dtheta_dx = np.array(
                [
                    circuit.placement.centerline_dtheta_dx(
                        x,
                        self.contour,
                        flow_radius,
                        dr_dx,
                    )
                    for x, flow_radius, dr_dx in zip(
                        xs,
                        flow_radii,
                        flow_dr_dx,
                    )
                ],
                dtype=float,
            )

            # integrate θ(x): θ0=0 for representative lane
            theta = np.zeros_like(xs, dtype=float)
            theta[1:] = np.cumsum(
                0.5
                * (dtheta_dx[1:] + dtheta_dx[:-1])
                * (xs[1:] - xs[:-1])
            )

            if hasattr(circuit.placement, "helix_angle"):
                gamma = circuit.placement.helix_angle
                circuit.helix_angle_vals = np.array(
                    [gamma(x) if callable(gamma) else float(gamma) for x in xs],
                    dtype=float,
                )
            else:
                circuit.helix_angle_vals = np.zeros_like(xs)
            circuit.is_helical = bool(
                np.any(np.abs(circuit.helix_angle_vals) > 1e-12)
            )

            circuit.set_centerline(
                [np.column_stack([xs, section_radii, theta])]
            )
            circuit.set_coolant_centerline(
                [np.column_stack([xs, flow_radii, theta])]
            )


    def build_channel_widths(self):
        """
        Compute wedge angle (theta) arrays used by cross-section analytics.

        Priority:
        1) If placement provides channel_width(x), use it.
        2) If placement is 'internal' (non-occluding), use height-as-width (simple radial-stack proxy).
        3) Otherwise, default to an even packing among surface-occluding channels:
           theta(x) = 2π / n_occ(x), where n_occ is provided by the circuit group.
        """
        for circuit in self.cooling_circuits:
            xs = circuit.x_domain
            p  = circuit.placement

            if getattr(p, "channel_width", None) is not None:
                widths = np.array([p.channel_width(x) for x in xs], dtype=float)
                local_sector = widths.copy()

            elif getattr(p, "occludes", True) is False:
                # internal / non-occluding: use height as an angular proxy (sim only)
                widths = np.array([circuit.channel_height(x) for x in xs], dtype=float)
                local_sector = widths.copy()

            else:
                # uniform angular share among all occluding channels present at x
                n_occ = np.array(
                    [self.number_of_channels(x, occluding_only=True) for x in xs],
                    dtype=float
                )
                gamma = np.array(
                    [
                        p.helix_angle(x)
                        if callable(p.helix_angle)
                        else float(p.helix_angle)
                        for x in xs
                    ],
                    dtype=float,
                )
                #print
                widths = np.where(n_occ > 0.0, 2.0 * np.pi / n_occ, 0.0) 
                local_sector = np.where(n_occ > 0.0, 2.0 * np.pi * np.cos(gamma) / n_occ, 0.0) # multiplying with cos gamma for helical implementation

            circuit.set_channel_width(widths)
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
                #total_channels += circuit.placement.n_channel_positions
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
    xyz = np.column_stack(
        [x, radius * np.cos(theta), radius * np.sin(theta)]
    )

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
        2.0
        * np.linalg.norm(np.cross(chord_before, chord_after), axis=1)
        / denominator
    )

    # The change in unit chord direction points toward the local centre of
    # curvature.  Only its radial projection is needed for the RL10 sign.
    tangent_change = (
        chord_after / length_after[:, None]
        - chord_before / length_before[:, None]
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
    radial_component = magnitude * np.einsum(
        "ij,ij->i", principal_normal, radial_unit
    )
    direction = np.sign(radial_component)

    # A locally axial curvature vector can have a vanishing radial projection.
    # Retain the established meridional sign convention as a fallback.
    slope_before = (
        radius[middle] - radius[middle - 1]
    ) / (
        x[middle] - x[middle - 1]
    )
    slope_after = (
        radius[middle + 1] - radius[middle]
    ) / (
        x[middle + 1] - x[middle]
    )
    d2r_dx2 = 2.0 * (slope_after - slope_before) / (
        x[middle + 1] - x[middle - 1]
    )
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


