from __future__ import annotations

import numpy as np
from math import gcd
from functools import reduce
from dataclasses import dataclass
from typing import Optional

from pyskyfire.regen.cross_section import SectionProfiles
from pyskyfire.regen.channel_placement import SurfacePlacement, InternalPlacement
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

    Notes
    -----
    - No OCC / geometry generation here.
    - Does NOT compute true local frames; those are a visualization concern.
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
    ):
        self.name = name
        self.contour = contour
        self.cross_section = cross_section
        self.placement = placement
        self.channel_height = channel_height
        self.coolant_transport = coolant_transport
        self.walls = walls
        self._roughness = roughness

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
        self.x_domain = None
        self.channel_heights = None
        self.channel_width = None


    # ---------------- basic getters ----------------
    def roughness(self, x):
        return self._roughness(x) if callable(self._roughness) else float(self._roughness)

    def total_thickness(self, x):
        parts = [wall.thickness(x) for wall in self.walls]
        return np.sum(parts, axis=0, dtype=float)

    # ---------------- profiles bundle ----------------
    def _prof(self, centerline: np.ndarray) -> "SectionProfiles":
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

        t_wall_tot = self.total_thickness(x)

        return SectionProfiles(
            h=np.asarray(self.channel_heights, float),
            theta=np.asarray(self.channel_local_sector, float),
            #theta=np.asarray(self.channel_width, float),
            t_wall=np.asarray(t_wall_tot, float),
            centerline=np.asarray(centerline, float),
            #blockage_ratio=br_arr,
        )

    # ---------------- preprocessing ----------------
    def precompute_thermal_properties(self):
        """
        Precompute A, Dh, and thermal perimeters along a representative centerline.
        """
        if self.centerlines is None:
            raise RuntimeError("centerlines not set")
        centerline = self.centerlines[0]
        x_vals = centerline[:, 0]
        r_vals = centerline[:, 1]
        th_vals = centerline[:, 2]

        dr_dx     = np.gradient(r_vals, x_vals)
        dtheta_dx = np.gradient(th_vals, x_vals)
        # arc-length scale factor along the space-curve in cylindrical coords
        ds_dx = np.sqrt(1.0 + dr_dx**2 + (r_vals * dtheta_dx)**2) # TODO: revisit this and think about whether it is correct still or not
        
        self.ds_dx_vals = ds_dx # added as part of the helical stuff

        prof = self._prof(centerline)

        hot_perimeter  = self.cross_section.P_thermal(prof)
        cold_perimeter = self.cross_section.P_coolant(prof)

        self.dA_dx_thermal_exhaust_vals = hot_perimeter  * ds_dx
        self.dA_dx_thermal_coolant_vals = cold_perimeter * ds_dx

        A_vals = self.cross_section.A_coolant(prof)
        self.A_coolant_vals = A_vals
        self.dA_dx_coolant_vals = np.gradient(A_vals, x_vals)

        self.Dh_coolant_vals = self.cross_section.Dh_coolant(prof)

        # simple planar curvature proxy in meridional plane (optional)
        # :math:`\kappa` ≈ |r''| / (1 + r'^2)^(3/2)
        r_pp = np.gradient(dr_dx, x_vals)
        self.radius_of_curvature_vals = 1.0 / np.clip(
            np.abs(r_pp) / (1.0 + dr_dx**2) ** 1.5, 1e-12, None
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
        t_wall = float(np.interp(x, self.x_domain, self.total_thickness(self.x_domain)))

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

        # Representative wall conductivity (use liner/wall[0]; consistent with your wall stack usage)
        k_wall = float(self.walls[0].material.get_k(T_wall_rep))

        # Cross-section provides R per unit channel length
        R_s = float(self.cross_section.R_coolant_per_len(prof_x, np.array([h_c]), k_wall)[0])

        # Convert to per unit axial length using ds/dx from your already-defined geometry
        ds_dx = float(self.ds_dx(x))
        ds_dx = max(ds_dx, 1e-30)

        var = R_s / ds_dx
        #print(f"enabled: {var}")
        return var


    def compute_volume(self):
        """
        Total circuit volume (all physical lanes) via ∫ A ds, scaled by lane count.
        """
        centerline = self.centerlines[0]
        x_vals = centerline[:, 0]
        r_vals = centerline[:, 1]
        th_vals = centerline[:, 2]

        dr_dx     = np.gradient(r_vals, x_vals)
        dtheta_dx = np.gradient(th_vals, x_vals)
        ds_dx = np.sqrt(1.0 + dr_dx**2 + (r_vals * dtheta_dx)**2)

        volume_per_lane = np.trapezoid(self.A_coolant_vals * ds_dx, x_vals)

        # if placement exposes n_channel_positions, scale; else assume 1
        n_lanes = getattr(self.placement, "n_channel_positions", 1)
        self.volume = float(volume_per_lane) * int(n_lanes)

    # ---------------- setters ----------------
    def set_centerline(self, centerline_list):
        """
        Provide one or more centerlines as arrays of shape (N, 3): [x, r, theta].
        For simulation we use the first as the representative path.
        """
        self.centerlines = [np.asarray(cl, float) for cl in centerline_list]

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
        return np.interp(x, self.x_domain, self.radius_of_curvature_vals)

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
    • Build ONE representative centerline per circuit: [x, r(x), θ_azimuth=0]
    • Compute per-circuit wedge angles (theta) and heights (h)
    • Trigger each circuit's precompute/finalize (A, Dh, perimeters, volume)
    • (Optionally) trigger combustion aerothermodynamics

    No visualization, no OCC/GMsh imports, no multiple angular copies.
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
        self.build_channel_centerlines()
        self.build_channel_widths()
        self.build_channel_heights()

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

            rs = np.array([circuit.placement.compute_centerline_radius(x, self.contour) for x in xs], float)
            dth = np.array([circuit.placement.dtheta_dx(x, self.contour) for x in xs], float)

            # integrate θ(x): θ0=0 for representative lane
            theta = np.zeros_like(xs, dtype=float)
            theta[1:] = np.cumsum(0.5 * (dth[1:] + dth[:-1]) * (xs[1:] - xs[:-1]))

            centerline = np.column_stack([xs, rs, theta])
            circuit.set_centerline([centerline])


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
                gamma = np.array([p.helix_angle(x) for x in xs], dtype=float)
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






def radius_of_curvature(
    points: np.ndarray,
    axis: str = "x",
    eps: float = 1e-12,
) -> np.ndarray:
    """
    Signed radius of curvature for a curve expressed in cylindrical coordinates
    [x, r, θ].

    • Positive  → curve bends *away* from the symmetry axis  
    • Negative  → curve bends *toward* the symmetry axis  
    • np.inf    → locally straight (|:math:`\kappa`| below `eps`)

    Parameters
    ----------
    points : (N, 3) ndarray
        [[x, r, theta], …] ordered along the curve.
    axis   : {'x', 'y', 'z'}, optional
        Which coordinate is the symmetry axis.  Default 'x'.
    eps    : float, optional
        Curvature values with |:math:`\kappa`| < eps are treated as zero (straight).

    Returns
    -------
    R : (N,) ndarray
        Signed radius of curvature at each sample.
    """
    # -------- 1. unpack and move the chosen axis to coordinate 0 ---------- #
    # cylindrical input is always (x, r, θ) → (axis, radial, θ)
    if axis != "x":
        raise NotImplementedError("Only a cylindrical x-axis is supported for now.")

    x = points[:, 0]           # axial coordinate (monotonic ordering recommended)
    r = points[:, 1]           # radial distance

    # -------- 2. first & second derivatives w.r.t. the axis -------------- #
    dr_dx   = np.gradient(r, x)         #   r′(x)
    d2r_dx2 = np.gradient(dr_dx, x)     #   r″(x)

    # -------- 3. curvature (:math:`\kappa`) and signed radius (R) --------------------- #
    kappa = d2r_dx2 / np.power(1.0 + dr_dx**2, 1.5)   # :math:`\kappa` = r″ / (1+r′²)³ᐟ²

    # treat tiny :math:`\kappa` as straight line → infinite radius
    with np.errstate(divide="ignore"):
        R = np.where(np.abs(kappa) < eps, np.inf, 1.0 / kappa)

    return R


