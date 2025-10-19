from __future__ import annotations
from pyskyfire.regen.cross_section import SectionProfiles
import numpy as np
from abc import ABC, abstractmethod
from math import gcd
from functools import reduce

class Contour:
    """Inner hot-gas contour of a bell-type rocket engine.

    Represents the geometric wall line from the start of the chamber to
    the nozzle exit, providing local areas, slopes, and radii.

    Parameters
    ----------
    xs : array_like
        Axial coordinates of the contour [m], strictly increasing.
    rs : array_like
        Corresponding wall radii [m].
    name : str, optional
        Identifier for this contour.

    Attributes
    ----------
    xs : ndarray
        Axial coordinates defining the wall line [m].
    rs : ndarray
        Radial coordinates corresponding to `xs` [m].
    _dr_dx : ndarray
        Precomputed derivative `dr/dx` for fast interpolation.
    name : str or None
        Optional descriptive name.
    x_t, r_t, A_t : float
        Axial position, radius, and area of the throat.
    r_e, A_e : float
        Radius and area at the exit plane.
    r_c, A_c : float
        Radius and area at the chamber start.
    eps, eps_c : float
        Nozzle and contraction area ratios.

    Notes
    -----
    This class is used by higher-level objects such as
    :class:`~pyskyfire.regen.cooling.CoolingCircuit` and
    :class:`~pyskyfire.regen.thrust.ThrustChamber`.

    See Also
    --------
    ContourToroidalAerospike : Dual-wall variant for toroidal aerospikes.
    """

    def __init__(self, xs, rs, name = None):

        self.xs = xs
        self.rs = rs
        self._dr_dx = np.gradient(rs, xs)
        self.name = name

    def __setattr__(self, name, value):
        # If the user tries to set 'xs' or 'rs', we need to recalculate self.dr_dx
        if name == "xs" and hasattr(self, "rs"):
            self._dr_dx = np.gradient(self.rs, value)

        elif name == "rs" and hasattr(self, "xs"):
            self._dr_dx = np.gradient(value, self.xs)

        super(Contour, self).__setattr__(name, value)
    
    @property
    def x_t(self):
        return self.xs[np.argmin(self.rs)]

    @property
    def r_t(self):
        return min(self.rs)

    @property
    def A_t(self):
        return np.pi * self.r_t**2

    @property
    def r_e(self):
        return self.rs[-1]
    
    @property
    def r_c(self):
        return self.rs[0]
    
    @property
    def A_e(self):
        return np.pi * self.r_e**2
    
    @property
    def A_c(self):
        return np.pi * self.r_c**2
    
    @property
    def eps(self):
        return self.A_e/self.A_t

    @property
    def eps_c(self): 
        """ Get the contraction ratio for the chamber
        Returns: 
            (float): contraction ratio """
        return self.A_c/self.A_t

    def r(self, x):
        """
        Return the local radius at axial position `x`.

        Parameters
        ----------
        x : float
            Axial coordinate [m].

        Returns
        -------
        float
            Distance from engine centerline to wall [m].
        """
        r1 = np.interp(x, self.xs, self.rs)
        #print(f"r at that point: {r1}, input x: {x} ")
        return r1
    
    def dr_dx(self, x):
        """
        Return the local slope of the wall, :math:`dr/dx`.

        Parameters
        ----------
        x : float
            Axial coordinate [m].

        Returns
        -------
        float
            Radial slope at position `x`.
        """
        return np.interp(x, self.xs, self._dr_dx)

    def A(self, x):
        """
        Return the local flow area.

        Parameters
        ----------
        x : float
            Axial coordinate [m].

        Returns
        -------
        float
            Flow area at that section [m²].
        """
        Area = np.pi * self.r(x)**2
        #print(f"Area: {Area}")
        return Area
    
    def normal_angle(self, x): #TODO: check this function af
        """
        Return the local wall normal angle with respect to the vertical plane.

        Parameters
        ----------
        x : float
            Axial coordinate [m].

        Returns
        -------
        float
            Angle between the outward normal and the plane normal to the x-axis [rad].
        """
        slope = self.dr_dx(x)
        # Dot product: N • V = 1
        # |N| = sqrt(slope^2 + 1), |V| = 1
        cos_angle = 1.0 / np.sqrt(slope**2 + 1.0)
        # Numerically clamp in case of tiny floating errors
        cos_angle = max(-1.0, min(1.0, cos_angle))

        angle = np.arccos(cos_angle)
        # arccos(...) is already in [0..pi]. Because slope^2 >= 0,
        # cos_angle is in (0..1], so angle is in [0..pi/2].
        return angle
    

class ContourToroidalAerospike:
    """Axisymmetric contour describing a toroidal-aerospike geometry.

    Defines both the inner and outer wall surfaces, allowing evaluation of
    annular areas and local geometric derivatives.

    Parameters
    ----------
    xs_outer, rs_outer : array_like
        Axial and radial coordinates of the outer wall [m].
    xs_inner, rs_inner : array_like
        Axial and radial coordinates of the inner wall [m].
    name : str, optional
        Identifier for this contour.

    Attributes
    ----------
    xs_outer, rs_outer : ndarray
        Outer wall geometry arrays [m].
    xs_inner, rs_inner : ndarray
        Inner wall geometry arrays [m].
    _dr_dx_outer, _dr_dx_inner : ndarray
        Local wall slope arrays for each surface.
    name : str or None
        Descriptive label.
    A_t, A_c, A_e : float
        Annular areas at throat, chamber, and exit [m²].
    eps, eps_c : float
        Expansion and contraction ratios (dimensionless).

    Notes
    -----
    The inner radius is validated to always remain below the outer radius.

    See Also
    --------
    Contour : Single-wall version used in standard bell-nozzle engines.
    """

    # -------------------------------------------------------------------------
    # Constructor & validation
    # -------------------------------------------------------------------------
    def __init__(self, xs_outer, rs_outer, xs_inner, rs_inner, *, name=None):
        # -- convert & basic checks --------------------------------------------------
        self.xs_outer = np.asarray(xs_outer, dtype=float)
        self.rs_outer = np.asarray(rs_outer, dtype=float)
        self.xs_inner = np.asarray(xs_inner, dtype=float)
        self.rs_inner = np.asarray(rs_inner, dtype=float)

        for tag, xs in ("xs_outer", self.xs_outer), ("xs_inner", self.xs_inner):
            if xs.ndim != 1:
                raise ValueError(f"{tag} must be 1-D")
            if xs.size < 2:
                raise ValueError(f"{tag} needs at least two points")
            if not np.all(np.diff(xs) > 0):
                raise ValueError(f"{tag} must be strictly increasing")

        if self.xs_outer.shape != self.rs_outer.shape:
            raise ValueError("xs_outer and rs_outer must have the same length")
        if self.xs_inner.shape != self.rs_inner.shape:
            raise ValueError("xs_inner and rs_inner must have the same length")
        if np.any(self.rs_inner > np.interp(self.xs_inner, self.xs_outer, self.rs_outer, left=np.inf, right=np.inf)):
            raise ValueError("Inner radius exceeds outer radius somewhere")

        # -- derivatives along each wall -------------------------------------------
        self._dr_dx_outer = np.gradient(self.rs_outer, self.xs_outer)
        self._dr_dx_inner = np.gradient(self.rs_inner, self.xs_inner)

        self.name = name  # last so that __setattr__ doesn’t re-enter validation

    # -------------------------------------------------------------------------
    # Helpers – interpolation on each wall
    # -------------------------------------------------------------------------
    def _interp_outer(self, x):
        """Linear interpolation of outer radius at *x*."""
        return np.interp(x, self.xs_outer, self.rs_outer)

    def _interp_inner(self, x):
        """Linear interpolation of inner radius at *x*."""
        return np.interp(x, self.xs_inner, self.rs_inner)

    # -------------------------------------------------------------------------
    # Outward-facing geometric API (outer wall when ambiguous)
    # -------------------------------------------------------------------------
    @property
    def x_t(self):
        """Axial location of the **outer** throat (m)."""
        return self.xs_outer[np.argmin(self.rs_outer)]

    @property
    def r_t(self):
        """Throat radius of the **outer** wall (m)."""
        return np.min(self.rs_outer)

    @property
    def r_e(self):
        """Outer radius at exit (m)."""
        return self.rs_outer[-1]

    @property
    def r_c(self):
        """Outer radius at chamber start (m)."""
        return self.rs_outer[0]

    # -------------------------------------------------------------------------
    # Areas (annular)
    # -------------------------------------------------------------------------
    def A(self, x):
        """
        Return the local annular flow area.

        Parameters
        ----------
        x : float
            Axial coordinate [m].

        Returns
        -------
        float
            Annular cross-sectional area [m²].
        """
        r_o = self._interp_outer(x)
        r_i = self._interp_inner(x)
        return np.pi * (r_o**2 - r_i**2)

    @property
    def A_t(self):
        """Annular area at the throat (m²)."""
        idx = np.argmin(self.rs_outer)
        r_o = self.rs_outer[idx]
        # need r_i at *same axial position* of throat (outer x). interpolate inner
        x_throat = self.xs_outer[idx]
        r_i = self._interp_inner(x_throat)
        return np.pi * (r_o**2 - r_i**2)

    @property
    def A_e(self):
        r_o = self.r_e
        # inner radius at outer-contour exit station (assume xs_outer[-1])
        r_i = self._interp_inner(self.xs_outer[-1])
        return np.pi * (r_o**2 - r_i**2)

    @property
    def A_c(self):
        r_o = self.r_c
        r_i = self._interp_inner(self.xs_outer[0])
        return np.pi * (r_o**2 - r_i**2)

    # -------------------------------------------------------------------------
    # Ratios
    # -------------------------------------------------------------------------
    @property
    def eps(self):
        """Area ratio exit / throat."""
        return self.A_e / self.A_t

    @property
    def eps_c(self):
        """Chamber contraction ratio."""
        return self.A_c / self.A_t

    # -------------------------------------------------------------------------
    # Slopes & normals
    # -------------------------------------------------------------------------
    def dr_dx(self, x, which="outer"):
        """
        Return the wall slope :math:`dr/dx` at `x` for either wall.

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        which : {'outer', 'inner'}, optional
            Which wall to evaluate. Default is 'outer'.

        Returns
        -------
        float
            Local slope of the selected wall.
        """
        if which == "outer":
            return np.interp(x, self.xs_outer, self._dr_dx_outer)
        elif which == "inner":
            return np.interp(x, self.xs_inner, self._dr_dx_inner)
        else:
            raise ValueError("which must be 'outer' or 'inner'")

    def normal_angle(self, x, which="outer"):
        """
        Return the normal-plane angle for either wall.

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        which : {'outer', 'inner'}, optional
            Wall selector. Default 'outer'.

        Returns
        -------
        float
            Angle between outward normal and vertical plane [rad].
        """
        slope = self.dr_dx(x, which=which)
        cos_ang = 1.0 / np.sqrt(slope**2 + 1.0)
        return np.arccos(np.clip(cos_ang, -1.0, 1.0))

    # -------------------------------------------------------------------------
    # Compatibility helpers / aliases (outer wall)
    # -------------------------------------------------------------------------
    def r(self, x):
        """Return **outer** radius at *x* (m). Provided for API compatibility."""
        return self._interp_outer(x)

    # -------------------------------------------------------------------------
    # Self-updating gradients if arrays mutate after construction
    # -------------------------------------------------------------------------
    def __setattr__(self, key, value):
        # Keep gradients coherent when arrays are replaced
        if key == "xs_outer":
            object.__setattr__(self, key, value)
            if hasattr(self, "rs_outer"):
                self._dr_dx_outer = np.gradient(self.rs_outer, value)
            return
        if key == "rs_outer":
            object.__setattr__(self, key, value)
            if hasattr(self, "xs_outer"):
                self._dr_dx_outer = np.gradient(value, self.xs_outer)
            return
        if key == "xs_inner":
            object.__setattr__(self, key, value)
            if hasattr(self, "rs_inner"):
                self._dr_dx_inner = np.gradient(self.rs_inner, value)
            return
        if key == "rs_inner":
            object.__setattr__(self, key, value)
            if hasattr(self, "xs_inner"):
                self._dr_dx_inner = np.gradient(value, self.xs_inner)
            return

        super().__setattr__(key, value)

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

    See Also
    --------
    WallGroup : Container combining multiple walls.
    """
    def __init__(self, material, thickness, name=None):

        self.name = name
        self.material = material
        self._thickness = thickness

        assert type(thickness) is float or type(thickness) is int or callable(thickness), "'thickness' input must be a float, int or callable"

    def thickness(self, x):
        """
        Return the wall thickness at position `x`.

        Parameters
        ----------
        x : float
            Axial coordinate [m].

        Returns
        -------
        float
            Wall thickness [m].
        """
        if callable(self._thickness):
            return self._thickness(x)

        else:
            return self._thickness

class WallGroup:
    """Container holding multiple :class:`Wall` layers.

    Provides cumulative quantities such as total wall thickness along the
    engine contour.

    Parameters
    ----------
    walls : list[Wall], optional
        Collection of wall layers. Defaults to an empty list.

    Attributes
    ----------
    walls : list[Wall]
        Managed list of wall objects.

    See Also
    --------
    ThrustChamber : Uses a `WallGroup` to compute wall-conduction resistance.
    """
    def __init__(self, walls=None):
        if walls is None:
            walls = []
        self.walls = walls

    def total_thickness(self, x):
        """
        Compute total wall thickness at axial position `x`.

        Parameters
        ----------
        x : float
            Axial coordinate [m].

        Returns
        -------
        float
            Sum of all wall thicknesses [m].
        """
        return sum(wall.thickness(x) for wall in self.walls)


class CoolingCircuit:
    """Representation of a cooling circuit following the chamber contour.

    Each circuit defines the geometry of one group of coolant channels,
    including their cross-sectional shape, span, and placement relative to
    the hot wall.

    Parameters
    ----------
    name : str
        Identifier of the circuit.
    contour : Contour
        Hot-gas wall geometry defining the outer boundary.
    cross_section : ChannelSection
        Cross-sectional geometry model providing A, Dh, and perimeters.
    span : tuple[float, float]
        Normalized start and end of circuit (-1 = chamber inlet, +1 = nozzle exit).
    placement : ChannelPlacement
        Strategy describing how channel centerlines are positioned.
    channel_height : callable
        Function returning the local channel height [m].
    coolant_transport : object
        Object providing coolant thermophysical properties.
    blockage_ratio : float or array_like, optional
        Fraction of the sector blocked by solid wall or rib.

    Attributes
    ----------
    name : str
        Circuit name.
    contour : Contour
        Reference hot-gas contour.
    cross_section : ChannelSection
        Shape model used for thermal and hydraulic quantities.
    placement : ChannelPlacement
        Placement rule describing angular and radial positioning.
    channel_height : callable
        Function returning height as a function of x.
    coolant_transport : object
        Provides coolant properties (k, μ, Cp, ρ, etc.).
    span : (float, float)
        Normalized start/end bounds.
    direction : int
        +1 if flow is forward (increasing x), −1 if reversed.
    A_coolant_vals, Dh_coolant_vals : ndarray
        Precomputed local coolant area and hydraulic diameter.
    dA_dx_thermal_exhaust_vals, dA_dx_thermal_coolant_vals : ndarray
        Surface-area differentials on hot and cold sides.
    radius_of_curvature_vals : ndarray
        Local curvature radius for each axial node.

    Notes
    -----
    The circuit precomputes local geometric properties for fast interpolation
    during steady-state or transient simulations.

    See Also
    --------
    CoolingCircuitGroup : Groups multiple cooling circuits.
    SectionProfiles : Bundles local geometry inputs for cross-section methods.
    """
    def __init__(self, name, contour, cross_section, span, placement, channel_height, coolant_transport, blockage_ratio=None): 
        self.name = name
        self.contour = contour
        self.cross_section = cross_section
        self.placement = placement
        self.channel_height = channel_height
        self.coolant_transport = coolant_transport
        self.blockage_ratio = blockage_ratio

        if span[0] > span[1]:
            self.span = [span[1], span[0]]
            self.direction = -1
        else:
            self.span = span
            self.direction = 1

    # --- minimal helper: wrap the SAME params you used to pass before -----------
    def _prof(self, centerline, local_coords):
        """
        Assemble a `SectionProfiles` object for a given centerline.

        Parameters
        ----------
        centerline : ndarray, shape (N, 3)
            Channel centerline coordinates (x, r, θ).
        local_coords : ndarray, shape (N, 3, 3)
            Local coordinate frames.

        Returns
        -------
        SectionProfiles
            Ready-to-use profile bundle for cross-section routines.
        """
        N = centerline.shape[0]
        br = getattr(self, "blockage_ratio", None)

        if br is None:
            br_arr = np.full(N, 0.5, dtype=float)  # default preserves current behavior
        else:
            br = np.asarray(br, dtype=float)
            br_arr = np.full(N, float(br), dtype=float) if br.ndim == 0 else br
            if br_arr.shape[0] != N:
                raise ValueError(f"blockage_ratio length {br_arr.shape[0]} != N {N}")

        return SectionProfiles(
            h=np.asarray(self.channel_heights, float),
            theta=np.asarray(self.channel_width, float),
            t_wall=np.asarray(self.t_wall_tot, float),
            centerline=np.asarray(centerline, float),
            local_coords=np.asarray(local_coords, float),
            blockage_ratio=br_arr
        )

    def precompute_thermal_properties(self):
        """
        Precompute all cross-section-dependent thermal geometry arrays.

        Calculates effective surface-area derivatives, hydraulic diameter, and
        radius of curvature for interpolation during simulation.
        """
        centerline   = self.centerlines[0]
        local_coords = self.local_coords_list[0]
        x_vals       = centerline[:, 0]
        r_vals       = centerline[:, 1]

        dx_dx_val     = self.centerline_deriv_list[0][:, 0]
        dr_dx_val     = self.centerline_deriv_list[0][:, 1]
        dtheta_dx_val = self.centerline_deriv_list[0][:, 2]

        # keep your current simplified ds/dx (your TODO remains)
        ds_dx = np.sqrt(1.0 + dr_dx_val**2)

        # pass EXACTLY the same inputs as before, but wrapped in prof
        prof = self._prof(centerline, local_coords)

        hot_perimeter  = self.cross_section.P_thermal(prof)
        cold_perimeter = self.cross_section.P_coolant(prof)

        self.dA_dx_thermal_exhaust_vals = hot_perimeter  * ds_dx
        self.dA_dx_thermal_coolant_vals = cold_perimeter * ds_dx

        A_coolant_vals = self.cross_section.A_coolant(prof)
        self.A_coolant_vals = A_coolant_vals
        self.dA_dx_coolant_vals = np.gradient(A_coolant_vals, x_vals)

        self.Dh_coolant_vals = self.cross_section.Dh_coolant(prof)

        self.radius_of_curvature_vals = radius_of_curvature(centerline)

    def compute_volume(self):
        """
        Compute total circuit volume by integrating local area along its centerline.

        Returns
        -------
        float
            Total coolant volume [m³].
        """
        centerline = self.centerlines[0]
        x_vals = centerline[:, 0]
        r_vals = centerline[:, 1]

        dr_dx     = self.centerline_deriv_list[0][:, 1]
        dtheta_dx = self.centerline_deriv_list[0][:, 2]
        ds_dx = np.sqrt(1.0 + dr_dx**2 + (r_vals * dtheta_dx)**2)

        volume_per_channel = np.trapezoid(self.A_coolant_vals * ds_dx, x_vals)
        total_volume = volume_per_channel * self.placement.n_channel_positions
        self.volume = total_volume

    def compute_single_centerline(self):
        """
        Generate OCC wire objects for each station along the first centerline.

        Intended for CAD or meshing visualization.
        """
        list_of_wires = []
        centerl = self.centerlines[0]
        for i in range(len(centerl)):
            local_coords = self.local_coords_list[i]
            prof_i = self._prof(centerl, local_coords)
            wire = self.cross_section.compute_cross_section(prof_i, i)
            list_of_wires.append(wire)
        self.wires = list_of_wires

    def compute_geometry(self):
        """
        Generate full 3D point-cloud representations for all channel centerlines.

        Each point cloud corresponds to one physical cooling channel.
        """

        all_point_clouds = []
        for i, centerline in enumerate(self.centerlines):
            local_coords = self.local_coords_list[i]
            prof_i = self._prof(centerline, local_coords)
            # minimal change: pass just prof
            point_cloud = self.cross_section.compute_point_cloud(prof_i)
            all_point_clouds.append(point_cloud)
        self.point_cloud = all_point_clouds

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
    
    def set_centerline_test(self, centerline_list):
        self.centerlines = centerline_list
        self.local_coords_list = []
        self.centerline_deriv_list = []

        for centerline in centerline_list:
            x_vals = centerline[:, 0]
            r_vals = centerline[:, 1]
            theta_vals = centerline[:, 2]

            points_3d = np.column_stack((
                x_vals,
                r_vals * np.cos(theta_vals),
                r_vals * np.sin(theta_vals),
            ))

            tangent_vectors = np.zeros_like(points_3d)
            for i in range(3):
                tangent_vectors[:, i] = np.gradient(points_3d[:, i], x_vals)
            norms = np.linalg.norm(tangent_vectors, axis=1, keepdims=True)
            tangent_vectors = tangent_vectors / np.clip(norms, 1e-12, None)
        
            # Vector from each point to the nearest point on the x-axis (i.e., to (x, 0, 0))
            delta_to_axis = np.column_stack([
                np.zeros_like(points_3d[:, 0]),   # no x component
                -points_3d[:, 1],                 # -y
                -points_3d[:, 2],                 # -z
            ])

            # Project "delta_to_axis" onto the plane perpendicular to the tangent
            dot_dt = np.sum(delta_to_axis * tangent_vectors, axis=1, keepdims=True)
            n = delta_to_axis - dot_dt * tangent_vectors

            # Normalize; handle near-degenerate cases with sensible fallbacks
            n_norm = np.linalg.norm(n, axis=1, keepdims=True)
            tiny = 1e-12
            bad = (n_norm[:, 0] < 1e-10)  # where projection nearly vanished

            # Fallback 1: start from pure inward radial (yz) direction, then orthogonalize to t
            rad = np.column_stack([
                np.zeros_like(points_3d[:, 0]),
                -points_3d[:, 1],
                -points_3d[:, 2]
            ])
            rad_norm = np.linalg.norm(rad, axis=1, keepdims=True)
            rad_unit = np.divide(rad, np.clip(rad_norm, tiny, None))
            n_fb1 = rad_unit - np.sum(rad_unit * tangent_vectors, axis=1, keepdims=True) * tangent_vectors
            n_fb1_norm = np.linalg.norm(n_fb1, axis=1, keepdims=True)

            # Fallback 2: if still bad (e.g., point on axis and awkward tangent), use a double-cross with e_x
            ex = np.array([1.0, 0.0, 0.0])
            t_cross_ex = np.cross(tangent_vectors, ex)
            n_fb2 = np.cross(tangent_vectors, t_cross_ex)
            n_fb2_norm = np.linalg.norm(n_fb2, axis=1, keepdims=True)

            # Choose the best available normal at each point
            use_fb1 = bad & (n_fb1_norm[:, 0] >= 1e-10)
            use_fb2 = bad & ~use_fb1

            n[use_fb1] = n_fb1[use_fb1]
            n_norm[use_fb1] = n_fb1_norm[use_fb1]

            n[use_fb2] = n_fb2[use_fb2]
            n_norm[use_fb2] = n_fb2_norm[use_fb2]

            # Final normalize
            n = np.divide(n, np.clip(np.linalg.norm(n, axis=1, keepdims=True), tiny, None))

            # Ensure orientation truly points toward the x-axis
            # (flip if the angle to delta_to_axis is obtuse)
            flip = (np.sum(n * delta_to_axis, axis=1) < 0.0)
            n[flip] *= -1.0

            normal_vectors = n  # unit normals, perpendicular to tangents, pointing inward toward x-axis

            binormal_vectors = np.cross(tangent_vectors, normal_vectors)
            local_coords = np.stack((tangent_vectors, normal_vectors, binormal_vectors), axis=1)
            self.local_coords_list.append(local_coords)

            dr_dx = np.gradient(r_vals, x_vals)
            dtheta_dx = np.gradient(theta_vals, x_vals)
            centerline_deriv = np.column_stack((x_vals, dr_dx, dtheta_dx))
            self.centerline_deriv_list.append(centerline_deriv)
        # Example call (after computing tangent_vectors and normal_vectors for one centerline):
        #plot_local_frames(points_3d, tangent_vectors, normal_vectors, vec_len=0.05)


    def set_centerline(self, centerline_list):
        self.centerlines = centerline_list
        self.local_coords_list = []
        self.centerline_deriv_list = []

        for centerline in centerline_list:
            x_vals = centerline[:, 0]
            r_vals = centerline[:, 1]
            theta_vals = centerline[:, 2]

            points_3d = np.column_stack((
                x_vals,
                r_vals * np.cos(theta_vals),
                r_vals * np.sin(theta_vals),
            ))

            tangent_vectors = np.zeros_like(points_3d)
            for i in range(3):
                tangent_vectors[:, i] = np.gradient(points_3d[:, i], x_vals)
            norms = np.linalg.norm(tangent_vectors, axis=1, keepdims=True)
            tangent_vectors = tangent_vectors / np.clip(norms, 1e-12, None)

            normals = np.zeros_like(points_3d)
            for i, (P, t) in enumerate(zip(points_3d, tangent_vectors)):
                x, y, z = P
                t_x, t_y, t_z = t
                if abs(t_x) < 1e-6:
                    candidate = np.array([0, -y, -z])
                else:
                    candidate = np.array([(y * t_y + z * t_z) / t_x, -y, -z])
                nrm = np.linalg.norm(candidate)
                normals[i] = candidate / nrm if nrm > 1e-6 else np.array([0.0, 0.0, 0.0])

            binormal_vectors = np.cross(tangent_vectors, normals)
            local_coords = np.stack((tangent_vectors, normals, binormal_vectors), axis=1)
            self.local_coords_list.append(local_coords)

            dr_dx = np.gradient(r_vals, x_vals)
            dtheta_dx = np.gradient(theta_vals, x_vals)
            centerline_deriv = np.column_stack((x_vals, dr_dx, dtheta_dx))
            self.centerline_deriv_list.append(centerline_deriv)

    def set_channel_width(self, widths_rad):
        self.channel_width = widths_rad
    
    def set_channel_height(self, heights):
        self.channel_heights = heights

    def set_t_wall_tot(self, t_wall_tot):
        self.t_wall_tot = t_wall_tot
    
    def set_blockage_ratio(self, blockage_ratio):
        """blockage_ratio can be scalar or length-N array over x-domain."""
        self.blockage_ratio = blockage_ratio

    def set_x_domain(self, x_domain):
        self.x_domain = x_domain

    def finalize(self):
        self.precompute_thermal_properties()
        self.compute_volume()
        # self.compute_geometry()

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
    • np.inf    → locally straight (|κ| below `eps`)

    Parameters
    ----------
    points : (N, 3) ndarray
        [[x, r, theta], …] ordered along the curve.
    axis   : {'x', 'y', 'z'}, optional
        Which coordinate is the symmetry axis.  Default 'x'.
    eps    : float, optional
        Curvature values with |κ| < eps are treated as zero (straight).

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

    # -------- 3. curvature (κ) and signed radius (R) --------------------- #
    kappa = d2r_dx2 / np.power(1.0 + dr_dx**2, 1.5)   # κ = r″ / (1+r′²)³ᐟ²

    # treat tiny κ as straight line → infinite radius
    with np.errstate(divide="ignore"):
        R = np.where(np.abs(kappa) < eps, np.inf, 1.0 / kappa)

    return R

class CoolingCircuitGroup:
    """Collection of :class:`CoolingCircuit` objects forming the full cooling system.

    Used by :class:`ThrustChamber` to coordinate multi-segment cooling
    layouts, manage overlap, and query active channels along the contour.

    Parameters
    ----------
    circuit_list : list[CoolingCircuit]
        List of circuit instances.
    configuration : str, optional
        Optional label or configuration identifier.

    Attributes
    ----------
    circuits : list[CoolingCircuit]
        All active cooling circuits.

    See Also
    --------
    CoolingCircuit
    ThrustChamber
    """
    def __init__(self, circuit_list, configuration=None):
        self.circuits = circuit_list
        # TODO: implement checks for validity of circuits 
        # TODO: implement helical cooling channels
        # TODO: make function to return number of cooling channels at position x
    
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
        for circuit in self.circuits:
            # Assume circuit.x_domain is a numpy array of x values for this circuit.
            x_start = min(circuit.x_domain[0], circuit.x_domain[-1])
            x_end = max(circuit.x_domain[0], circuit.x_domain[-1])
            if x_start <= x <= x_end:
                if occluding_only and not circuit.placement.occludes:
                    continue
                #total_channels += circuit.placement.n_channel_positions
                total_channels += circuit.placement.channel_count()
        return total_channels




class ChannelPlacement(ABC):
    """Abstract base class defining how coolant channels are positioned.

    Subclasses implement different placement strategies (surface, internal)
    that compute the radial coordinate of the channel centerline.

    Parameters
    ----------
    n_channel_positions : int
        Number of circumferential channel locations (“leaves”).
    channel_width : callable or None, optional
        Function returning angular width [rad]. May be `None` if uniform.
    occludes : bool, optional
        Whether this placement blocks part of the wall from hot gas view.

    Attributes
    ----------
    n_channel_positions : int
        Number of circumferential positions.
    channel_width : callable or None
        Angular width function or constant.
    occludes : bool
        Whether the placement occludes hot surface area.

    See Also
    --------
    SurfacePlacement
    InternalPlacement
    """
    def __init__(self, n_channel_positions: int, channel_width=None, occludes: bool = True):
        self.n_channel_positions = n_channel_positions   # << new home
        self.channel_width = channel_width
        self.occludes = occludes

    @abstractmethod
    def compute_centerline_radius(self,
                                  x: float,
                                  contour,
                                  wall_group) -> float:
        """
        Return the radial coordinate of the coolant-channel centerline.

        Parameters
        ----------
        x : float
            Axial position [m].
        contour : Contour
            Hot-gas contour object.
        wall_group : WallGroup
            Wall stack describing total thickness.

        Returns
        -------
        float
            Radius of the channel centerline [m].
        """

    def channel_count(self) -> int:
        return self.n_channel_positions

class SurfacePlacement(ChannelPlacement):
    """Placement model for surface-mounted cooling channels.

    Channels are positioned just outside the main wall stack,
    offset by total thickness and corrected for local contour angle.

    Parameters
    ----------
    n_channel_positions : int
        Number of channels around the circumference.

    Attributes
    ----------
    n_channel_positions : int
        Number of circumferential leaves.
    n_channels_per_leaf : int
        Fixed to 1 for surface channels.
    occludes : bool
        Always True — these channels block hot-side area.

    See Also
    --------
    ChannelPlacement
    InternalPlacement
    """
    def __init__(self, n_channel_positions: int):
        self.n_channels_per_leaf = 1
        super().__init__(n_channel_positions, channel_width = None, occludes=True)

    def compute_centerline_radius(self, x, contour, wall_group):
        """
        Compute the centerline radius for a surface-mounted channel.

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        contour : Contour
            Hot-gas contour of the chamber/nozzle.
        wall_group : WallGroup
            Wall stack through which the channel is offset.

        Returns
        -------
        float
            Channel centerline radius [m].
        """
        r_hot   = contour.r(x)
        alpha       = contour.normal_angle(x)
        t_total = wall_group.total_thickness(x)
        return r_hot + t_total/np.cos(alpha)

class InternalPlacement(ChannelPlacement):
    """Placement model for in-wall or in-chamber heat-exchanger channels.

    Allows multiple stacked channels per angular leaf and optional
    user-defined width laws.

    Parameters
    ----------
    n_channel_positions : int
        Number of circumferential leaves.
    n_channels_per_leaf : int
        Channels stacked radially per leaf.
    channel_width : callable
        Function returning angular spacing between channel rows [rad].
    occludes : bool, optional
        Whether this placement occludes the hot-side wall. Default False.

    Attributes
    ----------
    n_channel_positions : int
        Number of circumferential leaves.
    n_channels_per_leaf : int
        Channels per leaf.
    channel_width : callable
        Angular spacing function.
    occludes : bool
        Occlusion flag.

    See Also
    --------
    SurfacePlacement
    ChannelPlacement
    """

    def __init__(self,
                    n_channel_positions: int,          # leaves
                    n_channels_per_leaf: int,          # radial stack in each leaf
                    *,
                    channel_width,                          # row-to-row θ-spacing
                    occludes: bool=False):
            
        super().__init__(n_channel_positions,
                        channel_width=channel_width,
                        occludes=occludes)
        
        self.n_channels_per_leaf = n_channels_per_leaf

    def compute_centerline_radius(self, x, contour, wall_group):
        return None

    # override: total = leaves ×  channels / leaf
    """def channel_count(self) -> int:
        return self.n_channel_positions# * self.n_channels_per_leaf"""

class ThrustChamber:
    """Full thrust-chamber assembly combining geometry, cooling, and combustion models.

    Acts as the top-level container linking the hot-gas contour, wall stack,
    and cooling circuits into a coherent physical representation.

    Parameters
    ----------
    contour : Contour
        Hot-gas contour defining inner geometry.
    wall_group : WallGroup
        Structural wall stack.
    cooling_circuit_group : CoolingCircuitGroup
        Collection of cooling circuits.
    combustion_transport : object
        Provides combustion-gas properties and flow variables.
    optimal_values : dict, optional
        Optional dictionary of reference operating conditions.
    roughness : float or callable, optional
        Effective roughness height of coolant walls [m].
    K_factor : float, optional
        Curvature-loss coefficient.
    n_nodes : int, optional
        Number of discrete axial samples.
    h_gas_corr, h_cold_corr : float, optional
        Empirical correction factors for gas- and coolant-side correlations.

    Attributes
    ----------
    contour : Contour
        Geometric shape of the chamber/nozzle.
    wall_group : WallGroup
        Walls through which conduction occurs.
    cooling_circuit_group : CoolingCircuitGroup
        All defined cooling circuits.
    combustion_transport : object
        Hot-gas property model.
    n_nodes : int
        Number of discretization points.
    K_factor : float
        Curvature loss coefficient.
    h_gas_corr, h_cold_corr : float
        Correction multipliers.
    _roughness : float or callable
        Underlying roughness definition.

    Notes
    -----
    The `ThrustChamber` automatically initializes derived circuit geometry
    and may call combustion-transport property generation at construction.

    See Also
    --------
    CoolingCircuit
    WallGroup
    Contour
    """

    def __init__(self, contour, wall_group, cooling_circuit_group, combustion_transport, optimal_values=None, roughness=0.015e-3, K_factor=0.3, n_nodes=50, h_gas_corr=1.0, h_cold_corr=1.0):

        self.contour = contour
        self.wall_group = wall_group
        self.cooling_circuit_group = cooling_circuit_group   # e.g. cooling_circuits.circuits -> list of CoolingCircuit
        self.combustion_transport = combustion_transport
        self.n_nodes = n_nodes
        self.optimal_values = optimal_values

        self.h_gas_corr = h_gas_corr 
        self.h_cold_corr = h_cold_corr

        self._roughness = roughness
        self.K_factor = K_factor

        self.build_circuit_x_domain()
        self.build_channel_centerlines()
        self.build_channel_widths()
        self.build_channel_heights()
        self.build_t_wall_tot()

        for circuit in self.cooling_circuit_group.circuits: 
            circuit.finalize()

        #self.combustion_transport.compute_transport(self.contour) 
        #try:
        self.combustion_transport.compute_aerothermodynamics(self.contour)
        #except Exception:
        #    self.combustion_transport.compute_transport(self.contour)
        # TODO: consider wheather this should be "automated" here or that compute_transport should be done by the user. 
        # I can imagine some scenarioes where simulations are bogged down because the transport properties are automatically
        # computed whenever a thrust chamber is initialised. 

    def build_circuit_x_domain(self):
        """
        Build the x-domain for each cooling circuit by converting its fractional span
        into actual x-values. The sign and ordering of the span determine the coolant flow direction.
        This function uses the overall engine x-range from the contour.
        """
        x_min = self.contour.xs[0]
        x_max = self.contour.xs[-1]
        
        for circuit in self.cooling_circuit_group.circuits:
            f_start, f_end = circuit.span[0], circuit.span[1]
            # For non-negative fractions, multiply by x_max; for negatives, multiply by -x_min.
            x_start = f_start * x_max if f_start >= 0 else f_start * (-x_min)
            x_end   = f_end   * x_max if f_end   >= 0 else f_end   * (-x_min)
            
            # Ensure the domain runs from the lower to the higher x-value.
            if x_start > x_end:
                # If the span is reversed, create the linspace accordingly and then reverse it to maintain the intended flow order.
                x_domain = np.linspace(x_end, x_start, self.n_nodes)[::-1]
            else:
                x_domain = np.linspace(x_start, x_end, self.n_nodes)
            
            circuit.set_x_domain(x_domain)
    # the above one can stay in thrust chamber

    def build_channel_centerlines(self, mode="sim"):
        """
        Build centerline splines for each CoolingCircuit.
        For each circuit, use its pre-built x_domain.
        Each circuit is assigned angles in an interleaved fashion.
        """
        # Determine interleaving of channel angles across circuits.
        if mode == "sim":
            circuit_counts = [1 for _ in self.cooling_circuit_group.circuits]
        elif mode == "plot":
            circuit_counts = [c.placement.n_channel_positions for c in self.cooling_circuit_group.circuits]
        else:
            raise ValueError("Mode must be either 'sim' or 'plot'.")
        
        owners = interleaved_indices(circuit_counts) # TODO: do I need interleaved indecies with this config?
        total_channels = sum(circuit_counts)
        all_angles = np.linspace(0, 2*np.pi, total_channels, endpoint=False)
    
        # For each circuit, build its centerlines.
        for circuit_index, circuit in enumerate(self.cooling_circuit_group.circuits):
            xs_circuit = circuit.x_domain  # use pre-built x-domain
            my_indices = np.where(owners == circuit_index)[0]
            my_angles = all_angles[my_indices]
            local_centerlines = []
            
            for theta_ in my_angles:
                single_centerline = []
                """for x_ in xs_circuit:
                    r_ctr = circuit.placement.compute_centerline_radius(
                                x_, self.contour, self.wall_group)
                    single_centerline.append([x_, r_ctr, theta_])

                single_centerline = np.array(single_centerline)
                local_centerlines.append(single_centerline)"""
                if isinstance(circuit.placement, InternalPlacement):
                    # radial stack
                    
                    for j in range(circuit.placement.n_channels_per_leaf):
                        chain = []
                        for x_ in xs_circuit:
                            r_wall = self.contour.r(x_)
                            h      = circuit.placement.channel_width(x_)
                            r_ctr  = r_wall - (j + 0.5)*h      # inward
                            chain.append([x_, r_ctr, theta_])
                        local_centerlines.append(np.asarray(chain))
                else:
                    chain = []
                    for x_ in xs_circuit:
                        r_ctr = circuit.placement.compute_centerline_radius(
                                    x_, self.contour, self.wall_group)
                        chain.append([x_, r_ctr, theta_])
                    local_centerlines.append(np.asarray(chain))
            circuit.set_centerline(local_centerlines)

    def build_channel_widths(self):
        """
        Compute the channel widths (in radians) for each cooling circuit.
        Uses each circuit's pre-built x_domain and the new number_of_channels(x)
        function to determine the total active channels at each x position.
        """
        for circuit in self.cooling_circuit_group.circuits:
            xs_circuit = circuit.x_domain
            p = circuit.placement

            # --- Case 1: user supplied a width function -----------------------
            if p.channel_width is not None:
                widths = np.array([p.channel_width(x_val) for x_val in xs_circuit])

            elif isinstance(p, InternalPlacement):
                # width = height – rotated 90°
                widths = np.array([circuit.channel_height(x_val)
                                   for x_val in xs_circuit])

            # --- Case 2: default uniform distribution -------------------------
            else:
                widths = []
                for x_val in xs_circuit:
                    n_occ = self.cooling_circuit_group.number_of_channels(
                                x_val, occluding_only=True)
                    width = 2 * np.pi / n_occ if n_occ > 0 else 0.0
                    widths.append(width)
                widths = np.array(widths)

            circuit.set_channel_width(widths)

    # build channel widths is tricky, because it should only evenly distribute if the placement class is surface. 

    def build_channel_heights(self):
        """
        Compute the channel heights for each cooling circuit along its pre-built x_domain.
        Evaluate the channel height function at each x in the circuit's domain.
        """
        for circuit in self.cooling_circuit_group.circuits:
            xs_circuit = circuit.x_domain
            heights = np.array([circuit.channel_height(x_val) for x_val in xs_circuit])
            circuit.set_channel_height(heights)

    def build_t_wall_tot(self):
        """
        Build an array of total wall thicknesses along each circuit's x-domain and
        assign it to the corresponding cooling circuit using set_t_wall_tot.
        """
        for circuit in self.cooling_circuit_group.circuits:
            # Get the x-domain for this circuit
            xs = circuit.x_domain
            
            # Build the wall thickness array along the x-domain.
            # Assumes wall_group.total_thickness(x) returns the thickness at x.
            t_wall_tot_array = np.array([self.wall_group.total_thickness(x) for x in xs])
            
            # Set the computed wall thickness array to the circuit.
            circuit.set_t_wall_tot(t_wall_tot_array)


    def roughness(self, x):
        """
        Get the channel roughness, at a position, x.
        """
        if callable(self._roughness):
            return self._roughness(x)
        else:
            return self._roughness


        

# ----------------------------------------------------------------------
# Helper functions for interleaving channel distribution
# ----------------------------------------------------------------------
def interleaved_indices(circuit_counts):
    """
    Given a list of circuit_counts = [n0, n1, ..., nK], produce an array
    'owners' of length sum(circuit_counts), where each index i is assigned
    to exactly one circuit in an interleaved ratio of n0 : n1 : ... : nK.

    Example: circuit_counts = [30, 60].
    Then we have total=90, ratio=1:2.  The owners array might look like
      [0,1,1, 0,1,1, 0,1,1, ...]
    So circuit #0 gets 30 slots, circuit #1 gets 60 slots, interleaved 1:2.
    """
    # Compute gcd
    g = reduce(gcd, circuit_counts)  # e.g. gcd(30, 60) = 30
    # ratio array, e.g. [1, 2]
    ratios = [c // g for c in circuit_counts]
    block_size = sum(ratios)
    total = sum(circuit_counts)
    owners = np.empty(total, dtype=int)

    # Fill 'owners' block by block
    pos = 0
    for i in range(total):
        offset_in_block = i % block_size
        # figure out which circuit this offset belongs to
        # e.g. if ratios = [1,2], then offset < 1 => circuit0,
        #      if offset < 3 => circuit1, etc.
        circuit_id = 0
        rsum = 0
        for c_idx, r_ in enumerate(ratios):
            rsum_next = rsum + r_
            if offset_in_block < rsum_next:
                circuit_id = c_idx
                break
            rsum = rsum_next
        owners[i] = circuit_id

    return owners
