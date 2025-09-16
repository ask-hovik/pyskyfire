from __future__ import annotations
from pyskyfire.regen.cross_section import SectionProfiles
import numpy as np
from abc import ABC, abstractmethod
from math import gcd
from functools import reduce

class Contour:
    def __init__(self, xs, rs, name = None):
        """Class for representing the inner contour of a rocket engine, from the beginning of the combustion chamber to the nozzle exit.

        Args:
            xs (list): Array of x-positions, that the 'y' list corresponds to (m). Must be increasing values of x.
            rs (list): Array, containing local engine radius (m).

        Attributes:
            x_t (float): x-position of the throat (m)
            r_t (float): Throat radius (m)
            A_t (float): Throat area (m2)
            r_e (float): Exit radius (m)
            A_e (float): Exit area (m2)
            r_curvature_t (float): Radius of curvature at the throat (m) """
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
        """Get the distance from the centreline to the inner wall of the engine.
        Args:
            x (float): x position (m)
        Returns:
            float: Distance from engine centreline to edge of inner wall (m)"""
        r1 = np.interp(x, self.xs, self.rs)
        #print(f"r at that point: {r1}, input x: {x} ")
        return r1
    
    def dr_dx(self, x):
        """Get the slope of the engine wall, dr/dx.
        Args:
            x (float): Axial position (m).
        Returns:
            float: Rate of change of contour radius with respect to position, dr/dx """
        return np.interp(x, self.xs, self._dr_dx)

    def A(self, x):
        """Get the flow area for the exhaust gas
        Args:
            x (float): x position (m)
        Returns:
            float: Flow area (m2)"""
        Area = np.pi * self.r(x)**2
        #print(f"Area: {Area}")
        return Area
    
    def normal_angle(self, x): #TODO: check this function af
        """
        Return the smaller angle [0..pi/2] between the outward normal to the
        contour and the plane perpendicular to the x-axis (i.e. the 'vertical'
        direction in this 2D cross-section).

        Geometrically:
          - Tangent: T = (1, dr/dx)
          - Outward normal: N = (-dr/dx, 1)
          - Plane perpendicular to x-axis ~ vertical direction, V = (0, 1)
          - cos(theta) = (N dot V) / (|N| * |V|) = 1 / sqrt((dr/dx)^2 + 1).
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
    """
    Minimum-viable toroidal-aerospike contour class.

    * All **single-radius** properties (`r_t`, `r_c`, etc.) refer to the **outer** wall so
      existing bell-nozzle code keeps working.
    * All **areas** are **annular**:  :math:`A = \pi (r_{o}^{2} - r_{i}^{2})`.

    Parameters
    ----------
    xs_outer : array-like
        Axial sample points of the **outer** contour (m).
    rs_outer : array-like
        Outer radius at the given ``xs_outer`` (m).
    xs_inner : array-like
        Axial sample points of the **inner** contour (m).
    rs_inner : array-like
        Inner radius at the given ``xs_inner`` (m).
    name : str, optional
        Human-readable identifier.
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
        """Annular flow area at *x* (m²)."""
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
        """Radial slope ``dr/dx`` at *x* on chosen wall."""
        if which == "outer":
            return np.interp(x, self.xs_outer, self._dr_dx_outer)
        elif which == "inner":
            return np.interp(x, self.xs_inner, self._dr_dx_inner)
        else:
            raise ValueError("which must be 'outer' or 'inner'")

    def normal_angle(self, x, which="outer"):
        """Angle between outward normal and vertical plane (rad)."""
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
    def __init__(self, material, thickness, name=None):
        """Object for representing an engine wall.

        Args:
            name (str): Name of layer, for example "Chrome Coating" or "Main Wall"
            material (Material): Material object to define the material the wall is made of.
            thickness (float or callable): Thickness of the wall (m). Can be a constant float, or a function of position, i.e. t(x).
        """
        self.name = name
        self.material = material
        self._thickness = thickness

        assert type(thickness) is float or type(thickness) is int or callable(thickness), "'thickness' input must be a float, int or callable"

    def thickness(self, x):
        """Get the thickness of the wall at a position x.

        Args:
            x (float): Axial position along the engine (m)

        Returns:
            float: Wall thickness (m)
        """
        if callable(self._thickness):
            return self._thickness(x)

        else:
            return self._thickness

class WallGroup:
    def __init__(self, walls=None):
        """
        Class that manages multiple Wall objects.

        Args:
            walls (list): A list of Wall objects. If None, defaults to an empty list.
        """
        if walls is None:
            walls = []
        self.walls = walls

    def total_thickness(self, x):
        """
        Returns the sum of the thicknesses of all walls at the given position x.

        Args:
            x (float): Axial position (m)

        Returns:
            float: The total thickness of all walls (m)
        """
        return sum(wall.thickness(x) for wall in self.walls)

'''class CoolingCircuit:
    def __init__(self, name, contour, cross_section, span, placement, channel_height, coolant_transport): 
        """ defines a cooling circuit, which is a section of cooling channel
        
        Args: 
            contour (class): shape of combustion chamber hot wall
            cross_section (class): shape of cooling channel cross section
            n_channels_circuit (int): number of channels in circuit
            span (list): span over which the cooling channel works in x
        """
        self.name = name
        self.contour = contour
        self.cross_section = cross_section
        self.placement = placement
        self.channel_height = channel_height
        self.coolant_transport = coolant_transport

        if span[0] > span[1]:
            self.span = [span[1], span[0]]
            self.direction = -1
        else:
            self.span = span
            self.direction = 1

    def precompute_thermal_properties(self):
        centerline = self.centerlines[0]
        x_vals = centerline[:, 0]
        r_vals = centerline[:, 1]
        theta_vals = centerline[:, 2]
        
        # Interpolate the derivative components at the given x using np.interp
        dx_dx_val = self.centerline_deriv_list[0][:, 0]
        dr_dx_val = self.centerline_deriv_list[0][:, 1]
        dtheta_dx_val = self.centerline_deriv_list[0][:, 2]

        prof = self._make_profiles(centerline, local_coords)
        
        # Compute the stretching factor along the channel:
        #ds_dx = np.sqrt(dx_dx_val**2 + dr_dx_val**2 + (r_vals * dtheta_dx_val)**2) 
        ds_dx = np.sqrt(1 + dr_dx_val**2) # TODO: Update this to actually use everything from the centerline. Currently whack.

        # Compute the perimeter touching the hot exhaust gas
        widths = self.channel_width
        heights = self.channel_heights
        ts_wall = self.t_wall_tot

        hot_perimeter = self.cross_section.P_thermal(heights, widths, ts_wall, centerline) 
        dA_dx_thermal_exhaust_vals = hot_perimeter * ds_dx

        # ==== Hot gas thermal area ====
        self.dA_dx_thermal_exhaust_vals = dA_dx_thermal_exhaust_vals

        cold_perimiter = self.cross_section.P_coolant(heights, widths, ts_wall, centerline) 
        dA_dx_thermal_coolant_vals = cold_perimiter * ds_dx

        # ==== Coolant thermal area ====
        self.dA_dx_thermal_coolant_vals = dA_dx_thermal_coolant_vals

        A_coolant_vals = self.cross_section.A_coolant(heights, widths, ts_wall, centerline) 

        # ==== Coolant cross section ====
        self.A_coolant_vals = A_coolant_vals

        # ==== Coolant cross section derivative ====
        
        self.dA_dx_coolant_vals = np.gradient(A_coolant_vals, x_vals)
        

        Dh_coolant_vals = self.cross_section.Dh_coolant(heights, widths, ts_wall, centerline)

        # ==== Coolant hydraulic diameter ====
        self.Dh_coolant_vals = Dh_coolant_vals
        # TODO: reimplement radius of curvature somehow. 

        radii = radius_of_curvature(centerline)
        self.radius_of_curvature_vals = radii

    def compute_volume(self):
        """
        Calculate the total volume of the cooling channel by integrating
        the product of the cross-sectional area and the local stretching factor
        (i.e. the differential arc length) along the x-direction.

        Returns:
            Total volume (float) of the cooling circuit.
        """
        # Use the primary centerline (assuming all channels are identical)
        
        centerline = self.centerlines[0]
        x_vals = centerline[:, 0]
        r_vals = centerline[:, 1]

        # Extract the precomputed derivatives.
        # Note: In set_centerline(), self.centerline_deriv stores [x, dr/dx, dtheta/dx]
        dr_dx = self.centerline_deriv_list[0][:, 1]
        dtheta_dx = self.centerline_deriv_list[0][:, 2]
        
        # Compute the full stretching factor (ds/dx) accounting for curvature:
        ds_dx = np.sqrt(1 + dr_dx**2 + (r_vals * dtheta_dx)**2)
        
        # Integrate the coolant cross-sectional area multiplied by ds/dx.
        # This gives the volume per channel.
        volume_per_channel = np.trapezoid(self.A_coolant_vals * ds_dx, x_vals)
        
        # Multiply by the number of channels in the circuit.
        total_volume = volume_per_channel * self.placement.n_channel_positions
        self.volume = total_volume

    def compute_geometry(self):
        """
        Compute the discretized cross-sectional point clouds along each channel's centerline.
        For each centerline, the cross_section.compute_point_cloud method is used with the current channel
        geometrical parameters (channel heights, widths, wall thicknesses) and the corresponding local coordinate system.
        
        Results are stored in self.point_cloud as a list of arrays, one per centerline.
        """
        all_point_clouds = []  # list to store point clouds for each centerline

        # Assuming self.channel_heights, self.channel_width, and self.t_wall_tot are arrays
        # defined on the circuit's x-domain that apply to all centerlines.
        for i, centerline in enumerate(self.centerlines):
            # If you stored local coordinate systems for each centerline during set_centerline
            # e.g., self.local_coords_list; if not, consider storing them similarly as in set_centerline.
            local_coords = self.local_coords_list[i]

            # Compute the point cloud for this centerline.
            point_cloud = self.cross_section.compute_point_cloud(
                self.channel_heights,  # h values along the x-domain
                self.channel_width,    # channel angular width (theta) along the x-domain
                self.t_wall_tot,       # total wall thickness along the x-domain
                centerline,            # current centerline (Nx3 array)
                local_coords
            )
            all_point_clouds.append(point_cloud)

        # For backward compatibility, you could choose to also set self.point_cloud
        # to be that computed for the first centerline.
        self.point_cloud = all_point_clouds  # now a list of point clouds, one per channel

    
    def dA_dx_thermal_exhaust(self, x):
        """Return precomputed dA/dx thermal exhaust value at axial position x."""
        return np.interp(x, self.x_domain, self.dA_dx_thermal_exhaust_vals)
    
    def dA_dx_thermal_coolant(self, x):
        """Return precomputed dA/dx thermal coolant value at axial position x."""
        return np.interp(x, self.x_domain, self.dA_dx_thermal_coolant_vals)
    
    def A_coolant(self, x):
        """Return precomputed coolant channel cross-sectional area at axial position x."""
        return np.interp(x, self.x_domain, self.A_coolant_vals)
    
    def dA_dx_coolant(self, x):
        """Return precomputed derivative of coolant channel area at axial position x."""
        return np.interp(x, self.x_domain, self.dA_dx_coolant_vals)
    
    def Dh_coolant(self, x):
        """Return precomputed hydraulic diameter of the coolant channel at axial position x."""
        return np.interp(x, self.x_domain, self.Dh_coolant_vals)
    
    def radius_of_curvature(self, x):
        """Return precomputed local radius of curvature at axial position x."""
        return np.interp(x, self.x_domain, self.radius_of_curvature_vals)

    def set_centerline(self, centerline_list):
        """
        Store a list of centerlines (each given in cylindrical coordinates (x, r, theta))
        and precompute several quantities for each:
        - The 3D coordinates (converted using x, y = r*cos(theta), z = r*sin(theta)).
        - The local tangent vectors (by differentiating the 3D points).
        - The local "normal" vectors computed such that each is perpendicular to the tangent 
            and (if extended) points toward the x-axis.
        - The cylindrical derivatives (dr/dx, dtheta/dx).
        
        Args:
            centerline_list (list of numpy arrays): Each element is an Nx3 array 
                representing a single channel's (x, r, theta) in cylindrical coordinates.
        """
        self.centerlines = centerline_list

        # Create lists to hold computed properties for each centerline.
        #self.centerline_3d_list = []
        #self.
        #self.normal_vectors_list = []
        #self.binormal_vectors_list = []
        self.local_coords_list = []
        self.centerline_deriv_list = []

        # Loop over every centerline and compute properties.
        #tangent_vectors_list = []
        for centerline in centerline_list:
            x_vals = centerline[:, 0]
            r_vals = centerline[:, 1]
            theta_vals = centerline[:, 2]

            # Convert cylindrical to Cartesian coordinates.
            points_3d = np.column_stack((x_vals, r_vals * np.cos(theta_vals), r_vals * np.sin(theta_vals)))#np.column_stack((x_vals, r_vals, theta_vals)) #
            #self.centerline_3d_list.append(points_3d)

            # Compute tangent vectors in 3D via finite differences.
            tangent_vectors = np.zeros_like(points_3d)
            for i in range(3):
                tangent_vectors[:, i] = np.gradient(points_3d[:, i], x_vals)
            # Normalize the tangent vectors.
            norms = np.linalg.norm(tangent_vectors, axis=1, keepdims=True)
            tangent_vectors = tangent_vectors / norms
            #tangent_vectors_list.append(tangent_vectors)

            # Compute "normal" vectors for each point.
            normals = np.zeros_like(points_3d)
            for i, (P, t) in enumerate(zip(points_3d, tangent_vectors)):
                x, y, z = P
                t_x, t_y, t_z = t
                if abs(t_x) < 1e-6:
                    candidate = np.array([0, -y, -z])
                else:
                    candidate = np.array([(y * t_y + z * t_z) / t_x, -y, -z])
                candidate_norm = np.linalg.norm(candidate)
                if candidate_norm > 1e-6:
                    candidate = candidate / candidate_norm
                else:
                    candidate = np.array([0.0, 0.0, 0.0])
                normals[i] = candidate
            #self.normal_vectors_list.append(normals)

            # Compute binormal vectors and stack local coordinate systems.
            binormal_vectors = np.cross(tangent_vectors, normals)
            #self.binormal_vectors_list.append(binormal_vectors)
            local_coords = np.stack((tangent_vectors, normals, binormal_vectors), axis=1)
            self.local_coords_list.append(local_coords)

            # Compute cylindrical derivatives: dr/dx and dtheta/dx.
            dr_dx = np.gradient(r_vals, x_vals)
            dtheta_dx = np.gradient(theta_vals, x_vals)
            centerline_deriv = np.column_stack((x_vals, dr_dx, dtheta_dx))
            self.centerline_deriv_list.append(centerline_deriv)

    def set_channel_width(self, widths):
        """
        Set the channel width (in radians) for this cooling circuit.
        
        Args:
            widths (numpy.ndarray): An array of channel widths computed along the axial domain.
        """
        self.channel_width = widths
    
    def set_channel_height(self, heights):
        """
        Set the channel width (in radians) for this cooling circuit.
        
        Args:
            widths (numpy.ndarray): An array of channel widths computed along the axial domain.
        """
        self.channel_heights = heights

    def set_t_wall_tot(self, t_wall_tot):
        self.t_wall_tot = t_wall_tot

    def set_x_domain(self, x_domain):
        self.x_domain = x_domain

    def finalize(self):
        self.precompute_thermal_properties()
        self.compute_volume()
        #self.compute_geometry()'''

class CoolingCircuit:
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
        list_of_wires = []
        centerl = self.centerlines[0]
        for i in range(len(centerl)):
            local_coords = self.local_coords_list[i]
            prof_i = self._prof(centerl, local_coords)
            wire = self.cross_section.compute_cross_section(prof_i, i)
            list_of_wires.append(wire)
        self.wires = list_of_wires

    def compute_geometry(self):
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

        # === Debug plot of frames in 3D ===
        """fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        ax.plot(points_3d[:, 0], points_3d[:, 1], points_3d[:, 2], "k-", label="centerline")

        scale = 0.05 * np.max(r_vals)  # arrow length scaling
        for P, t, n, b in zip(points_3d, tangent_vectors, normals, binormal_vectors):
            ax.quiver(*P, *(t * scale), color="r", linewidth=0.5)
            ax.quiver(*P, *(n * scale), color="g", linewidth=0.5)
            ax.quiver(*P, *(b * scale), color="b", linewidth=0.5)

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("Local coordinate frames along centerline")
        ax.legend()
        ax.set_box_aspect([1,1,1])
        plt.show()"""




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


"""def plot_local_frames(points_3d, tangents, normals, vec_len=0.05):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y, Z = points_3d[:, 0], points_3d[:, 1], points_3d[:, 2]

    # Tangents
    ax.quiver(X, Y, Z,
              tangents[:, 0], tangents[:, 1], tangents[:, 2],
              length=vec_len, normalize=True, linewidth=0.5)

    # Normals
    ax.quiver(X, Y, Z,
              normals[:, 0], normals[:, 1], normals[:, 2],
              length=vec_len, normalize=True, linewidth=0.5)

    # Plot the x-axis for reference
    x_min, x_max = np.min(X), np.max(X)
    ax.plot([x_min, x_max], [0, 0], [0, 0], linestyle='--', linewidth=1)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_box_aspect((x_max - x_min + 1e-9,
                       np.ptp(Y) + 1e-9,
                       np.ptp(Z) + 1e-9))
    plt.tight_layout()
    plt.show()"""


def old_radius_of_curvature(points): # TODO: not sure where this function should live
    """
    Calculate the radius of curvature along a line given by an array of points.
    Points: numpy array of shape (N, 3) representing [x, r, theta]
    Returns:
        radii: numpy array of radius of curvature at each point (length N)
    """
    def circle_radius(p1, p2, p3):
        # Convert cylindrical (x, r, theta) to Cartesian coordinates (assuming symmetry around x-axis)
        def cyl_to_cart(p):
            x, r, theta = p
            return np.array([x, r * np.cos(theta), r * np.sin(theta)])

        a = cyl_to_cart(p1)
        b = cyl_to_cart(p2)
        c = cyl_to_cart(p3)

        # Calculate vectors
        ab = b - a
        bc = c - b
        ca = a - c

        # Triangle side lengths
        A = np.linalg.norm(bc)
        B = np.linalg.norm(ca)
        C = np.linalg.norm(ab)

        # Semi-perimeter
        s = (A + B + C) / 2

        # Area of triangle via Heron's formula
        area = np.sqrt(max(s * (s - A) * (s - B) * (s - C), 1e-12))  # small epsilon to avoid zero division

        # Circumradius formula
        radius = (A * B * C) / (4 * area)
        return radius

    N = len(points)
    radii = np.zeros(N)

    for i in range(N):
        if i == 0:
            p1, p2, p3 = points[i], points[i + 1], points[i + 2]
        elif i == N - 1:
            p1, p2, p3 = points[i - 2], points[i - 1], points[i]
        else:
            p1, p2, p3 = points[i - 1], points[i], points[i + 1]

        radii[i] = circle_radius(p1, p2, p3)

    return radii



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
    def __init__(self, circuit_list, configuration=None):
        self.circuits = circuit_list
        # TODO: implement checks for validity of circuits 
        # TODO: implement helical cooling channels
        # TODO: make function to return number of cooling channels at position x
    
    def number_of_channels(self, x, *, occluding_only=False):
        """
        Return the total number of cooling channels active at a given x position.
        
        Args:
            x (float): The axial x position along the engine.
        
        Returns:
            int: Total number of cooling channels active at the given x position.
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
        Given axial coordinate x, the hot-gas contour and the wall stack,
        return the r-coordinate of the coolant channel centerline.
        """

    def channel_count(self) -> int:
        return self.n_channel_positions

class SurfacePlacement(ChannelPlacement):
    def __init__(self, n_channel_positions: int):
        self.n_channels_per_leaf = 1
        super().__init__(n_channel_positions, channel_width = None, occludes=True)

    def compute_centerline_radius(self, x, contour, wall_group):
        r_hot   = contour.r(x)
        alpha       = contour.normal_angle(x)
        t_total = wall_group.total_thickness(x)
        return r_hot + t_total/np.cos(alpha)

class InternalPlacement(ChannelPlacement):
    """ In-chamber heat-exchanger channels.  May have their own width law. """
       
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
    def __init__(self, contour, wall_group, cooling_circuit_group, combustion_transport, optimal_values=None, roughness=0.015e-3, K_factor=0.3, n_nodes=50, h_gas_corr=1.0, h_cold_corr=1.0):
        """
        Args:
            contour (Contour): The hot-gas contour of the engine
            walls (WallCollection): Collection of walls, must have walls.total_thickness(x)
            cooling_circuits (CircuitMaster): Master container that holds individual CoolingCircuit objects
            channel_height (callable): A function returning the channel height h(x)
            n_nodes (int): Number of axial subdivisions to use for constructing centerlines
        """
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

    """def build_channel_centerlines(self, mode="sim"):
        
        Build centerline splines for each CoolingCircuit.
        For each circuit, use its pre-built x_domain.
        Each circuit is assigned angles in an interleaved fashion.
        
        # Determine interleaving of channel angles across circuits.
        if mode == "sim":
        # Create a list with one channel per circuit.
            circuit_counts = [1 for _ in self.cooling_circuit_group.circuits]
        elif mode == "plot":
            circuit_counts = [c.n_channels_circuit for c in self.cooling_circuit_group.circuits]
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
                for x_ in xs_circuit:
                    r_contour = self.contour.r(x_)
                    alpha = self.contour.normal_angle(x_)
                    t_wall = self.wall_group.total_thickness(x_)
                    r_centerline = r_contour + t_wall / np.cos(alpha)

                    # Each channel point: [x, radius, angle]
                    single_centerline.append([x_, r_centerline, theta_])
                local_centerlines.append(np.array(single_centerline))
            circuit.set_centerline(local_centerlines)"""

    # build channel centerlines obviously needs to move to the new class

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
