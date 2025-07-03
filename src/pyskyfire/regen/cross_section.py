import numpy as np
import matplotlib.pyplot as plt

class CrossSectionRectangular:
    def __init__(self):
        self.dummy = 0

    def A_coolant(self, h, theta, r, t_wall):
        """ Calculate the cross sectional area of the cooling channel

        Args: 
            h (float or np.array): Height of the cooling channel. A height of 0 represents the shortest the cooling channel can be, while still making geometrical sense
            theta (float or np.array): The angle representing the width of the cooling channel. 
            r (float or np.array): The radius of the thrust chamber at the location of the base of the cooling channel. The base of the cooling channel is the base of the height. 
        Returns: 
            A (float or np.array): Cross sectional area of the cooling channel
        """

        A = 0.5*theta*((r + h)**2 - r**2)
        return A
    
    def Dh_coolant(self, h, theta, r, t_wall): 
        """ Calculate the hydraulic diameter of a cooling channel 
        
        Args: 
            h (float or np.array): Height of the cooling channel. A height of 0 represents the shortest the cooling channel can be, while still making geometrical sense
            theta (float or np.array): The angle representing the width of the cooling channel. 
            r (float or np.array): The radius of the thrust chamber at the location of the base of the cooling channel. The base of the cooling channel is the base of the height. 
        Returns: 
            Dh (float or np.array): Hydraulic diameter of the cooling channel
        """

        A = self.A_coolant(h, theta, r, t_wall)
        P = r*theta + (r + h)*theta + 2*h # perimiter
        Dh = 4*A/P
        return Dh
    
    def P_thermal(self, h, theta, r, t_wall): 
        """ Calculate the perimiter of the cooling channel in contact with the hot combustion gas
        
        Args:
            theta (float or np.array): The angle representing the width of the cooling channel. 
            r (float or np.array): The radius of the thrust chamber at the location of the base of the cooling channel. The base of the cooling channel is the base of the height. 
            t_wall (float or np.array): The thickness of the wall of the cooling channel
        Returns: 
            P (float or np.array): Channel hot gas perimiter"""
        
        P = (r - t_wall)*theta
        return P

    def P_coolant(self, h, theta, r, t_wall): 
        """ Calculate the perimiter of the inside of the cooling channel facing the combustion chamber
        
        Args:
            theta (float or np.array): The angle representing the width of the cooling channel. 
            r (float or np.array): The radius of the thrust chamber at the location of the base of the cooling channel. The base of the cooling channel is the base of the height. 
        Returns: 
            P (float or np.array): Channel cold wall perimiter"""
        
        P = r*theta
        return P

class CrossSectionSquared:
    def __init__(self, n_points=8):
        self.n_points = n_points

    def A_coolant(self, h, theta, t_wall, centerline):

        r = centerline[:, 1]
        r_inner = r + t_wall
        r_outer = r + t_wall + h

        # TODO: update the side width ligic to take in rib width. For now: 
        theta_real = theta/2

        A_sector_outer = theta_real/2 * r_outer**2
        A_sector_inner = theta_real/2 * r_inner**2
        A_sector = A_sector_outer - A_sector_inner
        return A_sector
    
    def Dh_coolant(self, h, theta, t_wall, centerline):
        
        A = self.A_coolant(h, theta, t_wall, centerline)
        
        r = centerline[:, 1]
        r_inner = r + t_wall
        r_outer = r + t_wall + h

        theta_real = theta/2

        P_inner = r_inner*theta_real
        P_outer = r_outer*theta_real
        P = P_inner + P_outer + 2*h

        Dh = 4 * A / P
        return Dh
    
    def P_thermal(self, h, theta, t_wall, centerline):
        
        r = centerline[:, 1]
        theta_real = theta/2
        P_therm = r*theta_real
        return P_therm
    
    def P_coolant(self, h, theta, t_wall, centerline):

        r = centerline[:, 1]
        r_inner = r + t_wall
        theta_real = theta/2
        P_cool = r_inner*theta_real
        return P_cool
    
    def compute_point_cloud(self, h, theta, t_wall, centerline, local_coords):
        """
        Parameters
        ----------
        h, theta, t_wall : 1-D arrays of length N
            Channel height, total included angle, and wall thickness at each
            centre-line station.
        centerline : (N, 3) array
            Tangent-line positions (only the radius component is used here).
        local_coords : (N, 3, 3) array
            Orthonormal frame [tangent, normal, binormal] at every station.

        Returns
        -------
        inner_sections, outer_sections : (N, 2*n_points+4, 3) float32
            Point clouds describing each closed cross-section curve.
            *inner_sections* is the coolant-side surface (wall recessed by
            `t_wall`); *outer_sections* is the hot-gas wall.
        """
        N = centerline.shape[0]
        inner_sections, outer_sections = [], []

        for i in range(N):
            # ---- local scalars ------------------------------------------------
            r       = centerline[i, 1]          # chamber radius
            h_val   = h[i]
            th_val  = theta[i]
            tw      = t_wall[i]
            phi     = th_val * 0.5              # half-span angle

            n, b = local_coords[i, 2], local_coords[i, 1]   # normal, binormal

            # Chordal radii in the N–B plane
            r_in_hs  =  r              #* np.sin(phi)        # hot-side inner
            r_out_hs = (r + 2*tw + h_val)     #* np.sin(phi)        # hot-side outer
            r_in_cs  =  r + tw                    # coolant-side inner
            r_out_cs =  r + h_val + tw                      # coolant-side outer

            center = r * b

            # --- build inner arcs -------------------------------------------------
            a_in  = np.linspace(-th_val/2,  th_val/2,  self.n_points)   # −φ → +φ
            a_out = np.linspace( th_val/2, -th_val/2,  self.n_points)   # +φ → −φ  (reversed!)

            a_in_phi = np.linspace(-phi/2,  phi/2,  self.n_points) 
            a_out_phi = np.linspace( phi/2, -phi/2,  self.n_points)

            arc_in_hs = np.array([center + r_in_hs * (-np.cos(a)*b - np.sin(a)*n)
                                for a in a_in])
            arc_in_cs = np.array([center + r_in_cs * (-np.cos(a)*b - np.sin(a)*n)
                                for a in a_out_phi])

            # --- build outer arcs -------------------------------------------------
            arc_out_hs = np.array([center + r_out_hs * (-np.cos(a)*b - np.sin(a)*n)
                                for a in a_out])          # <-- use the *reversed* angle sweep
            arc_out_cs = np.array([center + r_out_cs * (-np.cos(a)*b - np.sin(a)*n)
                                for a in a_in_phi])           # counterpart for coolant side

            # (keep or drop the “* –1” flip you had before; it doesn’t affect the fix)
            arc_out_hs *= -1;  arc_out_cs *= -1;  arc_in_hs *= -1;  arc_in_cs *= -1

            # --- radial walls -----------------------------------------------------
            conn1_hs = np.linspace(arc_in_hs[-1],  arc_out_hs[0],  2)   # ends now coincide
            conn2_hs = np.linspace(arc_out_hs[-1], arc_in_hs[0],   2)
            conn1_cs = np.linspace(arc_in_cs[-1],  arc_out_cs[0],  2)
            conn2_cs = np.linspace(arc_out_cs[-1], arc_in_cs[0],   2)

            # --- closed curves ----------------------------------------------------
            curve_hs = np.vstack((arc_in_hs, arc_out_hs))
            curve_cs = np.vstack((arc_in_cs, arc_out_cs))

            outer_sections.append(curve_hs)
            inner_sections.append(curve_cs)

        
        #   Plots the last station (or change idx) in its local NB-plane.
        #import matplotlib.pyplot as plt

        idx = -1                 # which axial station to plot
        n_vec, b_vec = local_coords[idx, 1], local_coords[idx, 2]

        x_in  = inner_sections[idx] @ n_vec    # projection onto N
        y_in  = inner_sections[idx] @ b_vec    # projection onto B
        x_out = outer_sections[idx] @ n_vec
        y_out = outer_sections[idx] @ b_vec

        """plt.figure(figsize=(6, 6))
        plt.plot(x_in,  y_in,  'k-o', label='coolant wall')
        plt.plot(x_out, y_out, 'r-o', label='hot-gas wall')
        plt.gca().set_aspect('equal')
        plt.xlabel('local N')
        plt.ylabel('local B')
        plt.title('Squared channel cross-section')
        plt.grid(True)
        plt.legend()
        plt.show()"""
        

        return (np.asarray(inner_sections, dtype=np.float32),
                np.asarray(outer_sections, dtype=np.float32))
    
    """def compute_point_cloud(self, h, theta, t_wall, centerline, local_coords):
        
        N = centerline.shape[0]
        inner_sections = []  # to store the closed curve for each centerline point
        outer_sections = []

        for i in range(N):
            # Extract geometry at the centerline point.
            origin = centerline[i, :]          # 3D position (origin of inner arc)

            r_val = centerline[i, 1]             # the r-coordinate
            theta_val = theta[i]                 # channel angle at this point
            h_val = h[i]   
            t_wall_val = t_wall[i]                    

            theta_real_val = theta_val/2 

            # TODO: make functionality to choose closing jacket thickness, currently just decided as the wall thickness. 
            r_chamber = r_val
            r_chamber_coolant_side = r_val + t_wall_val
            r_coolant_roof = r_val + t_wall_val + h_val
            r_jacket = r_val + 2*t_wall_val + h_val

            normal   = local_coords[i, 1]      # unit vector toward chamber axis
            binormal = local_coords[i, 2]      # completes the right-handed basis

            # --- Compute arc1: inner arc ---
            # Use angles from -β to β.
            angles1 = np.linspace(-theta_real_val/2, theta_real_val/2, self.n_points)
            angles2 = np.linspace(-theta_val/2, theta_val/2, self.n_points)

            arc1_hs = np.array([r1_hs * (-np.cos(a) * normal - np.sin(a) * binormal) for a in angles1])
            arc1_cs = np.array([r1_cs * (-np.cos(a) * normal - np.sin(a) * binormal) for a in angles1])

            # --- Compute arc2: outer arc ---
            # Its center is shifted by h along the normal.
            center2 = h_val * normal
            #center2 = origin + h_val

            # Use angles from α down to -α so that the ordering goes from one end of the arc to the other.
            angles2 = np.linspace(-alpha, alpha, self.n_points)
            arc2_hs = np.array([center2 + r2_hs * (np.cos(a) * normal + np.sin(a) * binormal) for a in angles2])
            arc2_cs = np.array([center2 + r2_cs * (np.cos(a) * normal + np.sin(a) * binormal) for a in angles2])
            #arc2 = np.array([center2 + r2 * (-np.cos(a) - np.sin(a)) for a in angles2])

            # --- Build connector segments ---
            # Connector1: from the endpoint of arc1 (at angle β) to the corresponding point on arc2 (at angle α).
            connector1_hs = np.linspace(arc1_hs[-1], arc2_hs[0], 2)  # two points: start and end
            connector1_cs = np.linspace(arc1_cs[-1], arc2_cs[0], 2)

            # Connector2: from the endpoint of arc2 (at angle -α) to the endpoint of arc1 (at angle -β).
            connector2_hs = np.linspace(arc2_hs[-1], arc1_hs[0], 2)
            connector2_cs = np.linspace(arc2_cs[-1], arc1_cs[0], 2)

            # --- Concatenate segments to form the closed cross-sectional curve ---
            # Order: arc1 (inner arc: from -β to β) + connector1 + arc2 (outer arc: from α to -α) + connector2.
            closed_curve_hs = np.vstack((arc1_hs, connector1_hs, arc2_hs, connector2_hs))
            closed_curve_cs = np.vstack((arc1_cs, connector1_cs, arc2_cs, connector2_cs))

            outer_sections.append(closed_curve_hs)
            inner_sections.append(closed_curve_cs)

            plt.figure(figsize=(6,6))
            # Since the curve lies in the plane spanned by normal ([0,1,0]) and binormal ([0,0,1]),
            # we can plot using the y and z coordinates.
            plt.plot(closed_curve[:, 2], closed_curve[:, 1], 'b-', marker='o')
            plt.title("Closed Cross-Section Curve")
            plt.xlabel("y")
            plt.ylabel("z")
            plt.axis('equal')  # Ensure the axes are equally scaled.
            plt.grid(True)
            plt.show()
        
        # Return as an (N, M, 3) array.
        return np.array(inner_sections, dtype=np.float32), np.array(outer_sections, dtype=np.float32)"""
    
    
class CrossSectionRounded:
    def __init__(self, n_points=16):
        self.n_points = n_points

    def A_coolant(self, h, theta, t_wall, centerline):
        
        r = centerline[:, 1]
        r1 = r * np.sin(theta/2)
        r2 = (r + h) * np.sin(theta/2)
        beta = np.pi - theta/2 - np.pi/2
        alpha = np.pi - beta
        L_side = h * np.cos(theta/2)
        A1 = 0.5 * (r1 - t_wall)**2 * 2 * beta
        A2 = 0.5 * (r2 - t_wall)**2 * 2 * alpha
        S1 = (r1 - t_wall) * L_side
        S2 = (r2 - t_wall) * L_side
        A = A1 + A2 + S1 + S2
        return A
    
    def Dh_coolant(self, h, theta, t_wall, centerline):
        
        r = centerline[:, 1]
        A = self.A_coolant(h, theta, t_wall, centerline)
        r1 = r * np.sin(theta/2)
        r2 = (r + h) * np.sin(theta/2)
        beta = np.pi - theta/2 - np.pi/2
        alpha = np.pi - beta
        arc1 = (r1 - t_wall) * 2 * beta
        L_side = h * np.cos(theta/2)
        arc2 = (r2 - t_wall) * 2 * alpha
        P = arc1 + arc2 + 2 * L_side
        Dh = 4 * A / P
        return Dh
    
    def P_thermal(self, h, theta, t_wall, centerline):
        
        r = centerline[:, 1]
        r1 = r * np.sin(theta/2)
        angle = np.radians(112)
        P = r1 * angle
        return P

    def P_coolant(self, h, theta, t_wall, centerline):
        
        r = centerline[:, 1]        
        r1 = r * np.sin(theta/2)
        angle = np.radians(112)
        P = (r1 - t_wall) * angle
        return P

    def compute_point_cloud(self, h, theta, t_wall, centerline, local_coords):
        
        N = centerline.shape[0]
        inner_sections = []  # to store the closed curve for each centerline point
        outer_sections = []

        for i in range(N):
            # Extract geometry at the centerline point.
            origin = centerline[i, :]          # 3D position (origin of inner arc)
            r_val = centerline[i, 1]             # the r-coordinate
            theta_val = theta[i]                 # channel angle at this point
            h_val = h[i]   
            t_wall_val = t_wall[i]                    # channel height at this point

            # Compute beta and alpha.
            beta = (np.pi - theta_val) / 2.0     # inner arc angular half-span
            alpha = np.pi - beta                 # outer arc angular half-span

            # Compute the arc radii.
            r1_hs = r_val * np.sin(theta_val / 2.0)             # radius for arc1 (inner arc)
            r2_hs = (r_val + h_val) * np.sin(theta_val / 2.0)     # radius for arc2 (outer arc)
            

            r1_cs = r_val * np.sin(theta_val / 2.0) - t_wall_val
            r2_cs = (r_val + h_val) * np.sin(theta_val / 2.0) - t_wall_val     # radius for arc2 (outer arc)

            #print(f"2r = {2*r1}")
            #print(f"theta = {r_val*theta_val}")

            normal   = local_coords[i, 1]      # unit vector toward chamber axis
            binormal = local_coords[i, 2]      # completes the right-handed basis

            # --- Compute arc1: inner arc ---
            # Use angles from -β to β.
            angles1 = np.linspace(-beta, beta, self.n_points)
            arc1_hs = np.array([r1_hs * (-np.cos(a) * normal - np.sin(a) * binormal) for a in angles1])
            arc1_cs = np.array([r1_cs * (-np.cos(a) * normal - np.sin(a) * binormal) for a in angles1])

            # --- Compute arc2: outer arc ---
            # Its center is shifted by h along the normal.
            center2 = h_val * normal
            #center2 = origin + h_val

            # Use angles from α down to -α so that the ordering goes from one end of the arc to the other.
            angles2 = np.linspace(-alpha, alpha, self.n_points)
            arc2_hs = np.array([center2 + r2_hs * (np.cos(a) * normal + np.sin(a) * binormal) for a in angles2])
            arc2_cs = np.array([center2 + r2_cs * (np.cos(a) * normal + np.sin(a) * binormal) for a in angles2])
            #arc2 = np.array([center2 + r2 * (-np.cos(a) - np.sin(a)) for a in angles2])

            # --- Build connector segments ---
            # Connector1: from the endpoint of arc1 (at angle β) to the corresponding point on arc2 (at angle α).
            connector1_hs = np.linspace(arc1_hs[-1], arc2_hs[0], 2)  # two points: start and end
            connector1_cs = np.linspace(arc1_cs[-1], arc2_cs[0], 2)

            # Connector2: from the endpoint of arc2 (at angle -α) to the endpoint of arc1 (at angle -β).
            connector2_hs = np.linspace(arc2_hs[-1], arc1_hs[0], 2)
            connector2_cs = np.linspace(arc2_cs[-1], arc1_cs[0], 2)

            # --- Concatenate segments to form the closed cross-sectional curve ---
            # Order: arc1 (inner arc: from -β to β) + connector1 + arc2 (outer arc: from α to -α) + connector2.
            closed_curve_hs = np.vstack((arc1_hs, connector1_hs, arc2_hs, connector2_hs))
            closed_curve_cs = np.vstack((arc1_cs, connector1_cs, arc2_cs, connector2_cs))

            outer_sections.append(closed_curve_hs)
            inner_sections.append(closed_curve_cs)

            
            """plt.figure(figsize=(6, 6))
            # Project onto local NB-plane: use binormal (index 2) vs normal (index 1)
            plt.plot(closed_curve_cs[:, 2], closed_curve_cs[:, 1],
                     'k-o', label='coolant wall')
            plt.plot(closed_curve_hs[:, 2], closed_curve_hs[:, 1],
                     'r-o', label='hot-gas wall')
            plt.gca().set_aspect('equal')
            plt.xlabel('local B')
            plt.ylabel('local N')
            plt.title(f'Rounded channel cross-section (station {i})')
            plt.grid(True)
            plt.legend()
            plt.show()"""
            
        
        # Return as an (N, M, 3) array.
        return np.array(inner_sections, dtype=np.float32), np.array(outer_sections, dtype=np.float32)
    



"""def compute_point_cloud(self, h, theta, t_wall, centerline, local_coords):
    N = centerline.shape[0]
    inner_sections = []
    outer_sections = []
    for i in range(N):
        # ... [extract origin, r_val, theta_val, h_val, local_dirs] ...

        # inner arc (arc1)
        beta = (np.pi - theta_val) / 2.0
        r1 = r_val * np.sin(theta_val / 2.0)
        angles1 = np.linspace(-beta, beta, self.n_points)
        arc1 = np.array([
            r1 * (-np.cos(a) * normal - np.sin(a) * binormal)
            for a in angles1
        ])

        # outer arc (arc2)
        alpha = np.pi - beta
        r2 = (r_val + h_val) * np.sin(theta_val / 2.0)
        center2 = h_val * normal
        angles2 = np.linspace(-alpha, alpha, self.n_points)
        arc2 = np.array([
            center2 + r2 * (np.cos(a) * normal + np.sin(a) * binormal)
            for a in angles2
        ])

        inner_sections.append(arc1)
        outer_sections.append(arc2)

    # return a tuple of two (N, M, 3) arrays
    return np.array(inner_sections, dtype=np.float32), np.array(outer_sections, dtype=np.float32)"""


class CrossSectionRoundedInternal:
    # for this one, we pretend theta is just in units of length and not angle
    def __init__(self, n_points=16):
        self.n_points = n_points

    def A_coolant(self, h, theta, t_wall, centerline):
        

        A1 = h*theta
        A2 = np.pi*((theta-2*t_wall)/2)**2
        A = A1 + A2
        return A
    
    def Dh_coolant(self, h, theta, t_wall, centerline):
        
        A = self.A_coolant(h, theta, t_wall, centerline)
        P = 2*np.pi*((theta - 2*t_wall)/2) + h*2
        Dh = 4 * A / P
        return Dh
    
    def P_thermal(self, h, theta, t_wall, centerline):
        
        r1 = theta/2
        angle = np.radians(112)
        P = r1 * angle * 2 # *2 because this one has two sides facing combustion
        return P

    def P_coolant(self, h, theta, t_wall, centerline):
        
        r1 = (theta - 2*t_wall)/2
        angle = np.radians(112)
        P = r1 * angle * 2 # *2 because this one has two sides facing combustion
        return P
    
    def compute_point_cloud(self, h, theta, t_wall, centerline, local_coords):
        N          = centerline.shape[0]
        half_pts   = self.n_points          # points per semicircle
        inner, outer = [], []

        for i in range(N):
            w  = theta[i] / 2.0             # half-width  (linear, m)
            hh = h[i] / 2.0                 # half-height (linear, m)
            tw = t_wall[i]

            # local basis
            n  =  local_coords[i, 1]        # radial outward
            b  =  local_coords[i, 2]        # circumferential

            # centres of the two semicircles
            c_plus  =  b * w
            c_minus = -b * w

            # outer and inner radii of the curved ends
            r_o = w
            r_i = w - tw

            # parametric angles for semicircles (flat faces are parallel to n)
            ang = np.linspace(-np.pi/2, np.pi/2, half_pts)
            ang2 = np.linspace(np.pi/2, np.pi*3/2, half_pts)

            # build hot-side (outer) wall
            arc_plus  = [c_plus  + r_o*( np.cos(a)*b + np.sin(a)* n) for a in ang]
            arc_minus = [c_minus + r_o*( np.cos(a)*b - np.sin(a)* n) for a in ang2]
            side_top  = [c_plus  + n*hh,  c_minus + n*hh]
            side_bot  = [c_minus - n*hh,  c_plus  - n*hh]
            outer_sec = np.vstack([arc_plus, side_top, arc_minus[::-1], side_bot])

            # build coolant-side (inner) wall
            arc_plus_i  = [c_plus  + r_i*( np.cos(a)*b + np.sin(a)* n) for a in ang]
            arc_minus_i = [c_minus + r_i*( np.cos(a)*b - np.sin(a)* n) for a in ang2]
            side_top_i  = [c_plus  + (hh-tw)*n, c_minus + (hh-tw)*n]
            side_bot_i  = [c_minus - (hh-tw)*n, c_plus  - (hh-tw)*n]
            inner_sec   = np.vstack([arc_plus_i, side_top_i, arc_minus_i[::-1], side_bot_i])

            outer.append(outer_sec)
            inner.append(inner_sec)

        return np.asarray(inner, dtype=np.float32), np.asarray(outer, dtype=np.float32)
        
    def compute_point_cloud_ols(self, h, theta, t_wall, centerline, local_coords):
        
        N = centerline.shape[0]
        inner_sections = []  # to store the closed curve for each centerline point
        outer_sections = []

        for i in range(N):
            # Extract geometry at the centerline point.
            origin = centerline[i, :]          # 3D position (origin of inner arc)
            r_val = centerline[i, 1]             # the r-coordinate
            theta_val = theta[i]                 # channel angle at this point
            h_val = h[i]   
            t_wall_val = t_wall[i]                    # channel height at this point

            # Compute beta and alpha.
            beta = (np.pi - theta_val) / 2.0     # inner arc angular half-span
            alpha = np.pi - beta                 # outer arc angular half-span

            # Compute the arc radii.
            r1_hs = r_val * np.sin(theta_val / 2.0)             # radius for arc1 (inner arc)
            r2_hs = (r_val + h_val) * np.sin(theta_val / 2.0)     # radius for arc2 (outer arc)
            

            r1_cs = r_val * np.sin(theta_val / 2.0) - t_wall_val
            r2_cs = (r_val + h_val) * np.sin(theta_val / 2.0) - t_wall_val     # radius for arc2 (outer arc)

            #print(f"2r = {2*r1}")
            #print(f"theta = {r_val*theta_val}")

            binormal   = local_coords[i, 1]      # unit vector toward chamber axis
            normal = local_coords[i, 2]      # completes the right-handed basis

            # --- Compute arc1: inner arc ---
            # Use angles from -β to β.
            angles1 = np.linspace(-beta, beta, self.n_points)
            arc1_hs = np.array([r1_hs * (-np.cos(a) * normal - np.sin(a) * binormal) for a in angles1])
            arc1_cs = np.array([r1_cs * (-np.cos(a) * normal - np.sin(a) * binormal) for a in angles1])

            # --- Compute arc2: outer arc ---
            # Its center is shifted by h along the normal.
            center2 = h_val * normal
            #center2 = origin + h_val

            # Use angles from α down to -α so that the ordering goes from one end of the arc to the other.
            angles2 = np.linspace(-alpha, alpha, self.n_points)
            arc2_hs = np.array([center2 + r2_hs * (np.cos(a) * normal + np.sin(a) * binormal) for a in angles2])
            arc2_cs = np.array([center2 + r2_cs * (np.cos(a) * normal + np.sin(a) * binormal) for a in angles2])
            #arc2 = np.array([center2 + r2 * (-np.cos(a) - np.sin(a)) for a in angles2])

            # --- Build connector segments ---
            # Connector1: from the endpoint of arc1 (at angle β) to the corresponding point on arc2 (at angle α).
            connector1_hs = np.linspace(arc1_hs[-1], arc2_hs[0], 2)  # two points: start and end
            connector1_cs = np.linspace(arc1_cs[-1], arc2_cs[0], 2)

            # Connector2: from the endpoint of arc2 (at angle -α) to the endpoint of arc1 (at angle -β).
            connector2_hs = np.linspace(arc2_hs[-1], arc1_hs[0], 2)
            connector2_cs = np.linspace(arc2_cs[-1], arc1_cs[0], 2)

            # --- Concatenate segments to form the closed cross-sectional curve ---
            # Order: arc1 (inner arc: from -β to β) + connector1 + arc2 (outer arc: from α to -α) + connector2.
            closed_curve_hs = np.vstack((arc1_hs, connector1_hs, arc2_hs, connector2_hs))
            closed_curve_cs = np.vstack((arc1_cs, connector1_cs, arc2_cs, connector2_cs))

            outer_sections.append(closed_curve_hs)
            inner_sections.append(closed_curve_cs)

            """plt.figure(figsize=(6,6))
            # Since the curve lies in the plane spanned by normal ([0,1,0]) and binormal ([0,0,1]),
            # we can plot using the y and z coordinates.
            plt.plot(closed_curve[:, 2], closed_curve[:, 1], 'b-', marker='o')
            plt.title("Closed Cross-Section Curve")
            plt.xlabel("y")
            plt.ylabel("z")
            plt.axis('equal')  # Ensure the axes are equally scaled.
            plt.grid(True)
            plt.show()"""
        
        # Return as an (N, M, 3) array.
        return np.array(inner_sections, dtype=np.float32), np.array(outer_sections, dtype=np.float32)
    
    def compute_point_cloud_old(self, h, theta, t_wall, centerline, local_coords):
        
        N = centerline.shape[0]
        cross_sections = []  # to store the closed curve for each centerline point

        for i in range(N):
            # Extract geometry at the centerline point.
            origin = centerline[i, :]          # 3D position (origin of inner arc)
            r_val = centerline[i, 1]             # the r-coordinate
            theta_val = theta[i]                 # channel angle at this point
            h_val = h[i]                       # channel height at this point

            # Compute beta and alpha.
            beta = np.pi/2     # inner arc angular half-span
            alpha = np.pi - beta                 # outer arc angular half-span

            # Compute the arc radii.
            r1 = theta[i]/2             # radius for arc1 (inner arc)
            r2 = theta[i]/2   # radius for arc2 (outer arc)

            #print(f"2r = {2*r1}")
            #print(f"theta = {r_val*theta_val}")



            binormal = np.array([0.0, 1.0, 0.0])
            normal = np.array([0.0, 0.0, 1.0])

            # --- Compute arc1: inner arc ---
            # Use angles from -β to β.
            center1 = -h_val/2 * normal
            angles1 = np.linspace(-beta, beta, self.n_points)
            arc1 = np.array([center1 + r1 * (-np.cos(a) * normal - np.sin(a) * binormal) for a in angles1])
            #arc1 = np.array([[origin[0], origin[1] - r1*np.cos(a), origin[2] + r1*np.sin(a)] for a in angles1])
            #arc1 = np.array([origin + r1 * (np.cos(a) + np.sin(a)) for a in angles1])

            # --- Compute arc2: outer arc ---
            # Its center is shifted by h along the normal.
            center2 = h_val/2 * normal
            #center2 = origin + h_val

            # Use angles from α down to -α so that the ordering goes from one end of the arc to the other.
            angles2 = np.linspace(-alpha, alpha, self.n_points)
            arc2 = np.array([center2 + r2 * (np.cos(a) * normal + np.sin(a) * binormal) for a in angles2])
            #arc2 = np.array([center2 + r2 * (-np.cos(a) - np.sin(a)) for a in angles2])

            # --- Build connector segments ---
            # Connector1: from the endpoint of arc1 (at angle β) to the corresponding point on arc2 (at angle α).
            connector1 = np.linspace(arc1[-1], arc2[0], 2)  # two points: start and end

            # Connector2: from the endpoint of arc2 (at angle -α) to the endpoint of arc1 (at angle -β).
            connector2 = np.linspace(arc2[-1], arc1[0], 2)

            # --- Concatenate segments to form the closed cross-sectional curve ---
            # Order: arc1 (inner arc: from -β to β) + connector1 + arc2 (outer arc: from α to -α) + connector2.
            closed_curve = np.vstack((arc1, connector1, arc2, connector2))

            cross_sections.append(closed_curve)

            """plt.figure(figsize=(6,6))
            # Since the curve lies in the plane spanned by normal ([0,1,0]) and binormal ([0,0,1]),
            # we can plot using the y and z coordinates.
            plt.plot(closed_curve[:, 2], closed_curve[:, 1], 'b-', marker='o')
            plt.title("Closed Cross-Section Curve")
            plt.xlabel("y")
            plt.ylabel("z")
            plt.axis('equal')  # Ensure the axes are equally scaled.
            plt.grid(True)
            plt.show()"""
        
        # Return as an (N, M, 3) array.
        return np.array(cross_sections)