import numpy as np
from . import constants as const
from . import utils
from . import plot


class Impeller(object):
    """Execute the project of a centrifugal pump."""
    def __init__(self, Q, H, n, material=None):
        """Take input variables and execute the project.

        :param Q (float): flow rate [m^3/s]
        :param H (float) head [m]
        :param n (int): rotational speed [1/min]
        """
        self.Q = Q
        self.H = H
        self.n = n
        self.z = 6

        self.main_dimensions() # calculate the main dimensions and attributes of the pump
        self.geometry = self.compute_geometry() # calculate shroud and blade geometries

        if material is not None: # compute pump weight only in the case where the pump material is supplied
            self.material = material
            self.compute_mass()

        """
        TODO: 
        - Something is wrong with the function that sweeps the angles of the streamlines around along theta. 
        - Basic calculations for the suction eye. This needs to be added. 
        - Sdd routine to determine the number of blades.
        - Acid check the GPT written parts of the code
        - There is still something wrong with the diameter calculation somehow. Or there is some special rule from Gülich I have missed
        - Look through other examples from Gülich and compare
        - The impeller is currently being plotted upside down. Fix this. 
        - Add more traditional plots of the meridional projection and the plan projection. 
        - Add better return contents. Preferably write some stuff to file, so it can be used by other programs

        Future stuff:
        - Add input for type of pump, sewage etc. 
        - Maybe add a few standard examples?
        - Implement thickness of blades. 
        - Add input for resolution to calculate the impeller shape
        - Add function to export mesh for CFD?
        """

    def main_dimensions(self):
        """ Calculate the main dimensions of the impeller """
        
        self.n_q = specific_speed(self.n, self.Q, self.H)
        self.psi = pressure_coefficient(self.n_q, self.z)
        self.d_2 = outlet_diameter(self.psi, self.n, self.H)
        self.u_2 = utils.tangential_velocity(self.n, self.d_2)
        self.b_2 = outlet_width(self.n_q)*self.d_2 #It's supposed to be d_2a. I assume this is for times when it is semi axial
        self.d_1 = 0.171
        self.d_n = 0.07
        """
        results = {}
        for i in ["n_q", "psi", "u_2", "d_2", "b_2"]:
            results[i] = locals()[i]

        return results
        """
    def compute_geometry(self):
        """Compute meridional & and blade geometry"""
        streamlines_2D = meridional_streamlines(n_q=self.n_q, b_2=self.b_2, d_n=self.d_n, d_1=self.d_1, d_2=self.d_2) # this is the meridional section
        streamlines_3D = add_theta_to_streamlines(streamlines_2D, a=1, b=0, c=0, beta_1B=45, beta_2B=50) # actual 3d curvature of blades
        
        return {"meridional": streamlines_2D, "streamlines_3D": streamlines_3D}
    
    def compute_mass(self): 
        mass = self.H*self.material.rho/1000 # TODO: update this placeholder with actual estimate of impeller weight
        self.m_impeller = mass

    def __str__(self):
        """Print out main values in a readable manner"""
        return (
            f"\nImpeller:\n"
            f"  Flow rate, Q = {self.Q:.5f} m^3/s\n"
            f"  Head, H = {self.H:.2f} m\n"
            f"  rpm, n = {self.n} rpm\n"
            f"  Specific Speed, n_q = {self.n_q:.2f} \n"
            f"  Pressure coefficient, psi = {self.psi:.5f} \n"
            f"  Outer diameter, d_2 = {self.d_2:.5f} m\n"
            f"  Peripheral speed, u_2 = {self.u_2:3f} m/s\n"
            f"  Outlet height, b_2 = {self.b_2:.5f} m\n"
            
        )

    def plot_3d(self, a=1, b=0, c=0, beta_1B=45, beta_2B=50, num_blades=6):
        """Plot a 3D view of the impeller geometry."""
        # Do the geometry building & plotting here, rather than in __init__.
        dim = self.results[-1] if self.results else self.main_dimensions()
        meridionals = meridional_streamlines(**dim)
        new_meridionals = add_theta_to_streamlines(
            meridionals, a=a, b=b, c=c,
            beta_1B=beta_1B, beta_2B=beta_2B
        )
        plot.plot_streamlines_3d(new_meridionals, z=num_blades, show=True)


def specific_speed(n, Q, H):
    """Calculate centrifugal pump's specific speed.

    From Gülich Table D2.1

    :param n (float): rpm [1/min]
    :param Q (float): flow rate [m^3/s]
    :param H (float): head [m]
    :return n_q (float): specific speed
    """
    n_q = n*Q**0.5/H**0.75

    return n_q

def pressure_coefficient(n_q, z): 
    """ 
    From Gülich. Analytical expression for psi. Eq. 3.26, 4th edition
    """
    
    n_qRef = 100
    f_T = 1.1 #or 1.1, didn't quite get which is better suited

    psi_opt = 1.21*f_T*np.exp(-0.77*n_q/n_qRef) # 0.9261 (BEP)
    #psi_opt = 1.1*np.exp(-0.0087*n_q/n_qRef) # (sewage pump with 2 blades)
    #psi_0 = 1.31*np.exp(-0.3*n_q/n_qRef) # 1.13740 (diffuser pumps, Q = 0)
    #psi_0 = 1.25*np.exp(-0.3*n_q/n_qRef) # 1.08530 (volute pumps, Q = 0)

    return psi_opt

def outlet_diameter(psi, n, H):
    """ From Gülich table 7.1, 4th edition """
    
    #d_2 = 84.6/n*np.sqrt(H/psi)
    d_2 = 60/(np.pi*n)*np.sqrt(2*const.g*H/psi)

    return d_2

#def inlet_diameter():
    # put the function here

def outlet_width(n_q):
    """ Gülich, Eq. 7.1 
    """
    n_qRef = 100 # 100 or 1. Just had to guess. 
    b_2_star = 0.017 + 0.262*n_q/n_qRef - 0.08*(n_q/n_qRef)**2 + 0.0093*(n_q/n_qRef)**3
    return b_2_star

# some data. TODO: Need to move this later
z_a_star_norm = [1.0000, 0.9986, 0.9945, 0.9878, 0.9784, 0.9664, 0.9519, 0.9349, 0.9155, 0.8938, 0.8698, 0.8437, 0.8156, 0.7855, 0.7537, 0.7201, 0.6850, 0.6484, 0.6106, 0.5716, 0.5317, 0.4910, 0.4496, 0.4079, 0.3658, 0.3237, 0.2818, 0.2402, 0.1992, 0.1723, 0.1458, 0.1199, 0.0944, 0.0820, 0.0697, 0.0576, 0.0456, 0.0335, 0.0224, 0.0111, 0.0000]
r_a_star_norm = [1.0000, 0.9335, 0.8692, 0.8072, 0.7475, 0.6901, 0.6351, 0.5825, 0.5325, 0.4849, 0.4401, 0.3971, 0.3569, 0.3196, 0.2839, 0.2506, 0.2206, 0.1925, 0.1667, 0.1431, 0.1215, 0.1025, 0.0852, 0.0700, 0.0565, 0.0449, 0.0348, 0.0264, 0.0193, 0.0153, 0.0119, 0.0089, 0.0064, 0.0053, 0.0043, 0.0034, 0.0026, 0.0018, 0.0012, 0.0005, 0.0000]
z_i_star_norm = [1.0000, 0.9911, 0.9735, 0.9526, 0.9302, 0.9050, 0.8834, 0.8603, 0.8378, 0.8108, 0.7863, 0.7614, 0.7362, 0.7129, 0.6872, 0.6581, 0.6310, 0.6033, 0.5749, 0.5458, 0.5158, 0.4849, 0.4531, 0.4210, 0.3861, 0.3508, 0.3143, 0.2763, 0.2370, 0.2099, 0.1821, 0.1540, 0.1244, 0.1095, 0.0942, 0.0792, 0.0637, 0.0481, 0.0323, 0.0162, 0.0000]
r_i_star_norm = [1.0000, 0.8068, 0.6969, 0.6195, 0.5610, 0.5134, 0.4729, 0.4374, 0.4056, 0.3767, 0.3503, 0.3253, 0.3021, 0.2803, 0.2597, 0.2404, 0.2214, 0.2036, 0.1865, 0.1701, 0.1543, 0.1392, 0.1246, 0.1106, 0.0972, 0.0843, 0.0720, 0.0602, 0.0498, 0.0418, 0.0349, 0.0283, 0.0220, 0.0190, 0.0159, 0.0131, 0.0107, 0.0076, 0.0050, 0.0024, 0.0000]

def meridional_streamlines(n_q, b_2, d_n, d_1, d_2):    
    #d_n = 0.07 # temporary overrides
    #d_1 = 0.171
    
    n_qRef = 1.0
    resolution = 40

    # Radii
    r_1 = d_1/2
    r_n = d_n/2
    r_2 = d_2/2 

    r_ga_over_r_1 = 2.75*(n_qRef/n_q)**0.16 # Gülich, Table 7.8 (for sewage and xxx)
    r_gi_over_r_ga = 0.83*(n_q/n_qRef)**0.021

    r_ga = r_ga_over_r_1*r_1
    r_gi = r_gi_over_r_ga*r_ga

    # Axial lengths
    z_in = 0.12*(d_1 - d_n)
    z_E = (d_1 - d_n)/2*(n_q/n_qRef)**-0.05*0.75
    z_iSL = (z_E + b_2)*(0.2 + 0.002*(n_q/n_qRef))

    # Assuming eps_DS and eps_TS are 0. Can add expressions here if necessary
    z_ga = z_E
    z_gi = z_E + b_2 - z_iSL

    # streamlines
    z_a_star = []
    for z_star in z_a_star_norm:
        z_a_star.append(z_in + z_ga*z_star)

    z_i_star = []
    for z_star in z_i_star_norm:
        z_i_star.append(z_in + z_iSL + z_gi*z_star)

    r_a_star = []
    for r_star in r_a_star_norm:
        r_a_star.append(r_1 + (r_ga - r_1)*r_star)

    r_i_star = []
    for r_star in r_i_star_norm:
        r_i_star.append(r_n + (r_gi - r_n)*r_star)

    z_a_star.insert(0, z_in + z_ga)
    z_i_star.insert(0, z_in + z_gi + z_iSL)
    r_a_star.insert(0, r_2)
    r_i_star.insert(0, r_2)

    upper_streamline = list(zip(z_a_star, r_a_star)) # upper streamline
    upper_streamline = utils.interpolate_curve(upper_streamline, resolution)
    lower_streamline = list(zip(z_i_star, r_i_star)) # lower streamline
    lower_streamline = utils.interpolate_curve(lower_streamline, resolution)
    streamlines = create_inbetween_streamlines(upper_streamline, lower_streamline, 3)

    """
    # Plot the streamlines.
    plt.plot(*zip(*streamlines[0]))
    plt.plot(*zip(*streamlines[1]))
    plt.plot(*zip(*streamlines[2]))
    plt.plot(*zip(*streamlines[3]))
    plt.plot(*zip(*streamlines[4]))

    # Add labels and a legend.
    plt.xlabel("Axial coordinate (z)")
    plt.ylabel("Radial coordinate (r)")
    plt.legend()
    plt.title("Streamlines")
    plt.grid(True)
    plt.axis("equal")

    # Show the plot.
    plt.show()
    """

    return streamlines

#def inlet_blade_angle(c_1m, tau_1, u_1, c_u1, i_1): 
    # ok. That's a lot of info I don't have yet. 

def create_inbetween_streamlines(upper_points, lower_points, num_inbetween=2):
    """
    Given two lists of evenly spaced points (tuples) representing the upper and lower boundaries,
    generate a set of streamlines that includes the boundaries and a specified number of intermediate
    (in-between) streamlines.

    Parameters:
        upper_points (list of tuples): Points along the upper boundary.
        lower_points (list of tuples): Points along the lower boundary.
        num_inbetween (int): Number of intermediate streamlines to generate (excluding the boundaries).

    Returns:
        streamlines (list): A list of streamlines. Each streamline is a list of (z, r) tuples.
                            The first streamline corresponds to the lower boundary (alpha = 0)
                            and the last to the upper boundary (alpha = 1).
    """
    # Total number of curves includes the two boundaries plus the specified in-between curves.
    total_curves = num_inbetween + 2  
    streamlines = []
    
    # Generate each streamline by linearly interpolating between corresponding points
    # on the lower and upper boundaries.
    for i in range(total_curves):
        # alpha runs from 0 (lower boundary) to 1 (upper boundary)
        alpha = i / (total_curves - 1)
        streamline = [
            ((1 - alpha) * low[0] + alpha * up[0],
             (1 - alpha) * low[1] + alpha * up[1])
            for up, low in zip(upper_points, lower_points)
        ]
        streamlines.append(streamline)
    
    return streamlines

import numpy as np

def add_theta_to_streamlines(meridionals, a, b, c, beta_1B, beta_2B):
    
    """
    Given a list of streamlines (each a list of (z, r) tuples) in the meridional plane,
    compute a theta coordinate for each point based on a prescribed incremental rotation.
    
    For each point (except the first) along a streamline, the incremental angle is given by:
    
        beta_B = beta_B_star * (beta_2B - beta_1B) + beta_1B
        beta_B_star = a*y_star + b*y_star**2 + c*y_star**3
        y_star is a normalized coordinate from 0 to 1 along the streamline.
        
    Parameters:
        meridionals (list): List of streamlines, each a list of (z, r) tuples.
        subdivisions (int or None): If provided, each streamline is assumed to have this number
                                    of points. Otherwise, the function uses the length of each streamline.
        a, b, c (float): Coefficients for the polynomial defining beta_B_star.
        beta_1B (float): The lower bound of the incremental angle.
        beta_2B (float): The upper bound of the incremental angle.
        
    Returns:
        new_meridionals (list): A list of streamlines where each point is a (z, r, theta) tuple.
                                The first point of each streamline is assigned theta = 0.
    """
    new_meridionals = []
    
    for streamline in meridionals:
        # Determine the number of points.
        N = len(streamline)
        # Create a normalized coordinate along the streamline.
        # If subdivisions is provided but the streamline length is different,
        # you might consider interpolating; here we assume they are consistent.
        y_star = np.linspace(0, 1, N)
        
        new_streamline = []
        current_theta = 0  # starting theta
        
        # Loop over the points.
        # We assume the ordering in the input streamline corresponds to increasing y_star.
        for i, (z, r) in enumerate(streamline):
            if i == 0:
                # First point: assign theta = 0.
                new_streamline.append((z, r, 0))
            else:
                # Compute the normalized position for this point.
                # If the number of points in streamline isn't exactly N, use index-based normalization.
                norm_pos = y_star[i] if i < len(y_star) else i / (len(streamline)-1)
                beta_B_star = a * norm_pos + b * (norm_pos ** 2) + c * (norm_pos ** 3)
                delta_theta = beta_B_star * (beta_2B - beta_1B) + beta_1B
                current_theta += delta_theta/20
                new_streamline.append((z, r, current_theta))
                
        new_meridionals.append(new_streamline)
    
    return new_meridionals

