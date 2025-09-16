"""
contour.py

Functinos and methods that together create a thrust chamber contour (combustion chamber + nozzle).



This code takes in two json files theta_n.json and theta_e.json. These two files define the angles of the nozzle
near the throat and the exit of the nozzle. The data comes for discrete length fractions, with discrete
data points, so the code interpolates in the length fraction, area ratio, theta surface to yield the answer. 


References:
    - [1] - The Thrust Optimised Parabolic nozzle, AspireSpace, 
      http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
"""

import numpy as np
import warnings
from scipy.optimize import minimize_scalar, root_scalar

import os
import json
import bisect

from pyskyfire.regen.thrust_chamber import Contour


def get_theta_e_n(length_fraction, epsilon_value):
    """
    This function takes in length fraction and area ratio, and two internally defined json files, 
    theta_n.json and theta_e.json. These two files define the angles of the nozzle near the throat 
    and the exit of the nozzle. The code interpolates twice for each angle in the length fraction, 
    area ratio, theta space to yield the angles. 

    Args:
        length_fraction (float): The fractional length of the nozzle compared to a conical nozzle. A number between 0.60 and 1.00),
        epsilon_value (float): Area ratio of the nozzle
    Returns:
        tuple: A tuple (theta_e, theta_n), both in radians.
            theta_e (float): Interpolated exit angle in radians.
            theta_n (float): Interpolated throat angle in radians.

    Returns:
        (theta_e, theta_n) after performing two-stage 1D interpolation:
         1) For each bounding fraction dataset, interpolate across epsilon.
         2) Interpolate those results across the bounding fractions for the final answer.
    """

    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, 'data')
    theta_e_path = os.path.join(data_dir, 'theta_e.json')
    theta_n_path = os.path.join(data_dir, 'theta_n.json')

    # 1) Load JSON data
    with open(theta_e_path, 'r') as f:
        data_e = json.load(f)
    with open(theta_n_path, 'r') as f:
        data_n = json.load(f)

    # 2) The fraction array we have data for
    available_fractions = [0.60, 0.70, 0.80, 0.90, 1.00]

    def fraction_key(prefix, frac):
        # fraction_key('theta_e', 0.6) -> 'theta_e_60_percent'
        frac_int = int(round(frac * 100))
        return f"{prefix}_{frac_int}_percent"

    # 3) Find the two fractions that bracket length_fraction
    #    (similar to bisect_left but we want to guard boundaries)
    idx = bisect.bisect_left(available_fractions, length_fraction)
    if idx == 0:
        frac_lo = available_fractions[0]
        frac_hi = available_fractions[1]
    elif idx >= len(available_fractions):
        frac_lo = available_fractions[-2]
        frac_hi = available_fractions[-1]
    else:
        frac_hi = available_fractions[idx]
        frac_lo = available_fractions[idx - 1]

    # Helper function: for a given fraction key, do 1D interpolation in epsilon
    def interp_in_epsilon(eps, data_dict, prefix, frac):
        key = fraction_key(prefix, frac)
        eps_array = data_dict[key]['epsilon']
        theta_array = data_dict[key]['theta']
        # Use numpy.interp for linear interpolation
        val = np.interp(eps, eps_array, theta_array)
        return val

    # 4) Interpolate in epsilon for each bounding fraction
    theta_e_lo_val = interp_in_epsilon(epsilon_value, data_e, 'theta_e', frac_lo)
    theta_n_lo_val = interp_in_epsilon(epsilon_value, data_n, 'theta_n', frac_lo)
    theta_e_hi_val = interp_in_epsilon(epsilon_value, data_e, 'theta_e', frac_hi)
    theta_n_hi_val = interp_in_epsilon(epsilon_value, data_n, 'theta_n', frac_hi)

    # 5) Interpolate across the two fraction datasets
    #    If frac_lo == frac_hi, we can just return the lower-value result.
    if abs(frac_hi - frac_lo) < 1e-9:
        final_theta_e = theta_e_lo_val
        final_theta_n = theta_n_lo_val
    else:
        # Build a small array [frac_lo, frac_hi] and [value_lo, value_hi].
        # Then use np.interp on length_fraction
        final_theta_e = np.interp(length_fraction,
                                  [frac_lo, frac_hi],
                                  [theta_e_lo_val, theta_e_hi_val])
        final_theta_n = np.interp(length_fraction,
                                  [frac_lo, frac_hi],
                                  [theta_n_lo_val, theta_n_hi_val])

    return np.radians(final_theta_e), np.radians(final_theta_n)


def get_contour_internal(r_c, r_t, area_ratio, L_c, theta_conv, theta_div, nozzle, R_1f, R_2f, R_3f, length_fraction, export_tikz):
    """
    Generate the full nozzle contour coordinates based on geometric parameters.

    This function calculates the x-coordinates and radii (ys) forming the nozzle contour
    for a thrust chamber. It performs several operations:
      - Computes the entrant and exit throat curves using the provided curvature factors.
      - Constructs the chamber contour, with or without a fillet depending on R_2f.
      - Constructs the nozzle contour based on the specified nozzle type ("rao" or "conical").
      - Concatenates all segments and then processes them to ensure the x-values are strictly increasing.

    Args:
        r_c (float): Chamber radius.
        r_t (float): Throat radius.
        area_ratio (float): Nozzle area ratio.
        L_c (float): Chamber length.
        theta_conv (float): Convergence angle (in radians).
        theta_div (float): Divergence angle (in radians).
        nozzle (str): Specifies the nozzle type, either "rao" or "conical".
        R_1f (float): Scaling factor for the throat entrant curvature radius.
        R_2f (float or None): Scaling factor for the chamber fillet radius. Use 0 or None for a hard corner.
        R_3f (float): Scaling factor for the throat exit curvature radius.
        length_fraction (float): A value between 0.60 and 1.00 used for interpolation.

    Returns:
        tuple: A tuple (xs, ys) where:
            xs (numpy.ndarray): Array of x-coordinates for the nozzle contour.
            ys (numpy.ndarray): Array of corresponding radii for the nozzle contour.

    Raises:
        ValueError: If the nozzle type is not 'rao' or 'conical', or if contour processing fails due to 
                    non-monotonic (non-increasing) x-values.
    """
    
    # create containers for the chamber contour
    xs = []
    ys = []

    # Derived variables
    R1 = R_1f * r_t      # throat entrant curvature radius
    R2 = R_2f * r_t      # chamber fillet radius
    R3 = R_3f * r_t      # throat exit curvature radius
    
    # --------------------
    # What kind of nozzle?
    # --------------------
    if nozzle == "rao":
        try:
            theta_e, theta_n = get_theta_e_n(length_fraction, area_ratio)
        except ValueError as e:
            # need to update value errors for the new dataset
            if "The area ratio provided" in str(e): # TODO: Damn this is ugly af. Tired coding
                warnings.warn(f"{str(e)} Will use a 15 degree cone instead.", stacklevel=2)
                nozzle = "conical"
            else:
                raise ValueError("Some error idk i'm tired")
    elif nozzle == "conical": 
        theta_n = theta_div
    # --------------
    # Throat section
    # --------------

    # Converging (entrant) throat contour radius
    O_entr_x = 0                                    # Entrant radius circle center x-coordinate
    O_entr_y = r_t + R1                             # Entrant radius circle center x-coordinate
    x_tl = -R1*np.sin(theta_conv)                   # x-coordinate of end of arc to the left
    y_tl = (r_t + R1) - R1*np.cos(theta_conv)       # y-coordinate of end of arc to the left

    entrant_throat_angles = np.linspace(-theta_conv, 0, 100)                    # Parametrisation of the angles of the arc
    xs_entrant_throat_arc  = O_entr_x + R1*np.sin(entrant_throat_angles)        # creation of the x-coordinaates of the arc in cartesian
    ys_entrant_throat_arc  = O_entr_y - R1*np.cos(entrant_throat_angles)        # creation of the y-coordinaates of the arc in cartesian

    # Diverging (exit) throat contour radius
    O_exit_x = 0                                    # Exit radius circle center x-coordinate
    O_exit_y = r_t + R3                             # Exit radius circle center y-coordinate
    x_tr = O_exit_x + R3*np.sin(theta_n)                    # x-coordinate of end of arc to the right
    y_tr = O_exit_y - R3*np.cos(theta_n)                    # y-coordinate of end of arc to the right

    exit_throat_angles = np.linspace(0, theta_n, 100)                       # Parametrisatio of the angles of the arc
    xs_exit_throat_arc = O_exit_x + R3*np.sin(exit_throat_angles)           # creation of the x-coordinaates of the arc in cartesian
    ys_exit_throat_arc  = O_exit_y - R3*np.cos(exit_throat_angles)          # creation of the y-coordinaates of the arc in cartesian


    # -----------------
    # Construct chamber 
    # -----------------

    if R_2f == 0:
        # Hard corner (no chamber fillet)
        xs_chamber_fillet = np.array([]) # assign empty arrays to the chamber fillet, so that nothing gets added later
        ys_chamber_fillet = np.array([])

        x_c = -y_tl/np.tan(theta_conv) # x-position of the intersection of the converging section and the chamber cylinder
        xs_conv_section = [x_c, x_tl]
        ys_conv_section = [r_c, y_tl]

        xs_cylindrical = np.linspace(-L_c, x_c, 100)
        ys_cylindrical = np.linspace(r_c, r_c, 100)

    else:
        # Chamber fillet
        O_chamber_fillet_x = x_tl - ((r_c - R2) + R2 * np.cos(theta_conv) - y_tl) / np.tan(theta_conv) - R2 * np.sin(theta_conv) # Chamber fillet circle center x-coordinate
        O_chamber_fillet_y = r_c - R2                                                                                            # Chamber fillet circle center y-coordinate

        x_cr = O_chamber_fillet_x + R2*np.sin(theta_conv)   # x-coordinate of the chamber fillet to the right
        y_cr = O_chamber_fillet_y + R2*np.cos(theta_conv)   # y-coordinate of the chamber fillet to the right

        chamber_fillet_angles = np.linspace(0, theta_conv, 100)
        xs_chamber_fillet = O_chamber_fillet_x + R2 * np.sin(chamber_fillet_angles)
        ys_chamber_fillet = O_chamber_fillet_y + R2 * np.cos(chamber_fillet_angles)

        xs_conv_section = [x_cr, x_tl]
        ys_conv_section = [y_cr, y_tl]

        xs_cylindrical = np.linspace(-L_c, O_chamber_fillet_x, 100)
        ys_cylindrical = np.linspace(r_c, r_c, 100)


    # ----------------
    # Construct Nozzle
    # ----------------

    if nozzle == "conical":
        L_n = (r_t*(np.sqrt(area_ratio) - 1) + R3*(1/np.cos(theta_div) - 1))/np.tan(theta_div)
        y_e = r_t*np.sqrt(area_ratio)
        xs_nozzle = [x_tr, x_tr + L_n]
        ys_nozzle = [y_tr, y_e]

    elif nozzle == "rao":
        # constructing curve based on Rao
        Ex = length_fraction * ((area_ratio)**0.5 - 1) * r_t / np.tan(np.pi/12) # Length fraction times the ideal 15 degree (pi/12) nozzle length TODO: Can't remember why this is different than the one above?
        Ey = r_t * (area_ratio**0.5)
        m1 = np.tan(theta_n)
        m2 = np.tan(theta_e)
        C1 = y_tr - m1 * x_tr
        C2 = Ey - m2 * Ex
        Qx = (C2 - C1) / (m1 - m2)
        Qy = (m1 * C2 - m2 * C1) / (m1 - m2)

        xs_nozzle = []
        ys_nozzle = []
        for t in np.linspace(0, 1, 500):
            xs_nozzle.append((1 - t)**2 * x_tr + 2 * (1 - t) * t * Qx + t**2 * Ex)
            ys_nozzle.append((1 - t)**2 * y_tr + 2 * (1 - t) * t * Qy + t**2 * Ey)
    else: 
        raise ValueError("Nozzle must be either 'conical' or 'rao'.")
    
    # concatenate arrays
    if -L_c > O_chamber_fillet_x: 

        xs = np.concatenate([
            np.array(xs_chamber_fillet),
            np.array(xs_conv_section),
            np.array(xs_entrant_throat_arc),
            np.array(xs_exit_throat_arc), 
            np.array(xs_nozzle)
        ])
        ys = np.concatenate([
            np.array(ys_chamber_fillet),
            np.array(ys_conv_section),
            np.array(ys_entrant_throat_arc),
            np.array(ys_exit_throat_arc), 
            np.array(ys_nozzle)
        ])
    else: 

        xs = np.concatenate([
            np.array(xs_cylindrical),
            np.array(xs_chamber_fillet),
            np.array(xs_conv_section),
            np.array(xs_entrant_throat_arc),
            np.array(xs_exit_throat_arc), 
            np.array(xs_nozzle)
        ])
        ys = np.concatenate([
            np.array(ys_cylindrical),
            np.array(ys_chamber_fillet),
            np.array(ys_conv_section),
            np.array(ys_entrant_throat_arc),
            np.array(ys_exit_throat_arc), 
            np.array(ys_nozzle)
        ])
    idx = np.searchsorted(xs, -L_c)

    if 0 < idx < len(xs):
        x_lo, x_hi = xs[idx-1], xs[idx]
        y_lo, y_hi = ys[idx-1], ys[idx]
        y_cut = y_lo + (y_hi - y_lo)*( (-L_c - x_lo)/(x_hi - x_lo) )
        # now drop all points *before* idx and prepend the cutoff point:
        xs = np.concatenate(([-L_c], xs[idx:]))
        ys = np.concatenate(([y_cut],     ys[idx:]))

   
    # insert xs ys processing here
    processed_xs = []
    processed_ys = []

    for i, (x_val, y_val) in enumerate(zip(xs, ys)):
        if i == 0:
            processed_xs.append(x_val)
            processed_ys.append(y_val)
        else:
            if x_val == processed_xs[-1]:
                # Duplicate x-value found; skip this point.
                continue
            elif x_val > processed_xs[-1]:
                processed_xs.append(x_val)
                processed_ys.append(y_val)
            else:
                # x-value decreases – not strictly increasing.
                print(f"The x-values in the engine contour must be strictly increasing: found {x_val} at index {i} after {processed_xs[-1]}.")
                contour = Contour(xs, ys)
                #plot_contour([contour]) #TODO: reroute into the new viz module
                #plt.show()
                raise ValueError("Cannot continue simulation with current chamber contour")



    # Convert lists back to numpy arrays
    xs = np.array(processed_xs)
    ys = np.array(processed_ys)

    
    return xs, ys


"""def plot_contour(xs, ys):
    fig, ax = plt.subplots()
    line_top = ax.plot(xs, ys)
    color_top = line_top[0].get_color()
    ax.plot(xs, -np.array(ys), color=color_top, label="_nolegend_")
    ax.grid(True)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("Radius (m)")
    ax.set_aspect("equal", "box")
    ax.legend()
    plt.tight_layout()
    plt.show()"""

def compute_chamber_volume(xs, rs):
    """
    Compute the chamber volume by revolving the contour around the x-axis.

    This function calculates the volume of the chamber by integrating the square of the radii 
    (representing a circular cross-section) from the left boundary of the contour up to the throat, 
    which is defined as the point with the minimum radius.

    Args:
        xs (array-like): Sorted array of x-coordinates defining the contour (must be in ascending order).
        rs (array-like): Array of radii corresponding to the x-coordinates.

    Returns:
        float: The computed chamber volume, calculated as π times the integral of r² with respect to x.
    """
    
    xs = np.array(xs)
    rs = np.array(rs)

    # Identify the throat as the index with the minimum radius
    throat_idx = np.argmin(rs)
    
    # If the throat is already at the left boundary, there's no volume to compute
    if throat_idx == 0:
        return 0.0
    
    # Extract every point from x=xs[0] up to x=xs[throat_idx].
    # That is from left boundary (index=0) to the throat (index=throat_idx).
    x_sub = xs[0 : throat_idx + 1]
    r_sub = rs[0 : throat_idx + 1]
    
    # Compute the volume:  V = π * ∫ (r^2 dx)
    volume = np.pi * np.trapezoid(r_sub**2, x_sub)
    return volume


def get_contour(
    r_t, 
    area_ratio,
    r_c=None,
    L_c=None,
    V_c=None,
    eps_c=None,
    AR_c=None,
    theta_conv=45,
    theta_div=15,
    nozzle="rao",
    R_1f=1.5,
    R_2f=0.5,
    R_3f=0.382,
    length_fraction=0.8, 
    angle_input="degrees",
    export_tikz=False 
):
    """
    Generate the nozzle contour (xs, ys) using one of four valid input combinations.

    This function computes the nozzle contour for a thrust chamber using one of the following input methods:
      1. **Direct inputs**: Provide both chamber radius (r_c) and chamber length (L_c).
      2. **Volume & eps**: Provide chamber volume (V_c) and epsilon (eps_c). The chamber radius is computed from eps_c.
      3. **Volume & AR**: Provide chamber volume (V_c) and area ratio (AR_c). The chamber radius is computed 
         from the relation r_c = (L_c * AR_c) / 2. # TODO: I need to update the aspect ratio definition to something more sensible
      4. **Volume & chamber radius**: Provide chamber volume (V_c) and chamber radius (r_c); L_c is determined by minimization.

    The input angles (theta_conv and theta_div) are expected in degrees if `angle_input` is "degrees" 
    and are converted to radians internally.

    Args:
        r_t (float): Throat radius.
        area_ratio (float): Nozzle area ratio (epsilon).
        r_c (float, optional): Chamber radius.
        L_c (float, optional): Chamber length.
        V_c (float, optional): Chamber volume.
        eps_c (float, optional): Epsilon value used to compute the chamber radius.
        AR_c (float, optional): Area ratio used with V_c to compute dimensions.
        theta_conv (float, optional): Convergence angle in degrees (default is 45).
        theta_div (float, optional): Divergence angle in degrees (default is 15).
        nozzle (str, optional): Nozzle type; should be either "rao" or "conical" (default is "rao").
        R_1f (float, optional): Scaling factor for throat entrant curvature (default is 1.5).
        R_2f (float or None, optional): Scaling factor for chamber fillet curvature. Defaults to 0 if None.
        R_3f (float, optional): Scaling factor for throat exit curvature (default is 0.382).
        length_fraction (float, optional): A value between 0.60 and 1.00 used for interpolation (default is 0.8).
        angle_input (str, optional): Unit for theta_conv and theta_div ("degrees" or "radians"). Default is "degrees".

    Returns:
        tuple: A tuple (xs, ys) where:
            xs (numpy.ndarray): Array of x-coordinates for the nozzle contour.
            ys (numpy.ndarray): Array of corresponding radii for the nozzle contour.

    Raises:
        ValueError: If the provided input combination is invalid or if minimization fails 
                    for calculating L_c.
    """

    # TODO: Something big here. So the aspect ratio input is the best input for running an optimizer, because it gives
    # a good chamber shape no matter the input parameters. However, I have found it difficult to make a root finding 
    # scheme that will avoid passing invalid suggestions to the solver. I believe I could fix it if I ran with a "truncated"
    # method!

    if R_2f == None:
        R_2f = 0
    
    if angle_input == "degrees":
        theta_conv = np.radians(theta_conv)
        theta_div = np.radians(theta_div)

    direct_input = (r_c is not None and L_c is not None)
    volume_eps   = (V_c is not None and eps_c is not None)
    volume_AR    = (V_c is not None and AR_c is not None)
    volume_r_c   = (V_c is not None and r_c  is not None and L_c is None 
                    and eps_c is None and AR_c is None)
    
    # --------------------------------------------------------------------------
    # 1) Direct input: r_c, L_c
    # --------------------------------------------------------------------------
    if direct_input and (V_c is None and eps_c is None and AR_c is None):
        xs, ys = get_contour_internal(r_c, r_t, area_ratio, L_c, theta_conv, theta_div, nozzle, R_1f, R_2f, R_3f, length_fraction, export_tikz)
        #V = compute_chamber_volume(xs, ys) # TODO: yeah this should be just a property of the contour object. Need to rewrite this whole file
        #print(f"V: {V}")
        #input()
        return xs, ys
    
    # --------------------------------------------------------------------------
    # 2) Volume & eps: V_c, eps_c
    # --------------------------------------------------------------------------
    elif volume_eps and (r_c is None and L_c is None and AR_c is None):
        A_t = np.pi * r_t**2
        A_c = eps_c * A_t
        r_c = np.sqrt(A_c / np.pi)
        
        def volume_error(Lc_guess):
            xs_guess, ys_guess = get_contour_internal(
                r_c, r_t, area_ratio, Lc_guess, theta_conv, theta_div, nozzle, R_1f, R_2f, R_3f, length_fraction, export_tikz)
            vol = compute_chamber_volume(xs_guess, ys_guess)
            return abs(vol - V_c)
        
        sol = minimize_scalar(volume_error, bounds=(0.1, 10*r_t), method='bounded')
        if not sol.success:
            raise ValueError("Minimization failed for L_c with provided V_c and eps_c.")
        L_c = sol.x
        xs, ys = get_contour_internal(r_c, r_t, area_ratio, L_c, theta_conv, theta_div, nozzle, R_1f, R_2f, R_3f, length_fraction, export_tikz)
        return xs, ys
    

    # --------------------------------------------------------------------------
    # 3) Volume & AR via single‐loop on chamber radius r_c (root-finding on AR error)
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
        # This approach leverages the same volume integration as the old compute_chamber_volume(),
        # which computed chamber volume by rotating the 2D contour about the axis:
        #   V = ∫ π y(x)^2 dx
        # Here each segment is treated as a frustum (3-point slice):
        #   V_slice = π*(y0^2 + y0*y1 + y1^2)/3 * dx
        # so that we mirror the original compute_chamber_volume() method for consistency.
    # 3) Volume & AR via single‐loop on chamber radius r_c (root-finding on AR error)
    # --------------------------------------------------------------------------
    if volume_AR and not (r_c or L_c or eps_c):
        # compute radius bounds
        R1 = R_1f * r_t
        R2 = R_2f * r_t
        r_min = r_t + R1*(1 - np.sin(np.pi/2 - theta_conv)) \
                     + R2*(1 - np.sin(np.pi/2 - theta_conv))
        r_min = r_min + r_min*0.001
        r_max = r_min * 5.0

        L_c_trial = r_t*200

        def ar_constraint(r_c_guess): 
            
            xs_guess, ys_guess = get_contour_internal(
                r_c_guess, r_t, area_ratio, L_c_trial,
                theta_conv, theta_div, nozzle,
                R_1f, R_2f, R_3f, length_fraction, export_tikz)

            L_c = compute_cutoff_length(V_goal = V_c, xs_chamber = xs_guess, ys_chamber = ys_guess)

            S = 2* integrate_area(L_c, xs_guess, ys_guess)

            AR = L_c**2/S
            return AR - AR_c
        
        sol = root_scalar(ar_constraint, bracket=[r_min, r_max], method="bisect")
        r_c_opt = sol.root

        xs_opt, ys_opt = get_contour_internal(
            r_c_opt, r_t, area_ratio, L_c_trial,
            theta_conv, theta_div, nozzle,
            R_1f, R_2f, R_3f, length_fraction, export_tikz
            )
        
        L_c_final = compute_cutoff_length(V_goal=V_c,
                                      xs_chamber=xs_opt,
                                      ys_chamber=ys_opt)
        
        return get_contour_internal(
            r_c_opt, r_t, area_ratio, L_c_final,
            theta_conv, theta_div, nozzle,
            R_1f, R_2f, R_3f, length_fraction, export_tikz
        )





    
    # --------------------------------------------------------------------------
    # 4) Volume & chamber radius: V_c, r_c
    # --------------------------------------------------------------------------
    elif volume_r_c:
        def volume_error(Lc_guess):
            xs_guess, ys_guess = get_contour_internal(r_c, r_t, area_ratio, Lc_guess, theta_conv, theta_div, nozzle, R_1f, R_2f, R_3f, length_fraction, export_tikz)
            vol = compute_chamber_volume(xs_guess, ys_guess)
            return abs(vol - V_c)
        
        sol = minimize_scalar(volume_error, bounds=(0.1, 10*r_t), method='bounded')
        if not sol.success:
            raise ValueError("Minimization failed for L_c with provided V_c and r_c.")
        L_c = sol.x
        xs, ys = get_contour_internal(r_c, r_t, area_ratio, L_c, theta_conv, theta_div, nozzle, R_1f, R_2f, R_3f, length_fraction, export_tikz)
        return xs, ys
    
    else:
        raise ValueError(
            "Invalid input combination. Provide exactly one of the following:\n"
            "  (r_c, L_c), (V_c, eps_c), (V_c, AR_c), or (V_c, r_c).")
    

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d

def compute_cutoff_length(V_goal, xs_chamber, ys_chamber):
    """
    Given a target volume V_goal, and arrays xs_chamber, ys_chamber
    (where xs_chamber runs from some negatives up through positives),
    return L_c = |x_cutoff| such that the volume of revolution about
    the x-axis from x=0 out to x=x_cutoff just reaches V_goal.

    Parameters
    ----------
    V_goal : float
        Desired volume (same units as π * ∫ y^2 dx).
    xs_chamber : array_like, shape (N,)
        x-coordinates, must cover the range from negative up to positive.
    ys_chamber : array_like, shape (N,)
        y-values (assumed ≥0) corresponding to xs_chamber.

    Returns
    -------
    L_c : float
        The absolute distance |x_cutoff| from the origin where the
        cumulative volume first reaches V_goal.

    Raises
    ------
    ValueError
        If V_goal is negative, or larger than the total volume available.
    """
    xs = np.asarray(xs_chamber)
    ys = np.asarray(ys_chamber)
    if V_goal < 0:
        raise ValueError("V_goal must be non-negative.")
    
    # 1) Truncate to xs ≤ 0
    mask = xs <= 0
    if not np.any(mask):
        raise ValueError("No non-positive x values in xs_chamber.")
    x_neg = xs[mask]
    y_neg = ys[mask]

    # 2) Sort so that we start at x=0 and march outward to more negative
    #    (in case xs_chamber wasn't monotonic)
    sort_idx = np.argsort(x_neg)[::-1]  # descending: 0, -0.1, -0.2, ... 
    x_rev = x_neg[sort_idx]  
    y_rev = y_neg[sort_idx]

    # 3) Cumulative integral V(x) = π ∫0→x_rev[i] y^2 dx
    Vs = -np.pi * cumulative_trapezoid(y_rev**2, x_rev, initial=0)

    # 4) Check goal against total available volume
    V_total = Vs[-1]
    if V_goal > V_total:
        raise ValueError(f"V_goal ({V_goal:g}) exceeds total volume ({V_total:g}).")

    # 5) Build an interpolator x(V) and find x_cutoff
    interp_x = interp1d(Vs, x_rev, bounds_error=True)
    x_cutoff = float(interp_x(V_goal))

    # 6) Return the absolute value of that x
    return abs(x_cutoff)


def integrate_area(L_c, xs_chamber, ys_chamber):
    """
    Integrate the area under y(x) from x=0 down to x=-L_c.

    Parameters
    ----------
    L_c : float
        Positive cutoff length. The integration runs from x=0 to x=-L_c.
    xs_chamber : array_like, shape (N,)
        x-coordinates, must cover the range from negative up through positive.
    ys_chamber : array_like, shape (N,)
        y-values corresponding to xs_chamber.

    Returns
    -------
    area : float
        The area ∫_{0}^{-L_c} y(x) dx, returned as a positive number.

    Raises
    ------
    ValueError
        If L_c is negative, or if -L_c lies outside the negative portion of xs_chamber.
    """
    xs = np.asarray(xs_chamber)
    ys = np.asarray(ys_chamber)

    if L_c < 0:
        raise ValueError("L_c must be non-negative.")
    
    # 1) Keep only non-positive x
    mask = xs <= 0
    if not np.any(mask):
        raise ValueError("No non-positive x values in xs_chamber.")
    x_neg = xs[mask]
    y_neg = ys[mask]

    # 2) Sort descending (0, then -small, then more negative)
    sort_idx = np.argsort(x_neg)[::-1]
    x_rev = x_neg[sort_idx]
    y_rev = y_neg[sort_idx]

    # 3) Check that -L_c is within our span
    x_cutoff = -L_c
    x_min = x_rev[-1]
    if x_cutoff < x_min or x_cutoff > 0:
        raise ValueError(f"-L_c = {x_cutoff} lies outside the available x range [{x_min}, 0].")

    # 4) Find where to cut off, and interpolate if needed
    #    find first index i where x_rev[i] <= x_cutoff
    i = np.where(x_rev <= x_cutoff)[0][0]

    # if exact grid point, include it; else interpolate a new point
    if x_rev[i] == x_cutoff:
        x_int = x_rev[:i+1]
        y_int = y_rev[:i+1]
    else:
        x0, x1 = x_rev[i-1], x_rev[i]
        y0, y1 = y_rev[i-1], y_rev[i]
        # linear interpolation for y at x_cutoff
        y_cut = y0 + (y1 - y0) * (x_cutoff - x0) / (x1 - x0)
        x_int = np.concatenate([x_rev[:i], [x_cutoff]])
        y_int = np.concatenate([y_rev[:i], [y_cut]])

    # 5) Trapezoid rule (will be negative because x_int is descending)
    raw = np.trapezoid(y_int, x_int)
    return abs(raw)
