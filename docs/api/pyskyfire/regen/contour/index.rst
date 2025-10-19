pyskyfire.regen.contour
=======================

.. py:module:: pyskyfire.regen.contour

.. autoapi-nested-parse::

   contour.py

   Functinos and methods that together create a thrust chamber contour (combustion chamber + nozzle).



   This code takes in two json files theta_n.json and theta_e.json. These two files define the angles of the nozzle
   near the throat and the exit of the nozzle. The data comes for discrete length fractions, with discrete
   data points, so the code interpolates in the length fraction, area ratio, theta surface to yield the answer.


   .. rubric:: References

   - [1] - The Thrust Optimised Parabolic nozzle, AspireSpace,
     http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf





Module Contents
---------------

.. py:function:: get_theta_e_n(length_fraction, epsilon_value)

   This function takes in length fraction and area ratio, and two internally defined json files,
   theta_n.json and theta_e.json. These two files define the angles of the nozzle near the throat
   and the exit of the nozzle. The code interpolates twice for each angle in the length fraction,
   area ratio, theta space to yield the angles.

   :param length_fraction: The fractional length of the nozzle compared to a conical nozzle. A number between 0.60 and 1.00),
   :type length_fraction: float
   :param epsilon_value: Area ratio of the nozzle
   :type epsilon_value: float

   :returns:

             A tuple (theta_e, theta_n), both in radians.
                 theta_e (float): Interpolated exit angle in radians.
                 theta_n (float): Interpolated throat angle in radians.
   :rtype: tuple

   :returns:  1) For each bounding fraction dataset, interpolate across epsilon.
              2) Interpolate those results across the bounding fractions for the final answer.
   :rtype: (theta_e, theta_n) after performing two-stage 1D interpolation


.. py:function:: get_contour_internal(r_c, r_t, area_ratio, L_c, theta_conv, theta_div, nozzle, R_1f, R_2f, R_3f, length_fraction, export_tikz)

   Generate the full nozzle contour coordinates based on geometric parameters.

   This function calculates the x-coordinates and radii (ys) forming the nozzle contour
   for a thrust chamber. It performs several operations:
     - Computes the entrant and exit throat curves using the provided curvature factors.
     - Constructs the chamber contour, with or without a fillet depending on R_2f.
     - Constructs the nozzle contour based on the specified nozzle type ("rao" or "conical").
     - Concatenates all segments and then processes them to ensure the x-values are strictly increasing.

   :param r_c: Chamber radius.
   :type r_c: float
   :param r_t: Throat radius.
   :type r_t: float
   :param area_ratio: Nozzle area ratio.
   :type area_ratio: float
   :param L_c: Chamber length.
   :type L_c: float
   :param theta_conv: Convergence angle (in radians).
   :type theta_conv: float
   :param theta_div: Divergence angle (in radians).
   :type theta_div: float
   :param nozzle: Specifies the nozzle type, either "rao" or "conical".
   :type nozzle: str
   :param R_1f: Scaling factor for the throat entrant curvature radius.
   :type R_1f: float
   :param R_2f: Scaling factor for the chamber fillet radius. Use 0 or None for a hard corner.
   :type R_2f: float or None
   :param R_3f: Scaling factor for the throat exit curvature radius.
   :type R_3f: float
   :param length_fraction: A value between 0.60 and 1.00 used for interpolation.
   :type length_fraction: float

   :returns:

             A tuple (xs, ys) where:
                 xs (numpy.ndarray): Array of x-coordinates for the nozzle contour.
                 ys (numpy.ndarray): Array of corresponding radii for the nozzle contour.
   :rtype: tuple

   :raises ValueError: If the nozzle type is not 'rao' or 'conical', or if contour processing fails due to
       non-monotonic (non-increasing) x-values.


.. py:function:: compute_chamber_volume(xs, rs)

   Compute the chamber volume by revolving the contour around the x-axis.

   This function calculates the volume of the chamber by integrating the square of the radii
   (representing a circular cross-section) from the left boundary of the contour up to the throat,
   which is defined as the point with the minimum radius.

   :param xs: Sorted array of x-coordinates defining the contour (must be in ascending order).
   :type xs: array-like
   :param rs: Array of radii corresponding to the x-coordinates.
   :type rs: array-like

   :returns: The computed chamber volume, calculated as π times the integral of r² with respect to x.
   :rtype: float


.. py:function:: get_contour(r_t, area_ratio, r_c=None, L_c=None, V_c=None, eps_c=None, AR_c=None, theta_conv=45, theta_div=15, nozzle='rao', R_1f=1.5, R_2f=0.5, R_3f=0.382, length_fraction=0.8, angle_input='degrees', export_tikz=False)

   Generate the nozzle contour (xs, ys) using one of four valid input combinations.

   This function computes the nozzle contour for a thrust chamber using one of the following input methods:
     1. **Direct inputs**: Provide both chamber radius (r_c) and chamber length (L_c).
     2. **Volume & eps**: Provide chamber volume (V_c) and epsilon (eps_c). The chamber radius is computed from eps_c.
     3. **Volume & AR**: Provide chamber volume (V_c) and area ratio (AR_c). The chamber radius is computed
        from the relation r_c = (L_c * AR_c) / 2. # TODO: I need to update the aspect ratio definition to something more sensible
     4. **Volume & chamber radius**: Provide chamber volume (V_c) and chamber radius (r_c); L_c is determined by minimization.

   The input angles (theta_conv and theta_div) are expected in degrees if `angle_input` is "degrees"
   and are converted to radians internally.

   :param r_t: Throat radius.
   :type r_t: float
   :param area_ratio: Nozzle area ratio (epsilon).
   :type area_ratio: float
   :param r_c: Chamber radius.
   :type r_c: float, optional
   :param L_c: Chamber length.
   :type L_c: float, optional
   :param V_c: Chamber volume.
   :type V_c: float, optional
   :param eps_c: Epsilon value used to compute the chamber radius.
   :type eps_c: float, optional
   :param AR_c: Area ratio used with V_c to compute dimensions.
   :type AR_c: float, optional
   :param theta_conv: Convergence angle in degrees (default is 45).
   :type theta_conv: float, optional
   :param theta_div: Divergence angle in degrees (default is 15).
   :type theta_div: float, optional
   :param nozzle: Nozzle type; should be either "rao" or "conical" (default is "rao").
   :type nozzle: str, optional
   :param R_1f: Scaling factor for throat entrant curvature (default is 1.5).
   :type R_1f: float, optional
   :param R_2f: Scaling factor for chamber fillet curvature. Defaults to 0 if None.
   :type R_2f: float or None, optional
   :param R_3f: Scaling factor for throat exit curvature (default is 0.382).
   :type R_3f: float, optional
   :param length_fraction: A value between 0.60 and 1.00 used for interpolation (default is 0.8).
   :type length_fraction: float, optional
   :param angle_input: Unit for theta_conv and theta_div ("degrees" or "radians"). Default is "degrees".
   :type angle_input: str, optional

   :returns:

             A tuple (xs, ys) where:
                 xs (numpy.ndarray): Array of x-coordinates for the nozzle contour.
                 ys (numpy.ndarray): Array of corresponding radii for the nozzle contour.
   :rtype: tuple

   :raises ValueError: If the provided input combination is invalid or if minimization fails
       for calculating L_c.


.. py:function:: compute_cutoff_length(V_goal, xs_chamber, ys_chamber)

   Given a target volume V_goal, and arrays xs_chamber, ys_chamber
   (where xs_chamber runs from some negatives up through positives),
   return L_c = |x_cutoff| such that the volume of revolution about
   the x-axis from x=0 out to x=x_cutoff just reaches V_goal.

   :param V_goal: Desired volume (same units as π * ∫ y^2 dx).
   :type V_goal: float
   :param xs_chamber: x-coordinates, must cover the range from negative up to positive.
   :type xs_chamber: array_like, shape (N,)
   :param ys_chamber: y-values (assumed ≥0) corresponding to xs_chamber.
   :type ys_chamber: array_like, shape (N,)

   :returns: **L_c** -- The absolute distance |x_cutoff| from the origin where the
             cumulative volume first reaches V_goal.
   :rtype: float

   :raises ValueError: If V_goal is negative, or larger than the total volume available.


.. py:function:: integrate_area(L_c, xs_chamber, ys_chamber)

   Integrate the area under y(x) from x=0 down to x=-L_c.

   :param L_c: Positive cutoff length. The integration runs from x=0 to x=-L_c.
   :type L_c: float
   :param xs_chamber: x-coordinates, must cover the range from negative up through positive.
   :type xs_chamber: array_like, shape (N,)
   :param ys_chamber: y-values corresponding to xs_chamber.
   :type ys_chamber: array_like, shape (N,)

   :returns: **area** -- The area ∫_{0}^{-L_c} y(x) dx, returned as a positive number.
   :rtype: float

   :raises ValueError: If L_c is negative, or if -L_c lies outside the negative portion of xs_chamber.


