pyskyfire.pump.impeller
=======================

.. py:module:: pyskyfire.pump.impeller








Module Contents
---------------

.. py:class:: Impeller(Q, H, n, material=None)

   Bases: :py:obj:`object`


   Execute the project of a centrifugal pump.

   Take input variables and execute the project.

   :param Q (float): flow rate [m^3/s]
   :param H (float) head [m]
   :param n (int): rotational speed [1/min]


   .. py:attribute:: Q


   .. py:attribute:: H


   .. py:attribute:: n


   .. py:attribute:: z
      :value: 6



   .. py:attribute:: geometry


   .. py:method:: main_dimensions()

      Calculate the main dimensions of the impeller



   .. py:method:: compute_geometry()

      Compute meridional & and blade geometry



   .. py:method:: compute_mass()


   .. py:method:: plot_3d(a=1, b=0, c=0, beta_1B=45, beta_2B=50, num_blades=6)

      Plot a 3D view of the impeller geometry.



.. py:function:: specific_speed(n, Q, H)

   Calculate centrifugal pump's specific speed.

   From Gülich Table D2.1

   :param n (float): rpm [1/min]
   :param Q (float): flow rate [m^3/s]
   :param H (float): head [m]
   :return n_q (float): specific speed


.. py:function:: pressure_coefficient(n_q, z)

   From Gülich. Analytical expression for psi. Eq. 3.26, 4th edition


.. py:function:: outlet_diameter(psi, n, H)

   From Gülich table 7.1, 4th edition


.. py:function:: outlet_width(n_q)

   Gülich, Eq. 7.1



.. py:data:: z_a_star_norm
   :value: [1.0, 0.9986, 0.9945, 0.9878, 0.9784, 0.9664, 0.9519, 0.9349, 0.9155, 0.8938, 0.8698, 0.8437,...


.. py:data:: r_a_star_norm
   :value: [1.0, 0.9335, 0.8692, 0.8072, 0.7475, 0.6901, 0.6351, 0.5825, 0.5325, 0.4849, 0.4401, 0.3971,...


.. py:data:: z_i_star_norm
   :value: [1.0, 0.9911, 0.9735, 0.9526, 0.9302, 0.905, 0.8834, 0.8603, 0.8378, 0.8108, 0.7863, 0.7614,...


.. py:data:: r_i_star_norm
   :value: [1.0, 0.8068, 0.6969, 0.6195, 0.561, 0.5134, 0.4729, 0.4374, 0.4056, 0.3767, 0.3503, 0.3253,...


.. py:function:: meridional_streamlines(n_q, b_2, d_n, d_1, d_2)

.. py:function:: create_inbetween_streamlines(upper_points, lower_points, num_inbetween=2)

   Given two lists of evenly spaced points (tuples) representing the upper and lower boundaries,
   generate a set of streamlines that includes the boundaries and a specified number of intermediate
   (in-between) streamlines.

   :param upper_points: Points along the upper boundary.
   :type upper_points: list of tuples
   :param lower_points: Points along the lower boundary.
   :type lower_points: list of tuples
   :param num_inbetween: Number of intermediate streamlines to generate (excluding the boundaries).
   :type num_inbetween: int

   :returns:

             A list of streamlines. Each streamline is a list of (z, r) tuples.
                                 The first streamline corresponds to the lower boundary (alpha = 0)
                                 and the last to the upper boundary (alpha = 1).
   :rtype: streamlines (list)


.. py:function:: add_theta_to_streamlines(meridionals, a, b, c, beta_1B, beta_2B)

   Given a list of streamlines (each a list of (z, r) tuples) in the meridional plane,
   compute a theta coordinate for each point based on a prescribed incremental rotation.

   For each point (except the first) along a streamline, the incremental angle is given by:

       beta_B = beta_B_star * (beta_2B - beta_1B) + beta_1B
       beta_B_star = a*y_star + b*y_star**2 + c*y_star**3
       y_star is a normalized coordinate from 0 to 1 along the streamline.

   :param meridionals: List of streamlines, each a list of (z, r) tuples.
   :type meridionals: list
   :param subdivisions: If provided, each streamline is assumed to have this number
                        of points. Otherwise, the function uses the length of each streamline.
   :type subdivisions: int or None
   :param a: Coefficients for the polynomial defining beta_B_star.
   :type a: float
   :param b: Coefficients for the polynomial defining beta_B_star.
   :type b: float
   :param c: Coefficients for the polynomial defining beta_B_star.
   :type c: float
   :param beta_1B: The lower bound of the incremental angle.
   :type beta_1B: float
   :param beta_2B: The upper bound of the incremental angle.
   :type beta_2B: float

   :returns:

             A list of streamlines where each point is a (z, r, theta) tuple.
                                     The first point of each streamline is assigned theta = 0.
   :rtype: new_meridionals (list)


