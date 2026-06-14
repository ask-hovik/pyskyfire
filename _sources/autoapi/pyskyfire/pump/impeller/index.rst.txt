pyskyfire.pump.impeller
=======================

.. py:module:: pyskyfire.pump.impeller


Attributes
----------

.. autoapisummary::

   pyskyfire.pump.impeller.r_a_star_norm
   pyskyfire.pump.impeller.r_i_star_norm
   pyskyfire.pump.impeller.z_a_star_norm
   pyskyfire.pump.impeller.z_i_star_norm


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/pump/impeller/Impeller

.. autoapisummary::

   pyskyfire.pump.impeller.Impeller


Functions
---------

.. autoapisummary::

   pyskyfire.pump.impeller.add_theta_to_streamlines
   pyskyfire.pump.impeller.create_inbetween_streamlines
   pyskyfire.pump.impeller.meridional_streamlines
   pyskyfire.pump.impeller.outlet_diameter
   pyskyfire.pump.impeller.outlet_width
   pyskyfire.pump.impeller.pressure_coefficient
   pyskyfire.pump.impeller.specific_speed


Module Contents
---------------

.. py:function:: add_theta_to_streamlines(meridionals, a, b, c, beta_1B, beta_2B)

   
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















   ..
       !! processed by numpydoc !!

.. py:function:: create_inbetween_streamlines(upper_points, lower_points, num_inbetween=2)

   
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















   ..
       !! processed by numpydoc !!

.. py:function:: meridional_streamlines(n_q, b_2, d_n, d_1, d_2)

.. py:function:: outlet_diameter(psi, n, H)

   
   From Gülich table 7.1, 4th edition 
















   ..
       !! processed by numpydoc !!

.. py:function:: outlet_width(n_q)

   
   Gülich, Eq. 7.1 
















   ..
       !! processed by numpydoc !!

.. py:function:: pressure_coefficient(n_q, z)

   
   From Gülich. Analytical expression for psi. Eq. 3.26, 4th edition
















   ..
       !! processed by numpydoc !!

.. py:function:: specific_speed(n, Q, H)

   
   Calculate centrifugal pump's specific speed.

   From Gülich Table D2.1

   :param n (float): rpm [1/min]
   :param Q (float): flow rate [m^3/s]
   :param H (float): head [m]
   :return n_q (float): specific speed















   ..
       !! processed by numpydoc !!

.. py:data:: r_a_star_norm
   :value: [1.0, 0.9335, 0.8692, 0.8072, 0.7475, 0.6901, 0.6351, 0.5825, 0.5325, 0.4849, 0.4401, 0.3971,...


.. py:data:: r_i_star_norm
   :value: [1.0, 0.8068, 0.6969, 0.6195, 0.561, 0.5134, 0.4729, 0.4374, 0.4056, 0.3767, 0.3503, 0.3253,...


.. py:data:: z_a_star_norm
   :value: [1.0, 0.9986, 0.9945, 0.9878, 0.9784, 0.9664, 0.9519, 0.9349, 0.9155, 0.8938, 0.8698, 0.8437,...


.. py:data:: z_i_star_norm
   :value: [1.0, 0.9911, 0.9735, 0.9526, 0.9302, 0.905, 0.8834, 0.8603, 0.8378, 0.8108, 0.7863, 0.7614,...


