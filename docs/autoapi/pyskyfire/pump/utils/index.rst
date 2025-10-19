pyskyfire.pump.utils
====================

.. py:module:: pyskyfire.pump.utils


Functions
---------

.. autoapisummary::

   pyskyfire.pump.utils.interpolate_curve
   pyskyfire.pump.utils.tangential_velocity


Module Contents
---------------

.. py:function:: interpolate_curve(points, num_points=200)

   
   Given a list of points (tuples, e.g., [(x1, y1), (x2, y2), ...]) along a curve,
   interpolate the curve to produce a new list of points that are equally spaced 
   in arc-length.

   Parameters:
       points (list of tuples): Original points along the curve.
       num_points (int): Desired number of points in the output.

   Returns:
       new_points (list of tuples): Interpolated points with equal arc-length spacing.















   ..
       !! processed by numpydoc !!

.. py:function:: tangential_velocity(n, d)

   
   Convert rpm and diameter to tangential velocity

   :param n (float): rpm [1/min]
   :param d (float): diameter [m]
   :return u (float): tangential velocity [m/s]















   ..
       !! processed by numpydoc !!

