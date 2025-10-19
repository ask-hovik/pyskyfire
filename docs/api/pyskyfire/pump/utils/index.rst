pyskyfire.pump.utils
====================

.. py:module:: pyskyfire.pump.utils




Module Contents
---------------

.. py:function:: tangential_velocity(n, d)

   Convert rpm and diameter to tangential velocity

   :param n (float): rpm [1/min]
   :param d (float): diameter [m]
   :return u (float): tangential velocity [m/s]


.. py:function:: interpolate_curve(points, num_points=200)

   Given a list of points (tuples, e.g., [(x1, y1), (x2, y2), ...]) along a curve,
   interpolate the curve to produce a new list of points that are equally spaced
   in arc-length.

   :param points: Original points along the curve.
   :type points: list of tuples
   :param num_points: Desired number of points in the output.
   :type num_points: int

   :returns: Interpolated points with equal arc-length spacing.
   :rtype: new_points (list of tuples)


