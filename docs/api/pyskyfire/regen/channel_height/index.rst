pyskyfire.regen.channel_height
==============================

.. py:module:: pyskyfire.regen.channel_height




Module Contents
---------------

.. py:function:: make_channel_height_fn(contour, region_fractions, flat_heights, pinch_factors, transition_widths, logistic_k = 10.0)

   Like before, but now:
     - f = -1.0 maps to the chamber inlet,
     - f =  0.0 maps to the throat (min radius),
     - f = +1.0 maps to the nozzle exit.
   Values in-between interpolate linearly within the chamber or nozzle.


