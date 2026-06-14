pyskyfire.regen.thrust_chamber
==============================

.. py:module:: pyskyfire.regen.thrust_chamber


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/regen/thrust_chamber/CoolingCircuit
   /autoapi/pyskyfire/regen/thrust_chamber/FilmCooling
   /autoapi/pyskyfire/regen/thrust_chamber/ThrustChamber
   /autoapi/pyskyfire/regen/thrust_chamber/Wall

.. autoapisummary::

   pyskyfire.regen.thrust_chamber.CoolingCircuit
   pyskyfire.regen.thrust_chamber.FilmCooling
   pyskyfire.regen.thrust_chamber.ThrustChamber
   pyskyfire.regen.thrust_chamber.Wall


Functions
---------

.. autoapisummary::

   pyskyfire.regen.thrust_chamber.radius_of_curvature


Module Contents
---------------

.. py:function:: radius_of_curvature(points: numpy.ndarray, axis: str = 'x', eps: float = 1e-12) -> numpy.ndarray

   
   Signed radius of curvature for a curve expressed in cylindrical coordinates
   [x, r, θ].

   • Positive  → curve bends *away* from the symmetry axis  
   • Negative  → curve bends *toward* the symmetry axis  
   • np.inf    → locally straight (|κ| below `eps`)

   :Parameters:

       **points** : (:obj:`N`, 3) :obj:`ndarray <numpy.ndarray>`
           [[x, r, theta], …] ordered along the curve.

       **axis** : {'x', 'y', 'z'}, :obj:`optional`
           Which coordinate is the symmetry axis.  Default 'x'.

       **eps** : :class:`python:float`, :obj:`optional`
           Curvature values with |κ| < eps are treated as zero (straight).



   :Returns:

       **R** : (N,) :obj:`ndarray <numpy.ndarray>`
           Signed radius of curvature at each sample.











   ..
       !! processed by numpydoc !!

