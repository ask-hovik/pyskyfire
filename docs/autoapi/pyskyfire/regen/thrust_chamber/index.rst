pyskyfire.regen.thrust_chamber
==============================

.. py:module:: pyskyfire.regen.thrust_chamber


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/regen/thrust_chamber/ChannelPlacement
   /autoapi/pyskyfire/regen/thrust_chamber/Contour
   /autoapi/pyskyfire/regen/thrust_chamber/ContourToroidalAerospike
   /autoapi/pyskyfire/regen/thrust_chamber/CoolingCircuit
   /autoapi/pyskyfire/regen/thrust_chamber/CoolingCircuitGroup
   /autoapi/pyskyfire/regen/thrust_chamber/InternalPlacement
   /autoapi/pyskyfire/regen/thrust_chamber/SurfacePlacement
   /autoapi/pyskyfire/regen/thrust_chamber/ThrustChamber
   /autoapi/pyskyfire/regen/thrust_chamber/Wall
   /autoapi/pyskyfire/regen/thrust_chamber/WallGroup

.. autoapisummary::

   pyskyfire.regen.thrust_chamber.ChannelPlacement
   pyskyfire.regen.thrust_chamber.Contour
   pyskyfire.regen.thrust_chamber.ContourToroidalAerospike
   pyskyfire.regen.thrust_chamber.CoolingCircuit
   pyskyfire.regen.thrust_chamber.CoolingCircuitGroup
   pyskyfire.regen.thrust_chamber.InternalPlacement
   pyskyfire.regen.thrust_chamber.SurfacePlacement
   pyskyfire.regen.thrust_chamber.ThrustChamber
   pyskyfire.regen.thrust_chamber.Wall
   pyskyfire.regen.thrust_chamber.WallGroup


Functions
---------

.. autoapisummary::

   pyskyfire.regen.thrust_chamber.interleaved_indices
   pyskyfire.regen.thrust_chamber.radius_of_curvature


Module Contents
---------------

.. py:function:: interleaved_indices(circuit_counts)

   
   Given a list of circuit_counts = [n0, n1, ..., nK], produce an array
   'owners' of length sum(circuit_counts), where each index i is assigned
   to exactly one circuit in an interleaved ratio of n0 : n1 : ... : nK.

   Example: circuit_counts = [30, 60].
   Then we have total=90, ratio=1:2.  The owners array might look like
     [0,1,1, 0,1,1, 0,1,1, ...]
   So circuit #0 gets 30 slots, circuit #1 gets 60 slots, interleaved 1:2.















   ..
       !! processed by numpydoc !!

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

