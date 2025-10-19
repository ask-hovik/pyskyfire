pyskyfire.regen.thrust_chamber.Contour
======================================

.. py:class:: pyskyfire.regen.thrust_chamber.Contour(xs, rs, name=None)

   
   Inner hot-gas contour of a bell-type rocket engine.

   Represents the geometric wall line from the start of the chamber to
   the nozzle exit, providing local areas, slopes, and radii.

   :Parameters:

       **xs** : :term:`numpy:array_like`
           Axial coordinates of the contour [m], strictly increasing.

       **rs** : :term:`numpy:array_like`
           Corresponding wall radii [m].

       **name** : :class:`python:str`, :obj:`optional`
           Identifier for this contour.

   :Attributes:

       **xs** : :obj:`ndarray <numpy.ndarray>`
           Axial coordinates defining the wall line [m].

       **rs** : :obj:`ndarray <numpy.ndarray>`
           Radial coordinates corresponding to `xs` [m].

       **_dr_dx** : :obj:`ndarray <numpy.ndarray>`
           Precomputed derivative `dr/dx` for fast interpolation.

       **name** : :class:`python:str` or :data:`python:None`
           Optional descriptive name.

       **x_t, r_t, A_t** : :class:`python:float`
           Axial position, radius, and area of the throat.

       **r_e, A_e** : :class:`python:float`
           Radius and area at the exit plane.

       **r_c, A_c** : :class:`python:float`
           Radius and area at the chamber start.

       **eps, eps_c** : :class:`python:float`
           Nozzle and contraction area ratios.









   .. seealso::

       
       :obj:`ContourToroidalAerospike`
           Dual-wall variant for toroidal aerospikes.
       
       
   .. rubric:: Notes

   This class is used by higher-level objects such as
   :class:`~pyskyfire.regen.cooling.CoolingCircuit` and
   :class:`~pyskyfire.regen.thrust.ThrustChamber`.



   ..
       !! processed by numpydoc !!

   .. py:method:: A(x)

      
      Return the local flow area.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].



      :Returns:

          :class:`python:float`
              Flow area at that section [m²].











      ..
          !! processed by numpydoc !!


   .. py:method:: __setattr__(name, value)


   .. py:method:: dr_dx(x)

      
      Return the local slope of the wall, :math:`dr/dx`.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].



      :Returns:

          :class:`python:float`
              Radial slope at position `x`.











      ..
          !! processed by numpydoc !!


   .. py:method:: normal_angle(x)

      
      Return the local wall normal angle with respect to the vertical plane.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].



      :Returns:

          :class:`python:float`
              Angle between the outward normal and the plane normal to the x-axis [rad].











      ..
          !! processed by numpydoc !!


   .. py:method:: r(x)

      
      Return the local radius at axial position `x`.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].



      :Returns:

          :class:`python:float`
              Distance from engine centerline to wall [m].











      ..
          !! processed by numpydoc !!


   .. py:property:: A_c


   .. py:property:: A_e


   .. py:property:: A_t


   .. py:property:: eps


   .. py:property:: eps_c

      
      Get the contraction ratio for the chamber
      Returns: 
          (float): contraction ratio 
















      ..
          !! processed by numpydoc !!


   .. py:property:: r_c


   .. py:property:: r_e


   .. py:property:: r_t


   .. py:property:: x_t

