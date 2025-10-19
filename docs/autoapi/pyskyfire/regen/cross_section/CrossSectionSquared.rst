pyskyfire.regen.cross_section.CrossSectionSquared
=================================================

.. py:class:: pyskyfire.regen.cross_section.CrossSectionSquared(n_points: int = 8)

   Bases: :py:obj:`ChannelSection`

   .. autoapi-inheritance-diagram:: pyskyfire.regen.cross_section.CrossSectionSquared
      :parts: 1
      :private-bases:


   
   Simplified rectangular channel section (wedge-sector approximation).
















   ..
       !! processed by numpydoc !!

   .. py:method:: A_coolant(prof: SectionProfiles) -> numpy.ndarray

      
      Compute coolant cross-sectional area [m²].
















      ..
          !! processed by numpydoc !!


   .. py:method:: Dh_coolant(prof: SectionProfiles) -> numpy.ndarray

      
      Compute hydraulic diameter [m].
















      ..
          !! processed by numpydoc !!


   .. py:method:: P_coolant(prof: SectionProfiles) -> numpy.ndarray

      
      Return coolant-wetted perimeter [m].
















      ..
          !! processed by numpydoc !!


   .. py:method:: P_thermal(prof: SectionProfiles) -> numpy.ndarray

      
      Return thermal-contact perimeter [m].
















      ..
          !! processed by numpydoc !!


   .. py:method:: _theta_real(prof: SectionProfiles) -> numpy.ndarray

      
      Apply blockage ratio to effective included angle.
















      ..
          !! processed by numpydoc !!


   .. py:method:: compute_cross_section(prof: SectionProfiles, i: int)

      
      Construct a gmsh OCC wire representing the rectangular section.

      Builds a closed wire via arcs and straight walls positioned at the
      specified centerline station.

      :Parameters:

          **prof** : :obj:`SectionProfiles`
              Full section profile data.

          **i** : :class:`python:int`
              Station index to build.



      :Returns:

          :class:`python:int`
              gmsh OCC wire tag.








      .. rubric:: Notes

      Requires an initialized gmsh model. Only geometric primitives are
      created; meshing is up to the caller.



      ..
          !! processed by numpydoc !!

