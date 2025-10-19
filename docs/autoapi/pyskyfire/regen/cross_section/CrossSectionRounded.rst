pyskyfire.regen.cross_section.CrossSectionRounded
=================================================

.. py:class:: pyskyfire.regen.cross_section.CrossSectionRounded(n_points: int = 16)

   Bases: :py:obj:`ChannelSection`

   .. autoapi-inheritance-diagram:: pyskyfire.regen.cross_section.CrossSectionRounded
      :parts: 1
      :private-bases:


   
   Rounded cooling-channel geometry with curved sidewalls.
















   ..
       !! processed by numpydoc !!

   .. py:method:: A_coolant(prof: SectionProfiles) -> numpy.ndarray

      
      Compute effective coolant cross-sectional area [m²].
















      ..
          !! processed by numpydoc !!


   .. py:method:: Dh_coolant(prof: SectionProfiles) -> numpy.ndarray

      
      Compute hydraulic diameter [m].
















      ..
          !! processed by numpydoc !!


   .. py:method:: P_coolant(prof: SectionProfiles) -> numpy.ndarray

      
      Compute coolant-wetted perimeter [m].
















      ..
          !! processed by numpydoc !!


   .. py:method:: P_thermal(prof: SectionProfiles) -> numpy.ndarray

      
      Compute thermal-contact perimeter [m].
















      ..
          !! processed by numpydoc !!


   .. py:method:: _beta_alpha(theta: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]

      
      Return the inner (β) and outer (α) complementary angles.
















      ..
          !! processed by numpydoc !!


   .. py:method:: compute_cross_section(prof: SectionProfiles, i: int) -> int

      
      Construct a gmsh OCC wire representing the rounded section.


      :Parameters:

          **prof** : :obj:`SectionProfiles`
              Section profiles along the cooling circuit.

          **i** : :class:`python:int`
              Station index to build.



      :Returns:

          :class:`python:int`
              gmsh OCC wire tag.








      .. rubric:: Notes

      Requires that a gmsh model is active. The geometry is constructed
      from circle arcs and wall segments in a local coordinate frame,
      then transformed into global coordinates using the provided
      orthonormal basis at station ``i``.



      ..
          !! processed by numpydoc !!

