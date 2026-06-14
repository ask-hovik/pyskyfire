pyskyfire.regen.cross_section.CrossSectionRoundedInternal
=========================================================

.. py:class:: pyskyfire.regen.cross_section.CrossSectionRoundedInternal(n_points: int = 16)

   Bases: :py:obj:`ChannelSection`

   .. autoapi-inheritance-diagram:: pyskyfire.regen.cross_section.CrossSectionRoundedInternal
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

