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


   .. py:method:: R_coolant_per_len(prof: SectionProfiles, h_c: numpy.ndarray, k_wall: numpy.ndarray | float) -> numpy.ndarray

      
      Coolant-side thermal resistance per unit channel length [K m / W]. 
      In this function I included a rib on the backside of the cooling channel, 
      basically the same way as standard rib calculations with rib efficiency. 
      Unsure if this is completely appropriate. But it had only a minor effect
      on the RL10 validation case, which use low conductivity stainless. For 
      thin copper walls adding the rib has a large effect. 
















      ..
          !! processed by numpydoc !!


   .. py:method:: _beta_alpha(theta: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]

      
      Return the inner (β) and outer (α) complementary angles.
















      ..
          !! processed by numpydoc !!

