pyskyfire.regen.cross_section.CrossSectionSquared
=================================================

.. py:class:: pyskyfire.regen.cross_section.CrossSectionSquared(blockage_ratio: float, n_points: int = 8)

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


   .. py:method:: R_coolant_per_len(prof: SectionProfiles, h_c: numpy.ndarray, k_wall: numpy.ndarray | float) -> numpy.ndarray

      
      Coolant-side thermal resistance per unit *channel length* [K m / W],
      including rib sidewalls as fins (first-order model).

      Model:
      - Base perimeter = P_coolant(prof)  (what you already count)
      - Fin perimeter  = 2*h             (two side walls per channel)
      - Effective perimeter = P_base + eta_f * P_fin
      - R_s = 1 / (h_c * P_eff)















      ..
          !! processed by numpydoc !!


   .. py:method:: _theta_real(prof: SectionProfiles) -> numpy.ndarray

      
      Apply blockage ratio to effective included angle.
















      ..
          !! processed by numpydoc !!


   .. py:method:: blockage_ratio(x)

