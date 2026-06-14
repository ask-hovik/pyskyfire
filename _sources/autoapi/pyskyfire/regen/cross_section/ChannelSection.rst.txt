pyskyfire.regen.cross_section.ChannelSection
============================================

.. py:class:: pyskyfire.regen.cross_section.ChannelSection(n_points: int = 16)

   Bases: :py:obj:`abc.ABC`

   .. autoapi-inheritance-diagram:: pyskyfire.regen.cross_section.ChannelSection
      :parts: 1
      :private-bases:


   
   Abstract base class for cooling-channel cross-section definitions.

   Derived classes must provide analytical expressions for coolant area,
   hydraulic diameter, wetted perimeters, and optionally geometry export.

   :Parameters:

       **n_points** : :class:`python:int`, :obj:`optional`
           Resolution hint for discretized representations. Default is 16.














   ..
       !! processed by numpydoc !!

   .. py:method:: A_coolant(prof: SectionProfiles) -> numpy.ndarray
      :abstractmethod:



   .. py:method:: Dh_coolant(prof: SectionProfiles) -> numpy.ndarray
      :abstractmethod:



   .. py:method:: P_coolant(prof: SectionProfiles) -> numpy.ndarray
      :abstractmethod:



   .. py:method:: P_thermal(prof: SectionProfiles) -> numpy.ndarray
      :abstractmethod:



   .. py:method:: R_coolant_per_len(prof: SectionProfiles, h_c: numpy.ndarray, k_wall: numpy.ndarray | float) -> numpy.ndarray
      :abstractmethod:


      
      Effective coolant-side thermal resistance per unit axial length [K m / W].

      This is where ribs/fin efficiency/spreading conduction can be embedded.

      Default fallback (no rib model):
          R' = 1 / (h_c * A'_cool)
      where A'_cool is the geometric coolant thermal area per unit axial length.















      ..
          !! processed by numpydoc !!


   .. py:method:: R_wall_per_len(prof: SectionProfiles, walls: list, T_rep: numpy.ndarray, A_hot_per_len: numpy.ndarray) -> numpy.ndarray
      :abstractmethod:


      
      Effective wall-stack conduction resistance per unit axial length [K m / W].
















      ..
          !! processed by numpydoc !!

