pyskyfire.regen.ChannelSection
==============================

.. py:class:: pyskyfire.regen.ChannelSection(n_points: int = 16)

   Bases: :py:obj:`abc.ABC`

   .. autoapi-inheritance-diagram:: pyskyfire.regen.ChannelSection
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



   .. py:method:: compute_cross_section(prof: SectionProfiles, i: int) -> int
      :abstractmethod:


