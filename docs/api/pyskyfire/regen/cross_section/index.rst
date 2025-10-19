pyskyfire.regen.cross_section
=============================

.. py:module:: pyskyfire.regen.cross_section




Module Contents
---------------

.. py:class:: SectionProfiles

   All geometry profiles along the centerline (length N).


   .. py:attribute:: h
      :type:  numpy.ndarray


   .. py:attribute:: theta
      :type:  numpy.ndarray


   .. py:attribute:: t_wall
      :type:  numpy.ndarray


   .. py:attribute:: centerline
      :type:  numpy.ndarray


   .. py:attribute:: local_coords
      :type:  numpy.ndarray


   .. py:attribute:: blockage_ratio
      :type:  numpy.ndarray


.. py:class:: ChannelSection(n_points = 16)

   Bases: :py:obj:`abc.ABC`


   Common interface all cross-section shapes must implement.


   .. py:attribute:: n_points
      :value: 16



   .. py:method:: A_coolant(prof)
      :abstractmethod:



   .. py:method:: Dh_coolant(prof)
      :abstractmethod:



   .. py:method:: P_thermal(prof)
      :abstractmethod:



   .. py:method:: P_coolant(prof)
      :abstractmethod:



   .. py:method:: compute_cross_section(prof, i)
      :abstractmethod:



.. py:class:: CrossSectionSquared(n_points = 8)

   Bases: :py:obj:`ChannelSection`


   Common interface all cross-section shapes must implement.


   .. py:method:: A_coolant(prof)


   .. py:method:: Dh_coolant(prof)


   .. py:method:: P_thermal(prof)


   .. py:method:: P_coolant(prof)


   .. py:method:: compute_cross_section(prof, i)

      Build a closed OCC wire (int tag) for station i by:
      1) creating a circle EDGE in XY@origin,
      2) applying an affine transform to place it at (P_i, t_i, n_i, b_i),
      3) wrapping the transformed edge into a wire.
      NOTE: gmsh must be initialized and a model added by the caller.



   .. py:attribute:: n_points
      :value: 16



.. py:class:: CrossSectionRounded(n_points = 16)

   Bases: :py:obj:`ChannelSection`


   Common interface all cross-section shapes must implement.


   .. py:method:: A_coolant(prof)


   .. py:method:: Dh_coolant(prof)


   .. py:method:: P_thermal(prof)


   .. py:method:: P_coolant(prof)


   .. py:method:: compute_cross_section(prof, i)

      Build a closed OCC wire (int tag) for station i by:
      1) creating a circle EDGE in XY@origin,
      2) applying an affine transform to place it at (P_i, t_i, n_i, b_i),
      3) wrapping the transformed edge into a wire.
      NOTE: gmsh must be initialized and a model added by the caller.



   .. py:attribute:: n_points
      :value: 16



