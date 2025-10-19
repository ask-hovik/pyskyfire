pyskyfire.regen.InternalPlacement
=================================

.. py:class:: pyskyfire.regen.InternalPlacement(n_channel_positions: int, n_channels_per_leaf: int, *, channel_width, occludes: bool = False)

   Bases: :py:obj:`ChannelPlacement`

   .. autoapi-inheritance-diagram:: pyskyfire.regen.InternalPlacement
      :parts: 1
      :private-bases:


   
   Placement model for in-wall or in-chamber heat-exchanger channels.

   Allows multiple stacked channels per angular leaf and optional
   user-defined width laws.

   :Parameters:

       **n_channel_positions** : :class:`python:int`
           Number of circumferential leaves.

       **n_channels_per_leaf** : :class:`python:int`
           Channels stacked radially per leaf.

       **channel_width** : :func:`python:callable`
           Function returning angular spacing between channel rows [rad].

       **occludes** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           Whether this placement occludes the hot-side wall. Default False.

   :Attributes:

       **n_channel_positions** : :class:`python:int`
           Number of circumferential leaves.

       **n_channels_per_leaf** : :class:`python:int`
           Channels per leaf.

       **channel_width** : :func:`python:callable`
           Angular spacing function.

       **occludes** : :ref:`bool <python:bltin-boolean-values>`
           Occlusion flag.









   .. seealso::

       
       :obj:`SurfacePlacement`
           ..
       :obj:`ChannelPlacement`
           ..
       



   ..
       !! processed by numpydoc !!

   .. py:method:: compute_centerline_radius(x, contour, wall_group)

      
      Return the radial coordinate of the coolant-channel centerline.


      :Parameters:

          **x** : :class:`python:float`
              Axial position [m].

          **contour** : :obj:`Contour`
              Hot-gas contour object.

          **wall_group** : :obj:`WallGroup`
              Wall stack describing total thickness.



      :Returns:

          :class:`python:float`
              Radius of the channel centerline [m].











      ..
          !! processed by numpydoc !!

