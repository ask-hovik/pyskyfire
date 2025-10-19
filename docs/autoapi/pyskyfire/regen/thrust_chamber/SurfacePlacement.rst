pyskyfire.regen.thrust_chamber.SurfacePlacement
===============================================

.. py:class:: pyskyfire.regen.thrust_chamber.SurfacePlacement(n_channel_positions: int)

   Bases: :py:obj:`ChannelPlacement`

   .. autoapi-inheritance-diagram:: pyskyfire.regen.thrust_chamber.SurfacePlacement
      :parts: 1
      :private-bases:


   
   Placement model for surface-mounted cooling channels.

   Channels are positioned just outside the main wall stack,
   offset by total thickness and corrected for local contour angle.

   :Parameters:

       **n_channel_positions** : :class:`python:int`
           Number of channels around the circumference.

   :Attributes:

       **n_channel_positions** : :class:`python:int`
           Number of circumferential leaves.

       **n_channels_per_leaf** : :class:`python:int`
           Fixed to 1 for surface channels.

       **occludes** : :ref:`bool <python:bltin-boolean-values>`
           Always True — these channels block hot-side area.









   .. seealso::

       
       :obj:`ChannelPlacement`
           ..
       :obj:`InternalPlacement`
           ..
       



   ..
       !! processed by numpydoc !!

   .. py:method:: compute_centerline_radius(x, contour, wall_group)

      
      Compute the centerline radius for a surface-mounted channel.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].

          **contour** : :obj:`Contour`
              Hot-gas contour of the chamber/nozzle.

          **wall_group** : :obj:`WallGroup`
              Wall stack through which the channel is offset.



      :Returns:

          :class:`python:float`
              Channel centerline radius [m].











      ..
          !! processed by numpydoc !!

