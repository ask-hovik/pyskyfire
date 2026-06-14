pyskyfire.regen.channel_placement.SurfacePlacement
==================================================

.. py:class:: pyskyfire.regen.channel_placement.SurfacePlacement(n_channel_positions: int, helix_angle=lambda *args: 0.0)

   Bases: :py:obj:`ChannelPlacement`

   .. autoapi-inheritance-diagram:: pyskyfire.regen.channel_placement.SurfacePlacement
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

   .. py:method:: compute_centerline_radius(x, contour)

      
      Return the centerline radius r(x) for this placement.
















      ..
          !! processed by numpydoc !!


   .. py:method:: dtheta_dx(x, contour)

      
      Return dθ/dx [rad/m] for the representative lane centerline.
      Default: 0 (axial/vertical channels).
















      ..
          !! processed by numpydoc !!

