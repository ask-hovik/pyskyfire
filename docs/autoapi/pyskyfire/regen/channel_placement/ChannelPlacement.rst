pyskyfire.regen.channel_placement.ChannelPlacement
==================================================

.. py:class:: pyskyfire.regen.channel_placement.ChannelPlacement(n_channel_positions: int, channel_width=None, occludes: bool = True)

   Bases: :py:obj:`abc.ABC`

   .. autoapi-inheritance-diagram:: pyskyfire.regen.channel_placement.ChannelPlacement
      :parts: 1
      :private-bases:


   
   Abstract base class defining how coolant channels are positioned.

   Subclasses implement different placement strategies (surface, internal)
   that compute the radial coordinate of the channel centerline.

   :Parameters:

       **n_channel_positions** : :class:`python:int`
           Number of circumferential channel locations (“leaves”).

       **channel_width** : :func:`python:callable` or :data:`python:None`, :obj:`optional`
           Function returning angular width [rad]. May be `None` if uniform.

       **occludes** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           Whether this placement blocks part of the wall from hot gas view.

   :Attributes:

       **n_channel_positions** : :class:`python:int`
           Number of circumferential positions.

       **channel_width** : :func:`python:callable` or :data:`python:None`
           Angular width function or constant.

       **occludes** : :ref:`bool <python:bltin-boolean-values>`
           Whether the placement occludes hot surface area.









   .. seealso::

       
       :obj:`SurfacePlacement`
           ..
       :obj:`InternalPlacement`
           ..
       



   ..
       !! processed by numpydoc !!

   .. py:method:: channel_count() -> int


   .. py:method:: compute_centerline_radius(x: float, contour) -> float
      :abstractmethod:


      
      Return the centerline radius r(x) for this placement.
















      ..
          !! processed by numpydoc !!


   .. py:method:: dtheta_dx(x: float, contour) -> float

      
      Return dθ/dx [rad/m] for the representative lane centerline.
      Default: 0 (axial/vertical channels).
















      ..
          !! processed by numpydoc !!

