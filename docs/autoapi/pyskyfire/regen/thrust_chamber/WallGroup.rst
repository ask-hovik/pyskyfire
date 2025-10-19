pyskyfire.regen.thrust_chamber.WallGroup
========================================

.. py:class:: pyskyfire.regen.thrust_chamber.WallGroup(walls=None)

   
   Container holding multiple :class:`Wall` layers.

   Provides cumulative quantities such as total wall thickness along the
   engine contour.

   :Parameters:

       **walls** : :class:`python:list`\[:obj:`Wall`], :obj:`optional`
           Collection of wall layers. Defaults to an empty list.

   :Attributes:

       **walls** : :class:`python:list`\[:obj:`Wall`]
           Managed list of wall objects.









   .. seealso::

       
       :obj:`ThrustChamber`
           Uses a `WallGroup` to compute wall-conduction resistance.
       
       



   ..
       !! processed by numpydoc !!

   .. py:method:: total_thickness(x)

      
      Compute total wall thickness at axial position `x`.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].



      :Returns:

          :class:`python:float`
              Sum of all wall thicknesses [m].











      ..
          !! processed by numpydoc !!

