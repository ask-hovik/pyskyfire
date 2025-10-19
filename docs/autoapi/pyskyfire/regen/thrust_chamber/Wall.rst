pyskyfire.regen.thrust_chamber.Wall
===================================

.. py:class:: pyskyfire.regen.thrust_chamber.Wall(material, thickness, name=None)

   
   Single structural wall or coating layer in a thrust chamber.

   Encapsulates material and thickness information for heat-conduction
   calculations.

   :Parameters:

       **material** : :obj:`Material`
           Material object providing thermal conductivity and other data.

       **thickness** : :class:`python:float` or :func:`python:callable`
           Constant thickness [m] or function `t(x)` returning thickness.

       **name** : :class:`python:str`, :obj:`optional`
           Descriptive label of the layer (e.g. “Copper liner”).

   :Attributes:

       **material** : :obj:`Material`
           Thermal material definition.

       **_thickness** : :class:`python:float` or :func:`python:callable`
           Underlying storage for the thickness definition.

       **name** : :class:`python:str` or :data:`python:None`
           Optional descriptive name.









   .. seealso::

       
       :obj:`WallGroup`
           Container combining multiple walls.
       
       



   ..
       !! processed by numpydoc !!

   .. py:method:: thickness(x)

      
      Return the wall thickness at position `x`.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].



      :Returns:

          :class:`python:float`
              Wall thickness [m].











      ..
          !! processed by numpydoc !!

