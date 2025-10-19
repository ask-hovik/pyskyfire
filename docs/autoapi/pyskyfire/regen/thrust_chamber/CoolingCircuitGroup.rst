pyskyfire.regen.thrust_chamber.CoolingCircuitGroup
==================================================

.. py:class:: pyskyfire.regen.thrust_chamber.CoolingCircuitGroup(circuit_list, configuration=None)

   
   Collection of :class:`CoolingCircuit` objects forming the full cooling system.

   Used by :class:`ThrustChamber` to coordinate multi-segment cooling
   layouts, manage overlap, and query active channels along the contour.

   :Parameters:

       **circuit_list** : :class:`python:list`\[:obj:`CoolingCircuit`]
           List of circuit instances.

       **configuration** : :class:`python:str`, :obj:`optional`
           Optional label or configuration identifier.

   :Attributes:

       **circuits** : :class:`python:list`\[:obj:`CoolingCircuit`]
           All active cooling circuits.









   .. seealso::

       
       :obj:`CoolingCircuit`
           ..
       :obj:`ThrustChamber`
           ..
       



   ..
       !! processed by numpydoc !!

   .. py:method:: number_of_channels(x, *, occluding_only=False)

      
      Return the total number of active channels at position `x`.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].

          **occluding_only** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
              If True, count only circuits that occlude the wall surface.



      :Returns:

          :class:`python:int`
              Total number of channels currently active at `x`.











      ..
          !! processed by numpydoc !!

