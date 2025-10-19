pyskyfire.common.blocks.MassFlowMergerBlock
===========================================

.. py:class:: pyskyfire.common.blocks.MassFlowMergerBlock(name: str, st_in: list[str], st_out: str, medium)

   Bases: :py:obj:`FluidBlock`

   .. autoapi-inheritance-diagram:: pyskyfire.common.blocks.MassFlowMergerBlock
      :parts: 1
      :private-bases:


   
   Merge multiple inlet streams into a single outlet.

   Assumes identical fluid species and mixes by enthalpy at a reference
   pressure equal to the minimum inlet pressure.

   :Parameters:

       **name** : :class:`python:str`
           Block name.

       **st_in** : :class:`python:list`\[:class:`python:str`]
           Inlet station keys.

       **st_out** : :class:`python:str`
           Outlet station key.

       **medium** : :class:`python:str`
           CoolProp fluid string.

   :Attributes:

       **station_inputs, station_outputs, signal_inputs, signal_outputs** : :class:`python:list`\[:class:`python:str`]
           Network I/O metadata.













   ..
       !! processed by numpydoc !!

   .. py:method:: compute(stations, signals)

      
      Mix inlet streams by mass and enthalpy to produce the outlet.


      :Parameters:

          **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Network stations (must contain all in ``st_in``).

          **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Scalar signals (unused).



      :Returns:

          **stations_out** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Outlet station dict.

          **signals_out** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Empty dict.




      :Raises:

          :obj:`ValueError`
              If the merged total mass flow is non-positive.




      .. rubric:: Notes

      .. Note::
          Property calls use a small temperature offset (``-1e-3 K``) to avoid saturation issues in CoolProp.



      ..
          !! processed by numpydoc !!

