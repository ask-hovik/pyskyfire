pyskyfire.common.blocks.TransmissionBlock
=========================================

.. py:class:: pyskyfire.common.blocks.TransmissionBlock(name, sink_keys, out_key='P_required')

   Bases: :py:obj:`SignalBlock`

   .. autoapi-inheritance-diagram:: pyskyfire.common.blocks.TransmissionBlock
      :parts: 1
      :private-bases:


   
   Sum multiple shaft-power signals into one required-power signal.

   Useful for aggregating pump/compressor power demands prior to a turbine
   block that must supply the total shaft power.

   :Parameters:

       **name** : :class:`python:str`
           Block name.

       **sink_keys** : :class:`python:list`\[:class:`python:str`]
           Input power signal keys to sum.

       **out_key** : :class:`python:str`, :obj:`optional`
           Output key for the total required power (default ``"P_required"``).

   :Attributes:

       **station_inputs, station_outputs, signal_inputs, signal_outputs** : :class:`python:list`\[:class:`python:str`]
           Network I/O metadata.













   ..
       !! processed by numpydoc !!

   .. py:method:: compute(st, sg)

      
      Sum input power signals.


      :Parameters:

          **st** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Stations dict (unused).

          **sg** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Signals dict containing all ``sink_keys`` [W].



      :Returns:

          **stations_out** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Empty dict.

          **signals_out** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              ``{out_key: total_power}`` in watts.











      ..
          !! processed by numpydoc !!

