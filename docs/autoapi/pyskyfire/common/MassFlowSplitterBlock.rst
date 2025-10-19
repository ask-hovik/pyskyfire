pyskyfire.common.MassFlowSplitterBlock
======================================

.. py:class:: pyskyfire.common.MassFlowSplitterBlock(name: str, st_in: str, st_out: list[str], medium, fractions: list[float] | None = None, frac_keys: list[str] | None = None)

   Bases: :py:obj:`FluidBlock`

   .. autoapi-inheritance-diagram:: pyskyfire.common.MassFlowSplitterBlock
      :parts: 1
      :private-bases:


   
   Split a single inlet stream into multiple outlet branches.

   Fractions can be fixed (sum to 1.0) or read each pass from scalar
   signals to enable dynamic splits.

   :Parameters:

       **name** : :class:`python:str`
           Block name.

       **st_in** : :class:`python:str`
           Inlet station key.

       **st_out** : :class:`python:list`\[:class:`python:str`]
           Outlet station keys.

       **medium** : :class:`python:str`
           CoolProp fluid string.

       **fractions** : :class:`python:list`\[:class:`python:float`] | :data:`python:None`, :obj:`optional`
           Fixed fractions (sum to 1.0) if provided.

       **frac_keys** : :class:`python:list`\[:class:`python:str`] | :data:`python:None`, :obj:`optional`
           Signal keys providing fractions at runtime.

   :Attributes:

       **station_inputs, station_outputs, signal_inputs, signal_outputs** : :class:`python:list`\[:class:`python:str`]
           Network I/O metadata.






   :Raises:

       :obj:`ValueError`
           If both or neither of ``fractions`` and ``frac_keys`` are given,
           or if fixed ``fractions`` do not sum to 1.0.







   ..
       !! processed by numpydoc !!

   .. py:method:: compute(stations, signals)

      
      Create one outlet station per branch based on split fractions.


      :Parameters:

          **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Network stations (must contain ``st_in``).

          **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Scalar signals supplying fractions when ``frac_keys`` is used.



      :Returns:

          **stations_out** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Output stations for each branch.

          **signals_out** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Empty dict (no signals produced).




      :Raises:

          :obj:`ValueError`
              If dynamic fractions do not sum to 1.0.







      ..
          !! processed by numpydoc !!

