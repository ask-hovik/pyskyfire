pyskyfire.common.SimpleDuctBlock
================================

.. py:class:: pyskyfire.common.SimpleDuctBlock(name: str, st_in: str, st_out: str, eta: float, medium)

   Bases: :py:obj:`FluidBlock`

   .. autoapi-inheritance-diagram:: pyskyfire.common.SimpleDuctBlock
      :parts: 1
      :private-bases:


   
   Homogeneous duct loss model with fixed efficiency.

   Applies a fixed pressure efficiency :math:`\eta` via
   ``p_out = η * p_in`` while keeping temperature unchanged.

   :Parameters:

       **name** : :class:`python:str`
           Block name.

       **st_in** : :class:`python:str`
           Inlet station key.

       **st_out** : :class:`python:str`
           Outlet station key.

       **eta** : :class:`python:float`
           Pressure efficiency (``0 < η ≤ 1``).

       **medium** : :class:`python:str`
           CoolProp fluid string.

   :Attributes:

       **dp_key** : :class:`python:str`
           Name of the emitted pressure-drop signal.

       **station_inputs, station_outputs, signal_inputs, signal_outputs** : :class:`python:list`\[:class:`python:str`]
           Network I/O metadata.






   :Raises:

       :obj:`ValueError`
           If ``eta`` is not in ``(0, 1]``.







   ..
       !! processed by numpydoc !!

   .. py:method:: compute(stations: Dict[str, pyskyfire.common.engine_network.Station], signals: Dict[str, float]) -> tuple[Dict[str, pyskyfire.common.engine_network.Station], Dict[str, float]]

      
      Apply fixed loss and emit outlet station and Δp.


      :Parameters:

          **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Network stations (must contain ``st_in``).

          **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Scalar signals (unused).



      :Returns:

          **stations_out** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              ``{st_out: Station}`` with ``p_out = η * p_in`` and unchanged ``T``.

          **signals_out** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              ``{dp_key: float}`` with ``dp = p_in - p_out`` [Pa].











      ..
          !! processed by numpydoc !!

