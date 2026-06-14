pyskyfire.common.blocks.TurbineBlock
====================================

.. py:class:: pyskyfire.common.blocks.TurbineBlock(name: str, st_in: str, st_out: str, P_req_key: str, eta: float, medium: str)

   Bases: :py:obj:`FluidBlock`

   .. autoapi-inheritance-diagram:: pyskyfire.common.blocks.TurbineBlock
      :parts: 1
      :private-bases:


   
   Single-stage turbine model.

   Consumes a required shaft-power signal and computes an outlet state
   that provides that work at given efficiency.

   :Parameters:

       **name** : :class:`python:str`
           Block name.

       **st_in** : :class:`python:str`
           Inlet station key.

       **st_out** : :class:`python:str`
           Outlet station key.

       **P_req_key** : :class:`python:str`
           Name of the required power signal [W].

       **eta** : :class:`python:float`
           Turbine isentropic efficiency (0..1).

       **medium** : :class:`python:str`
           CoolProp fluid string.

   :Attributes:

       **dp_key** : :class:`python:str`
           Name of the emitted pressure-drop signal.

       **station_inputs, station_outputs, signal_inputs, signal_outputs** : :class:`python:list`\[:class:`python:str`]
           Network I/O metadata.









   .. seealso::

       
       :obj:`TransmissionBlock`
           Aggregates power demands into the required shaft power.
       
       



   ..
       !! processed by numpydoc !!

   .. py:method:: compute(stations: Dict[str, pyskyfire.common.engine_network.Station], signals: Dict[str, float]) -> tuple[Dict[str, pyskyfire.common.engine_network.Station], Dict[str, float]]

      
      Compute outlet state delivering the requested shaft power.


      :Parameters:

          **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Network stations (must contain ``st_in``).

          **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Scalar signals (must contain ``P_req_key``).



      :Returns:

          **stations_out** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              ``{st_out: Station}`` with updated pressure and temperature.

          **signals_out** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              ``{dp_key: float}`` pressure drop across the turbine [Pa].




      :Raises:

          :obj:`ValueError`
              If inlet mass flow is non-positive.







      ..
          !! processed by numpydoc !!

