pyskyfire.common.blocks.RegenBlock
==================================

.. py:class:: pyskyfire.common.blocks.RegenBlock(name: str, st_in: str, st_out: str, circuit_index: int, thrust_chamber, medium)

   Bases: :py:obj:`FluidBlock`

   .. autoapi-inheritance-diagram:: pyskyfire.common.blocks.RegenBlock
      :parts: 1
      :private-bases:


   
   Regenerative-cooling segment (single circuit).

   Runs a 1-D steady heating analysis for a specified cooling circuit and
   returns the outlet state and the pressure drop across the circuit.

   :Parameters:

       **name** : :class:`python:str`
           Block name; used to form ``dp_{name}``.

       **st_in** : :class:`python:str`
           Inlet station key (coolant entering the circuit).

       **st_out** : :class:`python:str`
           Outlet station key (coolant leaving the circuit).

       **circuit_index** : :class:`python:int`
           Index of the cooling circuit to simulate.

       **thrust_chamber**
           Geometry/physics object passed to the solver.

       **medium** : :class:`python:str`
           CoolProp fluid string for the coolant.

   :Attributes:

       **dp_key** : :class:`python:str`
           Name of the emitted pressure-drop signal.

       **station_inputs, station_outputs, signal_inputs, signal_outputs** : :class:`python:list`\[:class:`python:str`]
           Network I/O metadata.













   ..
       !! processed by numpydoc !!

   .. py:method:: compute(stations: dict[str, pyskyfire.common.engine_network.Station], signals: dict[str, float]) -> tuple[dict[str, pyskyfire.common.engine_network.Station], dict[str, float]]

      
      Run steady heating analysis and emit outlet state and Δp.


      :Parameters:

          **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Network stations (must contain ``st_in``).

          **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Scalar signals (unused).



      :Returns:

          **stations_out** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              ``{st_out: Station}`` at circuit exit.

          **signals_out** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              ``{dp_key: float}`` pressure loss across the circuit [Pa].











      ..
          !! processed by numpydoc !!


   .. py:method:: post_process(stations: dict[str, pyskyfire.common.engine_network.Station], signals: dict[str, float]) -> dict[str, any]

      
      Re-run the solver on a finer grid to collect detailed outputs.


      :Parameters:

          **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Final converged stations.

          **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Final converged scalars.



      :Returns:

          :class:`python:dict`\[:class:`python:str`, :obj:`Any`]
              Solver output dictionary (profiles and scalars) suitable
              for reporting/plotting.











      ..
          !! processed by numpydoc !!

