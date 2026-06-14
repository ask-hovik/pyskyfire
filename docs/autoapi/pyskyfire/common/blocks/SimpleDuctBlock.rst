pyskyfire.common.blocks.SimpleDuctBlock
=======================================

.. py:class:: pyskyfire.common.blocks.SimpleDuctBlock(name: str, st_in: str, st_out: str, pressure_ratio: float, medium: str)

   Bases: :py:obj:`FluidBlock`

   .. autoapi-inheritance-diagram:: pyskyfire.common.blocks.SimpleDuctBlock
      :parts: 1
      :private-bases:


   
   Adiabatic homogeneous duct loss model with fixed pressure ratio.

   Applies a fixed pressure ratio via

       p_out = pressure_ratio * p_in

   and computes the outlet temperature from constant specific enthalpy:

       h_out = h_in

   This is a better approximation for an adiabatic low-Mach duct or
   concentrated pressure-loss element than holding temperature constant.












   .. rubric:: Notes

   This is not an isentropic expansion. Entropy increases. For real fluids,
   especially near saturation or the critical point, T_out may differ strongly
   from T_in.



   ..
       !! processed by numpydoc !!

   .. py:method:: compute(stations: Dict[str, pyskyfire.common.engine_network.Station], signals: Dict[str, float]) -> tuple[Dict[str, pyskyfire.common.engine_network.Station], Dict[str, float]]

      
      Compute one pass of the block.


      :Parameters:

          **stations_in** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Input network stations keyed by name.

          **signals_in** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Input scalar signals keyed by name.



      :Returns:

          **stations_out** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Updated or newly created stations produced by the block.

          **signals_out** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Updated or newly created scalar signals produced by the block.








      .. rubric:: Notes

      .. Note::
          Implementations should treat missing inputs as an error. The engine network orchestrator is expected to provide the declared inputs.



      ..
          !! processed by numpydoc !!

