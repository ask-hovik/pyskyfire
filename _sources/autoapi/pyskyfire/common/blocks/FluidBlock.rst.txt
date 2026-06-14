pyskyfire.common.blocks.FluidBlock
==================================

.. py:class:: pyskyfire.common.blocks.FluidBlock(medium: str)

   Bases: :py:obj:`abc.ABC`

   .. autoapi-inheritance-diagram:: pyskyfire.common.blocks.FluidBlock
      :parts: 1
      :private-bases:


   
   Abstract base for blocks that transform **fluid** states.

   Each :class:`FluidBlock` consumes and/or produces :class:`Station`
   objects and may also consume or emit scalar *signals*
   (e.g., power, pressure drop).

   :Parameters:

       **medium** : :class:`python:str`
           CoolProp fluid identifier used for property calls
           (e.g., ``"N2O"`` or ``"Water"``).

   :Attributes:

       **station_inputs** : :class:`python:list`\[:class:`python:str`]
           Names of input station keys expected by the block.

       **station_outputs** : :class:`python:list`\[:class:`python:str`]
           Names of output station keys produced by the block.

       **signal_inputs** : :class:`python:list`\[:class:`python:str`]
           Names of scalar signal keys read by the block.

       **signal_outputs** : :class:`python:list`\[:class:`python:str`]
           Names of scalar signal keys written by the block.

       **medium** : :class:`python:str`
           CoolProp fluid string used by this instance.









   .. seealso::

       
       :class:`~pyskyfire.common.engine_network.Station`
           Thermodynamic state container used by the network.
       :class:`~pyskyfire.common.blocks.SignalBlock`
           Abstract base for blocks that operate on **scalar signals** only.
       
       
   .. rubric:: Notes

   .. note::
       Subclasses must implement :meth:`compute`. Missing declared inputs
       should be treated as an error; the engine network orchestrator is
       responsible for providing them.



   ..
       !! processed by numpydoc !!

   .. py:method:: compute(stations_in: dict[str, pyskyfire.common.engine_network.Station], signals_in: dict[str, float]) -> tuple[dict[str, pyskyfire.common.engine_network.Station], dict[str, float]]
      :abstractmethod:


      
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


   .. py:method:: post_process(stations: dict[str, pyskyfire.common.engine_network.Station], signals: dict[str, float]) -> dict[str, any]

      
      Optional finalization hook run **after convergence**.


      :Parameters:

          **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Final converged stations.

          **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Final converged scalar signals.



      :Returns:

          :class:`python:dict`\[:class:`python:str`, :obj:`Any`]
              Arbitrary post-processed results to be collected by
              the network (e.g., axial profiles, derived scalars). Default is empty.











      ..
          !! processed by numpydoc !!

