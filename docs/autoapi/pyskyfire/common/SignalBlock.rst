pyskyfire.common.SignalBlock
============================

.. py:class:: pyskyfire.common.SignalBlock

   Bases: :py:obj:`abc.ABC`

   .. autoapi-inheritance-diagram:: pyskyfire.common.SignalBlock
      :parts: 1
      :private-bases:


   
   Abstract base for blocks that operate on **scalar signals** only.



   :Attributes:

       **station_inputs** : :class:`python:list`\[:class:`python:str`]
           Always empty for pure signal blocks.

       **station_outputs** : :class:`python:list`\[:class:`python:str`]
           Always empty for pure signal blocks.

       **signal_inputs** : :class:`python:list`\[:class:`python:str`]
           Names of scalar signal keys read by the block.

       **signal_outputs** : :class:`python:list`\[:class:`python:str`]
           Names of scalar signal keys written by the block.









   .. seealso::

       
       :obj:`FluidBlock`
           Base class for blocks that read/write fluid *stations*.
       
       



   ..
       !! processed by numpydoc !!

   .. py:method:: compute(stations_in: dict[str, pyskyfire.common.engine_network.Station], signals_in: dict[str, float]) -> tuple[dict[str, pyskyfire.common.engine_network.Station], dict[str, float]]
      :abstractmethod:


      
      Compute one pass of the block (signals only).


      :Parameters:

          **stations_in** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Unused for pure signal blocks (present for interface uniformity).

          **signals_in** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Input scalar signals keyed by name.



      :Returns:

          **stations_out** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Always empty for pure signal blocks.

          **signals_out** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Signals produced by the block.











      ..
          !! processed by numpydoc !!


   .. py:method:: post_process(stations: dict[str, pyskyfire.common.engine_network.Station], signals: dict[str, float]) -> dict[str, any]

      
      Optional finalization hook run **after convergence**.


      :Parameters:

          **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Final converged stations (unused).

          **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Final converged scalar signals (unused).



      :Returns:

          :class:`python:dict`\[:class:`python:str`, :obj:`Any`]
              Arbitrary post-processed results. Default is empty.











      ..
          !! processed by numpydoc !!

