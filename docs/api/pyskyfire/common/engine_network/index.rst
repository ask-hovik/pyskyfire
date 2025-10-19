pyskyfire.common.engine_network
===============================

.. py:module:: pyskyfire.common.engine_network




Module Contents
---------------

.. py:class:: Station

   .. py:attribute:: p
      :type:  float


   .. py:attribute:: T
      :type:  float


   .. py:attribute:: mdot
      :type:  float


.. py:class:: EngineNetwork(stations, signals, blocks)

   Fixed-point network that stores *two* dictionaries:
       • stations : thermodynamic states        (Station objects)
       • signals  : scalar data (Δp, powers, targets, flags …)
   Every Block advertises four lists:
       station_inputs,  station_outputs,
       signal_inputs,   signal_outputs
   and its compute() returns (stations_out, signals_out).


   .. py:attribute:: stations


   .. py:attribute:: signals


   .. py:attribute:: blocks


   .. py:attribute:: residuals
      :type:  List[float]
      :value: []



   .. py:method:: update()

      Executes one full sweep over all blocks.
      Returns the maximum relative change (residual) seen.



   .. py:method:: run_fixed_point(tol = 1e-06, max_iter = 100)

      Iterate until the largest relative change among *all* advertised
      outputs is below *tol* or until *max_iter* is reached.



