pyskyfire.common.EngineNetwork
==============================

.. py:class:: pyskyfire.common.EngineNetwork(stations: Dict[str, Station], signals: Dict[str, float], blocks: List[Block])

   
   Fixed-point engine network orchestrating stations and signals.

   The network manages two dictionaries:

   - **stations**: thermodynamic states as :class:`Station` objects.
   - **signals**: scalar data (e.g., pressure drops, powers, targets, flags).

   Blocks added to the network advertise four metadata lists
   (``station_inputs``, ``station_outputs``, ``signal_inputs``, ``signal_outputs``)
   and implement ``compute(stations, signals) -> (stations_out, signals_out)``.
   The network executes blocks in a chosen order and applies a fixed-point
   iteration until convergence.

   :Parameters:

       **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
           Initial station map.

       **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
           Initial scalar signal map.

       **blocks** : :class:`python:list`\[:obj:`Block`]
           Block instances participating in the network.

   :Attributes:

       **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
           Current station map keyed by name.

       **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
           Current scalar signal map keyed by name.

       **blocks** : :class:`python:list`\[:obj:`Block`]
           Processing blocks in execution order.

       **residuals** : :class:`python:list`\[:class:`python:float`]
           History of per-iteration maximum relative residuals.

       **block_results** : :class:`python:dict`\[:class:`python:str`, :class:`python:dict`\[:class:`python:str`, :obj:`Any`]]
           Optional post-process results keyed by block name, populated after
           convergence.









   .. seealso::

       
       :obj:`FluidBlock`
           Base class for blocks that transform fluid stations.
       :obj:`SignalBlock`
           Base class for blocks that operate on scalar signals only.
       
       
   .. rubric:: Notes

   The default execution order preserves the given block list. A dependency-
   aware topological sort can replace it; see :meth:`_toposort`.



   ..
       !! processed by numpydoc !!

   .. py:method:: _merge_and_residual(store: Dict[str, Any], delta: Dict[str, Any], res: float) -> float

      
      Merge updates into a store and track the max relative residual.

      For each key in ``delta``:

      - If values are :class:`Station`, compute componentwise relative changes
        for ``p``, ``T`` and (if both present) ``mdot``; update ``res`` with
        the maximum.
      - If values are scalars, compute a relative change and update ``res``.
      - Overwrite/write the value in ``store``.

      :Parameters:

          **store** : :class:`python:dict`\[:class:`python:str`, :obj:`Any`]
              Destination dictionary to update in-place.

          **delta** : :class:`python:dict`\[:class:`python:str`, :obj:`Any`]
              New values produced by a block this pass.

          **res** : :class:`python:float`
              Current running maximum relative residual.



      :Returns:

          :class:`python:float`
              Updated maximum relative residual.








      .. rubric:: Notes

      To avoid division by zero, denominators are clamped with small epsilons:
      ``1e-10`` for stations and ``1e-8`` for scalars. If either ``mdot`` is
      ``NaN``, the mass-flow term is skipped for that key.



      ..
          !! processed by numpydoc !!


   .. py:method:: _toposort(blocks: Iterable[Block]) -> List[Block]

      
      Compute an execution order for blocks.

      The default implementation preserves the given order. Replace with a
      dependency-aware topological sort (e.g., Kahn’s algorithm) using each
      block’s ``station_*`` and ``signal_*`` metadata for strict ordering.

      :Parameters:

          **blocks** : :obj:`Iterable`\[:obj:`Block`]
              Blocks to order.



      :Returns:

          :class:`python:list`\[:obj:`Block`]
              Execution order for the fixed-point sweep.








      .. rubric:: Notes

      Sorting is currently not implemented; the user is responsible for
      providing a solvable order.



      ..
          !! processed by numpydoc !!


   .. py:method:: run_fixed_point(tol: float = 1e-06, max_iter: int = 100)

      
      Iterate fixed-point sweeps until convergence or iteration cap.

      Repeatedly calls :meth:`update` until the maximum relative residual falls
      below ``tol`` or until ``max_iter`` sweeps have been performed.

      :Parameters:

          **tol** : :class:`python:float`, :obj:`optional`
              Convergence threshold on max relative residual. Default is ``1e-6``.

          **max_iter** : :class:`python:int`, :obj:`optional`
              Maximum number of sweeps. Default is ``100``.







      :Raises:

          :obj:`RuntimeError`
              If the network does not converge within ``max_iter``.




      .. rubric:: Notes

      Appends per-iteration residuals to :attr:`residuals`. On convergence,
      performs a post-process sweep calling each block’s :meth:`post_process`
      and stores any non-empty results in :attr:`block_results`.



      ..
          !! processed by numpydoc !!


   .. py:method:: update() -> float

      
      Execute one full sweep over all blocks.

      Calls each block’s :meth:`compute` with the current ``stations`` /
      ``signals``, merges its outputs, and tracks the maximum relative change
      across all updated entries.




      :Returns:

          :class:`python:float`
              Maximum relative residual observed this sweep.











      ..
          !! processed by numpydoc !!

