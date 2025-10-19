pyskyfire.regen.solver
======================

.. py:module:: pyskyfire.regen.solver


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/regen/solver/BoundaryConditions
   /autoapi/pyskyfire/regen/solver/HeatExchangerPhysics

.. autoapisummary::

   pyskyfire.regen.solver.BoundaryConditions
   pyskyfire.regen.solver.HeatExchangerPhysics


Functions
---------

.. autoapisummary::

   pyskyfire.regen.solver.analyse_residuals
   pyskyfire.regen.solver.solve_heat_exchanger_euler
   pyskyfire.regen.solver.steady_heating_analysis


Module Contents
---------------

.. py:function:: analyse_residuals(residual_log, n_cells, p=2)

   
   Aggregate local solver residuals into global history and final per-cell vector.


   :Parameters:

       **residual_log** : :class:`python:list` or :data:`python:None`
           List of tuples ``(cell, iter, R1, R2)`` recorded during solves.
           If ``None`` or empty, returns ``(None, None)``.

       **n_cells** : :class:`python:int`
           Number of axial cells.

       **p** : :class:`python:int` or :class:`python:float`, :obj:`optional`
           Norm order for global residual history: ``2`` for RMS,
           ``np.inf`` for :math:`L_\infty`, etc. Default is 2.



   :Returns:

       **history** : :obj:`ndarray <numpy.ndarray>` or :data:`python:None`
           Global residual norm for iterations ``0..k_max``, or ``None``.

       **final_per_cell** : :obj:`ndarray <numpy.ndarray>` or :data:`python:None`
           Final residual magnitude per cell at its last local iteration, or ``None``.











   ..
       !! processed by numpydoc !!

.. py:function:: solve_heat_exchanger_euler(thrust_chamber, boundary_conditions, n_nodes, circuit_index, output, log_residuals=True)

   
   Solve 1-D steady heating with a marching Euler scheme.


   :Parameters:

       **thrust_chamber** : :obj:`Any`
           Chamber model exposing geometry and property methods.

       **boundary_conditions** : :obj:`BoundaryConditions`
           Inlet temperature/pressure/mass-flow conditions.

       **n_nodes** : :class:`python:int`
           Number of axial nodes.

       **circuit_index** : :class:`python:int`
           Which cooling circuit to simulate.

       **output** : :ref:`bool <python:bltin-boolean-values>`
           If True, print progress to stdout.

       **log_residuals** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           If True, record local residuals each Newton iteration per cell.



   :Returns:

       :class:`python:dict`
           Results with keys:
           
           - ``x`` : axial coordinates [m]
           - ``T`` : temperatures array, shape ``(n_nodes, 1 + n_walls + 1)``
             (coolant + reversed wall interfaces from cold→hot)
           - ``T_static`` : coolant static temperature [K]
           - ``T_stagnation`` : coolant stagnation temperature [K]
           - ``p_static`` : coolant static pressure [Pa]
           - ``p_stagnation`` : coolant stagnation pressure [Pa]
           - ``dQ_dA`` : local heat flux [W m⁻²]
           - ``velocity`` : coolant velocity [m s⁻¹]
           - ``residuals`` : tuple of (global history, final-per-cell) or (None, None)











   ..
       !! processed by numpydoc !!

.. py:function:: steady_heating_analysis(thrust_chamber, boundary_conditions, n_nodes=100, circuit_index=0, solver='newton', output=True)

   
   Run the steady heating analysis.


   :Parameters:

       **thrust_chamber** : :obj:`Any`
           Chamber model exposing the required geometry & property APIs.

       **boundary_conditions** : :obj:`BoundaryConditions`
           Coolant inlet boundary conditions.

       **n_nodes** : :class:`python:int`, :obj:`optional`
           Number of axial nodes. Default is 100.

       **circuit_index** : :class:`python:int`, :obj:`optional`
           Cooling-circuit index. Default is 0.

       **solver** : {'newton'}, :obj:`optional`
           Solver selector. Currently only ``'newton'`` is implemented.

       **output** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           If True, print progress. Default is True.



   :Returns:

       :class:`python:dict`
           See :func:`solve_heat_exchanger_euler` for keys.











   ..
       !! processed by numpydoc !!

