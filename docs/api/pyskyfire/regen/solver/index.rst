pyskyfire.regen.solver
======================

.. py:module:: pyskyfire.regen.solver






Module Contents
---------------

.. py:class:: BoundaryConditions(T_coolant_in, p_coolant_in, mdot_coolant)

   Object for storing boundary conditions for the solver.

   :param T_coolant_in: Static(?) temperature of coolant at cooling channel inlet (K)
   :type T_coolant_in: float
   :param p_coolant_in: Static(?) pressure of coolant at cooling channel inlet (Pa)
   :type p_coolant_in: float
   :param mdot_coolant: mass flow rate of coolant through cooling channel inlet (kg/s)
   :type mdot_coolant: float


   .. py:attribute:: T_coolant_in


   .. py:attribute:: p_coolant_in


   .. py:attribute:: mdot_coolant


.. py:class:: HeatExchangerPhysics(thrust_chamber, circuit_index)

   Encapsulates the physical calculations for the heat exchanger.
   This class is responsible for computing the heat fluxes and rates based
   on the given engine properties and operating conditions.


   .. py:attribute:: thrust_chamber


   .. py:attribute:: circuit_index


   .. py:attribute:: counter
      :value: 0



   .. py:method:: dQ_hot_dx(x, T_hw)

      Computes the heat transfer rate from the hot gas to the wall. Bartz equation, for example



   .. py:method:: dQ_cond_dx(x, T_hw, T_cw)


   .. py:method:: dQ_cold_dx(x, T_cw, T_cool)

      Computes the heat transfer rate from the wall to the coolant.



   .. py:method:: coolant_temperature_rate(T_cool, p_cool, dQ_cold_dx)

      Computes the rate of change of the coolant temperature.



   .. py:method:: coolant_pressure_rate(x, T_cool, p_cool)

      Computes the rate of change of the coolant pressure due to frictional losses.



   .. py:method:: interface_temperatures(x, T_hw, T_cw)

      Returns a list [T0, T1, ..., Tn] of wall-stack temperatures,
      where T0=T_hw, Tn=T_cw, and in between are each wall interface.



.. py:function:: solve_heat_exchanger_euler(thrust_chamber, boundary_conditions, n_nodes, circuit_index, output, log_residuals=True)

   Solve the 1D steady-state heat exchanger from x=0 to x=x_domain[-1].

   :param engine: your engine object containing geometry & property methods
   :param n_nodes: number of axial nodes along the thrust chamber

   :returns:    "x"          -> 1D array of axial positions
                "T"          -> 2D array of temperatures, shape (n_nodes, 3):
                                columns = [T_coolant, T_wall_cold_side, T_wall_hot_side]
                "T_coolant"  -> 1D array of coolant temperatures
                "p_coolant"  -> 1D array of coolant pressures
                "dQ_dA"      -> 1D array of local heat fluxes (W/m^2)
                "velocity"   -> 1D array of coolant velocities (m/s)
   :rtype: A dictionary with keys


.. py:function:: analyse_residuals(residual_log, n_cells, p=2)

   :param residual_log: The list returned by `solve_channel`.  If None, nothing happens.
   :type residual_log: list | None
   :param n_cells: Number of axial nodes in the simulation.
   :type n_cells: int
   :param p: Order of the global norm: 2 for RMS, np.inf for L∞, etc.
   :type p: int | float

   :returns: * **history** (*(n_iter,) ndarray | None*) -- Global residual norm per iteration 0..k_max.
             * **final_per_cell** (*(n_cells,) ndarray | None*) -- Residual magnitude in each cell at the last local iteration.


.. py:function:: steady_heating_analysis(thrust_chamber, boundary_conditions, n_nodes=100, circuit_index=0, solver='newton', output=True)

   Run the steady heating analysis.

   :param engine: Engine object with geometry, boundary conditions, and physics.
   :param n_nodes: Number of nodes (for Newton) or resolution for post-processing (for Radau).
   :param solver: String, currently only "newton" available

   :returns: A dictionary with simulation results.


