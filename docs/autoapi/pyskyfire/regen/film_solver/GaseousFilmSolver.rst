pyskyfire.regen.film_solver.GaseousFilmSolver
=============================================

.. py:class:: pyskyfire.regen.film_solver.GaseousFilmSolver(thrust_chamber, boundary_conditions, circuit_index: int = 0)

   .. py:method:: _convective_h_gaseous(x: float, Mbl: float) -> float


   .. py:method:: _derivatives(x: float, Mbl: float, T_aw: float) -> tuple[float, float]


   .. py:method:: _initial_conditions(x_i: float, Mc_total: float) -> tuple[float, float]


   .. py:method:: _local_gas_flux(x: float) -> float


   .. py:method:: _radiative_heat_flux(x: float) -> float


   .. py:method:: _recovery_temperature(x: float) -> float


   .. py:method:: solve(x_array: Sequence[float], x_dryout: float, Mc_total: float) -> GaseousFilmResults

