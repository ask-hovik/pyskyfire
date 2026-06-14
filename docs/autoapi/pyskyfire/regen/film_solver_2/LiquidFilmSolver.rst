pyskyfire.regen.film_solver_2.LiquidFilmSolver
==============================================

.. py:class:: pyskyfire.regen.film_solver_2.LiquidFilmSolver(thrust_chamber, boundary_conditions)

   .. py:method:: _evaporation_rate(x: float, Gamma: float, x_inj: float) -> tuple[float, float, float, float]


   .. py:method:: _h0_colburn(x: float, x_eff: float) -> float


   .. py:method:: _radiative_heat_flux(x: float) -> float


   .. py:method:: _transpiration_correction(x: float, h0: float, m_vap: float) -> float


   .. py:method:: solve(x_array: Sequence[float]) -> LiquidFilmResults

