pyskyfire.skycea.nozzle_solver
==============================

.. py:module:: pyskyfire.skycea.nozzle_solver




Module Contents
---------------

.. py:class:: EquilNozzleInputs

   .. py:attribute:: mech
      :type:  str


   .. py:attribute:: chamber_solution
      :type:  Optional[cantera.Solution]
      :value: None



   .. py:attribute:: chamber_T
      :type:  Optional[float]
      :value: None



   .. py:attribute:: chamber_p
      :type:  Optional[float]
      :value: None



   .. py:attribute:: chamber_X
      :type:  Optional[Union[str, Dict[str, float], numpy.ndarray]]
      :value: None



   .. py:attribute:: A_throat
      :type:  float
      :value: 0.0



   .. py:attribute:: A_exit
      :type:  Optional[float]
      :value: None



   .. py:attribute:: p_amb
      :type:  Optional[float]
      :value: None



   .. py:attribute:: fd_eps
      :type:  float
      :value: 0.0001



.. py:class:: EquilibriumNozzle(inp)

   .. py:attribute:: inp


   .. py:attribute:: gas


   .. py:attribute:: h0


   .. py:attribute:: s0


   .. py:attribute:: throat


   .. py:attribute:: mdot


   .. py:method:: solve_throat()

      Unknowns: (T*, p*) with constraints:
        1) s(T*,p*,eq) = s0
        2) h0 - h(T*,p*,eq) - a_eq(T*,p*)^2 / 2 = 0   [since M=1 ⇒ V = a_eq]



   .. py:method:: solve_exit_at_area_ratio(Ae_over_At)


   .. py:method:: thrust(exit_state)

      1-D thrust (no loss model): F = mdot*V_e + (p_e - p_amb) A_e



