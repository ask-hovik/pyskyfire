pyskyfire.skycea.EquilibriumNozzle
==================================

.. py:class:: pyskyfire.skycea.EquilibriumNozzle(inp: EquilNozzleInputs)

   .. py:method:: _a_eq(T, p)

      
      Equilibrium speed of sound a_eq = sqrt((∂p/∂ρ)_s,eq).
      Approximate by small symmetric isentropic TP steps:
      - Move slightly in p, adjust T to keep s = s0 (via 1D root), re-equilibrating each time;
      - Compute dp/dρ at constant s from two neighboring points.
















      ..
          !! processed by numpydoc !!


   .. py:method:: _equilibrate_TP(T, p)


   .. py:method:: _set_elements_by_mixing_rule()

      
      Construct an arbitrary species mixture that matches the requested element totals,
      then let equilibrium find a consistent composition. If you already have propellant
      streams, replace this with your mixing logic and equilibrium at (T_c, p_c).
















      ..
          !! processed by numpydoc !!


   .. py:method:: solve_exit_at_area_ratio(Ae_over_At: float)


   .. py:method:: solve_throat()

      
      Unknowns: (T*, p*) with constraints:
        1) s(T*,p*,eq) = s0
        2) h0 - h(T*,p*,eq) - a_eq(T*,p*)^2 / 2 = 0   [since M=1 ⇒ V = a_eq]
















      ..
          !! processed by numpydoc !!


   .. py:method:: thrust(exit_state)

      
      1-D thrust (no loss model): F = mdot*V_e + (p_e - p_amb) A_e
















      ..
          !! processed by numpydoc !!

