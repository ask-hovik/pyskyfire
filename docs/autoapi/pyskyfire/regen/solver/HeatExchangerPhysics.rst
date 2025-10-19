pyskyfire.regen.solver.HeatExchangerPhysics
===========================================

.. py:class:: pyskyfire.regen.solver.HeatExchangerPhysics(thrust_chamber, circuit_index)

   
   Encapsulate hot-side, wall, and coolant heat-transfer/pressure models.

   This helper evaluates local heat transfer and pressure-loss terms for a
   given thrust-chamber and cooling circuit.

   :Parameters:

       **thrust_chamber** : :obj:`Any`
           Object exposing geometry and property models used by the solver
           (e.g., ``contour``, ``combustion_transport``, ``cooling_circuit_group``,
           ``wall_group``).

       **circuit_index** : :class:`python:int`
           Index of the cooling circuit to use.














   ..
       !! processed by numpydoc !!

   .. py:method:: coolant_pressure_rate(x, T_cool, p_cool)

      
      Axial rates of static and stagnation pressure (friction + area change).


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].

          **T_cool** : :class:`python:float`
              Coolant static temperature [K].

          **p_cool** : :class:`python:float`
              Coolant static pressure [Pa].



      :Returns:

          :class:`python:tuple`\[:class:`python:float`, :class:`python:float`]
              ``(dp_static/dx, dp_stagnation/dx)`` in [Pa m⁻¹].











      ..
          !! processed by numpydoc !!


   .. py:method:: coolant_temperature_rate(T_cool, p_cool, dQ_cold_dx)

      
      Axial rate of change of coolant temperature.


      :Parameters:

          **T_cool** : :class:`python:float`
              Coolant static temperature [K].

          **p_cool** : :class:`python:float`
              Coolant static pressure [Pa].

          **dQ_cold_dx** : :class:`python:float`
              Heat removed by coolant per unit length [W m⁻¹].



      :Returns:

          :class:`python:float`
              ``dT_cool/dx`` [K m⁻¹].











      ..
          !! processed by numpydoc !!


   .. py:method:: dQ_cold_dx(x, T_cw, T_cool)

      
      Coolant-side heat removal per unit length.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].

          **T_cw** : :class:`python:float`
              Coolant-side wall temperature [K].

          **T_cool** : :class:`python:float`
              Bulk coolant temperature [K].



      :Returns:

          :class:`python:float`
              ``dQ_cw/dx`` [W m⁻¹].











      ..
          !! processed by numpydoc !!


   .. py:method:: dQ_cond_dx(x, T_hw, T_cw)

      
      Conduction heat flow through the wall stack per unit length.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].

          **T_hw** : :class:`python:float`
              Hot-side wall temperature [K].

          **T_cw** : :class:`python:float`
              Coolant-side wall temperature [K].



      :Returns:

          :class:`python:float`
              ``dQ_cond/dx`` [W m⁻¹].








      .. rubric:: Notes

      Treats each wall as a 1-D resistor in series:
      :math:`R_j = L_j / (k_j A)` with ``A = dA_dx_hot`` per unit length.



      ..
          !! processed by numpydoc !!


   .. py:method:: dQ_hot_dx(x, T_hw)

      
      Hot-side heat input per unit length using Bartz-style correlation.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].

          **T_hw** : :class:`python:float`
              Hot-side wall temperature [K].



      :Returns:

          :class:`python:float`
              ``dQ_hw/dx`` [W m⁻¹], positive when heating the wall.








      .. rubric:: Notes

      Implements an **enthalpy-driven** Bartz form. The gas-side coefficient
      is evaluated with property data drawn from the chamber model at ``x``.



      ..
          !! processed by numpydoc !!


   .. py:method:: interface_temperatures(x, T_hw, T_cw)

      
      Wall-interface temperatures across the stack at position ``x``.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].

          **T_hw** : :class:`python:float`
              Hot-side wall temperature [K].

          **T_cw** : :class:`python:float`
              Coolant-side wall temperature [K].



      :Returns:

          :class:`python:list`\[:class:`python:float`]
              Temperatures ``[T_hot, T_1, ..., T_cold]`` across interfaces.











      ..
          !! processed by numpydoc !!

