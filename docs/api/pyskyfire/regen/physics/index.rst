pyskyfire.regen.physics
=======================

.. py:module:: pyskyfire.regen.physics






Module Contents
---------------

.. py:function:: h_gas_bartz_enthalpy_driven(k_gr, D_hyd, Cp_gr, mu_gr, mdot_g, A_chmb, T_g, T_gr)

   Bartz correlation for hot gas side heat transfer coefficient

   :param k_gr is the:

   Subscripts:
       g denotes free stream gas properties
       r denotes reference enthalpy conditions which are averaged between the free stream and wall metal conditions


.. py:function:: h_gas_bartz(D_t, mu_g, cp_g, Pr_g, p_c, c_star, A_t, A_x, sigma)

   Bartz correlation from wikipedia without the curvature correction


.. py:function:: sigma(T_hw, T_c, gamma_g, M_g, omega)

   Dimensionless parameter accounting for variation of gas properties across the boundary layer


.. py:function:: h_coolant_colburn(k_cf, D_c, Cp_cr, mu_cf, mdot_c, A_c, phi_curv=1)

   Colburn correlation for coolant side heat transfer coefficient

   :param k_Cf is the conductivity in the coolant film:
   :param D_C is the coolant tube hydraulic diameter:
   :param cp_cr is then maybe the averaged cp between the film condition and the wall metal condition:
   :param mu_Cf is the viscosity at the film conditions:
   :param mdot_c is the coolant mass flow rate:
   :param A_c is the cross sectional area of coolant flow:

   Subscripts:
       Cf denotes film coolant film condition
       c denotes bulk coolant conditions
       r denotes reference enthalpy condition which are averaged between the free stream and wall metal conditions


.. py:function:: u_coolant(rho, mdot_c_single_channel, A_cool)

   coolant velocity as a function of x


.. py:function:: reynolds(rho, u, L, mu)

.. py:function:: phi_curv(Re_c, D_c, R_curv)

   Curvature effect on the cold side heat transfer coefficient


.. py:data:: ReDh_laminar
   :value: 2300


.. py:data:: ReDh_turbulent
   :value: 3500


.. py:function:: f_darcy_laminar(ReDh, Dh, x)

.. py:function:: f_darcy_turbulent(ReDh, Dh, x, roughness)

.. py:function:: f_darcy(ReDh, Dh, x, roughness)

.. py:function:: T_aw(gamma, M_inf, T_inf, Pr)

   Adiabatic wall temperature, sometimes called recovery temperature T_r in some sources
   TODO: Implement a more robust version that could also handle laminar flow, in case this function is suddenly
   used in a function where this makes sense.


