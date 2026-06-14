pyskyfire.regen.physics
=======================

.. py:module:: pyskyfire.regen.physics


Attributes
----------

.. autoapisummary::

   pyskyfire.regen.physics.ReDh_laminar
   pyskyfire.regen.physics.ReDh_turbulent


Functions
---------

.. autoapisummary::

   pyskyfire.regen.physics.T_aw
   pyskyfire.regen.physics.f_darcy
   pyskyfire.regen.physics.f_darcy_laminar
   pyskyfire.regen.physics.f_darcy_turbulent
   pyskyfire.regen.physics.h_coolant_colburn
   pyskyfire.regen.physics.h_gas_bartz
   pyskyfire.regen.physics.h_gas_bartz_enthalpy_driven
   pyskyfire.regen.physics.phi_curv
   pyskyfire.regen.physics.reynolds
   pyskyfire.regen.physics.sigma
   pyskyfire.regen.physics.u_coolant


Module Contents
---------------

.. py:function:: T_aw(gamma, M_inf, T_inf, Pr)

   
   Compute the adiabatic-wall (recovery) temperature.


   :Parameters:

       **gamma** : :class:`python:float`
           Ratio of specific heats.

       **M_inf** : :class:`python:float`
           Freestream Mach number.

       **T_inf** : :class:`python:float`
           Freestream static temperature [K].

       **Pr** : :class:`python:float`
           Prandtl number.



   :Returns:

       :class:`python:float`
           Adiabatic-wall temperature [K].








   .. rubric:: Notes

   The recovery factor ``r`` is approximated as ``Pr^(1/3)``,
   valid for turbulent flow. For laminar conditions, a smaller
   recovery factor should be used (≈ Pr⁰˙⁵).

   TODO: Implement a more robust version that could also handle laminar flow, in case this function is suddenly
   used in a function where this makes sense.



   ..
       !! processed by numpydoc !!

.. py:function:: f_darcy(ReDh, Dh, x, roughness)

   
   Composite laminar–turbulent Darcy friction factor with smooth transition.


   :Parameters:

       **ReDh** : :class:`python:float`
           Hydraulic Reynolds number.

       **Dh** : :class:`python:float`
           Hydraulic diameter [m].

       **x** : :class:`python:float`
           Axial coordinate.

       **roughness** : :func:`python:callable` or :data:`python:None`
           Function returning surface roughness at ``x`` [m], or ``None`` for smooth wall.



   :Returns:

       :class:`python:float`
           Darcy friction factor (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: f_darcy_laminar(ReDh, Dh, x)

   
   Laminar Darcy friction factor.


   :Parameters:

       **ReDh** : :class:`python:float`
           Hydraulic Reynolds number.

       **Dh** : :class:`python:float`
           Hydraulic diameter [m].

       **x** : :class:`python:float`
           Axial coordinate (unused, kept for API consistency).



   :Returns:

       :class:`python:float`
           Darcy friction factor (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: f_darcy_turbulent(ReDh, Dh, x, roughness)

   
   Turbulent Darcy friction factor using Putukhov or Colebrook–White.


   :Parameters:

       **ReDh** : :class:`python:float`
           Hydraulic Reynolds number.

       **Dh** : :class:`python:float`
           Hydraulic diameter [m].

       **x** : :class:`python:float`
           Axial coordinate.

       **roughness** : :func:`python:callable` or :data:`python:None`
           Function returning surface roughness at ``x`` [m], or ``None`` for smooth wall.



   :Returns:

       :class:`python:float`
           Darcy friction factor (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: h_coolant_colburn(k_cf, D_c, Cp_cr, mu_cf, mdot_c, A_c, phi_curv=1)

   
   Colburn correlation for turbulent coolant-side heat transfer.


   :Parameters:

       **k_cf** : :class:`python:float`
           Thermal conductivity of coolant film [W m⁻¹ K⁻¹].

       **D_c** : :class:`python:float`
           Coolant hydraulic diameter [m].

       **Cp_cr** : :class:`python:float`
           Mean specific heat (reference enthalpy) [J kg⁻¹ K⁻¹].

       **mu_cf** : :class:`python:float`
           Dynamic viscosity of coolant film [Pa s].

       **mdot_c** : :class:`python:float`
           Coolant mass-flow rate [kg s⁻¹].

       **A_c** : :class:`python:float`
           Coolant-flow area [m²].

       **phi_curv** : :class:`python:float`, :obj:`optional`
           Curvature-correction factor (dimensionless). Default is 1.



   :Returns:

       :class:`python:float`
           Coolant-side heat-transfer coefficient [W m⁻² K⁻¹].








   .. rubric:: Notes

   Subscripts: 
       Cf denotes film coolant film condition
       c denotes bulk coolant conditions
       r denotes reference enthalpy condition which are averaged between the free stream and wall metal conditions



   ..
       !! processed by numpydoc !!

.. py:function:: h_gas_bartz(D_t, mu_g, cp_g, Pr_g, p_c, c_star, A_t, A_x, sigma)

   
   Compute the Bartz hot-gas correlation without curvature correction.


   :Parameters:

       **D_t** : :class:`python:float`
           Throat diameter [m].

       **mu_g** : :class:`python:float`
           Gas viscosity [Pa s].

       **cp_g** : :class:`python:float`
           Gas specific heat capacity [J kg⁻¹ K⁻¹].

       **Pr_g** : :class:`python:float`
           Gas Prandtl number.

       **p_c** : :class:`python:float`
           Chamber pressure [Pa].

       **c_star** : :class:`python:float`
           Characteristic velocity [m s⁻¹].

       **A_t** : :class:`python:float`
           Throat area [m²].

       **A_x** : :class:`python:float`
           Local flow area [m²].

       **sigma** : :class:`python:float`
           Boundary-layer property-variation correction (dimensionless).



   :Returns:

       :class:`python:float`
           Heat-transfer coefficient [W m⁻² K⁻¹].











   ..
       !! processed by numpydoc !!

.. py:function:: h_gas_bartz_enthalpy_driven(k_gr, D_hyd, Cp_gr, mu_gr, mdot_g, A_chmb, T_g, T_gr)

   
   Compute the hot-gas-side heat-transfer coefficient (Bartz correlation, enthalpy-based).


   :Parameters:

       **k_gr** : :class:`python:float`
           Thermal conductivity at the reference enthalpy condition [W m⁻¹ K⁻¹].

       **D_hyd** : :class:`python:float`
           Hydraulic diameter of the gas passage [m].

       **Cp_gr** : :class:`python:float`
           Specific heat capacity at the reference enthalpy condition [J kg⁻¹ K⁻¹].

       **mu_gr** : :class:`python:float`
           Dynamic viscosity at the reference condition [Pa s].

       **mdot_g** : :class:`python:float`
           Total hot-gas mass flow rate [kg s⁻¹].

       **A_chmb** : :class:`python:float`
           Flow cross-sectional area in the chamber [m²].

       **T_g** : :class:`python:float`
           Free-stream gas temperature [K].

       **T_gr** : :class:`python:float`
           Reference-enthalpy temperature [K].



   :Returns:

       :class:`python:float`
           Hot-gas-side heat-transfer coefficient ``h_g`` [W m⁻² K⁻¹].








   .. rubric:: Notes

   This form of the Bartz correlation uses the ratio ``T_g/T_gr`` to
   account for variable-property effects between free-stream and reference conditions.

   subscript g denotes free stream gas properties
   subscript r denotes reference enthalpy conditions which are averaged between the free stream and wall metal conditions



   ..
       !! processed by numpydoc !!

.. py:function:: phi_curv(Re_c, D_c, R_curv)

   
   Curvature correction factor for the coolant-side heat-transfer coefficient.


   :Parameters:

       **Re_c** : :class:`python:float`
           Coolant Reynolds number.

       **D_c** : :class:`python:float`
           Hydraulic diameter [m].

       **R_curv** : :class:`python:float`
           Radius of curvature of the coolant passage [m].
           Use ``np.inf`` for straight sections.



   :Returns:

       :class:`python:float`
           Curvature factor ``φ`` (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: reynolds(rho, u, L, mu)

   
   Compute the Reynolds number.


   :Parameters:

       **rho** : :class:`python:float`
           Fluid density [kg m⁻³].

       **u** : :class:`python:float`
           Mean velocity [m s⁻¹].

       **L** : :class:`python:float`
           Characteristic length or hydraulic diameter [m].

       **mu** : :class:`python:float`
           Dynamic viscosity [Pa s].



   :Returns:

       :class:`python:float`
           Reynolds number (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: sigma(T_hw, T_c, gamma_g, M_g, omega)

   
   Dimensionless property-variation factor used in Bartz correlation.


   :Parameters:

       **T_hw** : :class:`python:float`
           Wall temperature [K].

       **T_c** : :class:`python:float`
           Core-flow (free-stream) temperature [K].

       **gamma_g** : :class:`python:float`
           Ratio of specific heats.

       **M_g** : :class:`python:float`
           Local Mach number.

       **omega** : :class:`python:float`
           Empirical property exponent (≈ 0.68 for diatomic gases).



   :Returns:

       :class:`python:float`
           ``σ`` correction factor (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: u_coolant(rho, mdot_c_single_channel, A_cool)

   
   Compute coolant velocity.


   :Parameters:

       **rho** : :class:`python:float`
           Coolant density [kg m⁻³].

       **mdot_c_single_channel** : :class:`python:float`
           Mass-flow rate through one channel [kg s⁻¹].

       **A_cool** : :class:`python:float`
           Coolant cross-sectional area [m²].



   :Returns:

       :class:`python:float`
           Coolant velocity [m s⁻¹].











   ..
       !! processed by numpydoc !!

.. py:data:: ReDh_laminar
   :value: 2300


.. py:data:: ReDh_turbulent
   :value: 3500


