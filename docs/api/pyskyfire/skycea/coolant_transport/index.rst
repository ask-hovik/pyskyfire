pyskyfire.skycea.coolant_transport
==================================

.. py:module:: pyskyfire.skycea.coolant_transport




Module Contents
---------------

.. py:class:: TransportProperties(Pr, mu, k, cp=None, rho=None, gamma_coolant=None)

   
   Container for specifying your transport properties. Each input can either be a function of temperature (K) and pressure (Pa) in that order, e.g. mu(T, p). Otherwise they can be constant floats.

   :param Pr: Prandtl number.
   :type Pr: float or callable
   :param mu: Absolute viscosity (Pa s).
   :type mu: float or callable
   :param k: Thermal conductivity (W/m/K).
   :type k: float or callable
   :param cp: Isobaric specific heat capacity (J/kg/K) - only required for coolants.
   :type cp: float or callable, optional
   :param rho: Density (kg/m^3) - only required for coolants.
   :type rho: float or callable, optional
   :param gamma_coolant: Ratio of specific heats (cp/cv) for a compressible coolant. If this is submitted, it is assumed that this object represents a compressible coolant.
   :type gamma_coolant: float or callable, optional

   .. attribute:: compressible_coolant

      Whether or not this TransportProperties object represents a compressible coolant.

      :type: bool


   .. py:attribute:: type


   .. py:method:: Pr(T, p)

      Prandtl number.

      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Prandtl number
      :rtype: float



   .. py:method:: mu(T, p)

      Absolute viscosity (Pa s)

      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Absolute viscosity (Pa s)
      :rtype: float



   .. py:method:: k(T, p)

      Thermal conductivity (W/m/K)

      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Thermal conductivity (W/m/K)
      :rtype: float



   .. py:method:: rho(T, p)

      Density (kg/m^3)
      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Density (kg/m^3)
      :rtype: float



   .. py:method:: cp(T, p)

      Isobaric specific heat capacity (J/kg/K)

      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Isobaric specific heat capacity (J/kg/K)
      :rtype: float



   .. py:method:: gamma_coolant(T, p)

      Ratio of specific heat capacities for a compressible coolant.

      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Ratio of specific heat capacities (cp/cv).
      :rtype: float



.. py:class:: CoolantTransport(fluid)

   .. py:attribute:: fluid


   .. py:method:: get_Pr(T, p)


   .. py:method:: get_mu(T, p)


   .. py:method:: get_k(T, p)


   .. py:method:: get_cp(T, p)


   .. py:method:: get_rho(T, p)


   .. py:method:: get_cv(T, p)


   .. py:method:: get_gamma(T, p)


