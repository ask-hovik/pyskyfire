pyskyfire.skycea.coolant_transport.TransportProperties
======================================================

.. py:class:: pyskyfire.skycea.coolant_transport.TransportProperties(Pr, mu, k, cp=None, rho=None, gamma_coolant=None)

   .. py:method:: Pr(T, p)

      
      Prandtl number.

      Args:
          T (float): Temperature (K)
          p (float): Pressure (Pa)

      Returns:
          float: Prandtl number















      ..
          !! processed by numpydoc !!


   .. py:method:: cp(T, p)

      
      Isobaric specific heat capacity (J/kg/K)

      Args:
          T (float): Temperature (K)
          p (float): Pressure (Pa)

      Returns:
          float: Isobaric specific heat capacity (J/kg/K)















      ..
          !! processed by numpydoc !!


   .. py:method:: gamma_coolant(T, p)

      
      Ratio of specific heat capacities for a compressible coolant.

      Args:
          T (float): Temperature (K)
          p (float): Pressure (Pa)

      Returns:
          float: Ratio of specific heat capacities (cp/cv).















      ..
          !! processed by numpydoc !!


   .. py:method:: k(T, p)

      
      Thermal conductivity (W/m/K)

      Args:
          T (float): Temperature (K)
          p (float): Pressure (Pa)

      Returns:
          float: Thermal conductivity (W/m/K)















      ..
          !! processed by numpydoc !!


   .. py:method:: mu(T, p)

      
      Absolute viscosity (Pa s)

      Args:
          T (float): Temperature (K)
          p (float): Pressure (Pa)

      Returns:
          float: Absolute viscosity (Pa s)















      ..
          !! processed by numpydoc !!


   .. py:method:: rho(T, p)

      
      Density (kg/m^3)
      Args:
          T (float): Temperature (K)
          p (float): Pressure (Pa)
      Returns:
          float: Density (kg/m^3)
















      ..
          !! processed by numpydoc !!

