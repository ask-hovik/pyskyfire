pyskyfire.skycea.TransportProperties
====================================

.. py:class:: pyskyfire.skycea.TransportProperties(Pr, mu, k, cp=None, rho=None, gamma_coolant=None)

   
   Container for transport-property models (constants or callables).

   This class stores transport properties that can be either constant scalars
   or callables of temperature and pressure in the order ``(T [K], p [Pa])``.
   It is suitable for **walls**, **coolants**, or generic working fluids. When
   ``gamma_coolant`` is supplied, the instance is treated as a *compressible
   coolant*.

   :Parameters:

       **Pr** : :class:`python:float` or :func:`python:callable`
           Prandtl number [-] or a function ``Pr(T, p) -> float``.

       **mu** : :class:`python:float` or :func:`python:callable`
           Dynamic (absolute) viscosity [Pa·s] or a function ``mu(T, p) -> float``.

       **k** : :class:`python:float` or :func:`python:callable`
           Thermal conductivity [W/(m·K)] or a function ``k(T, p) -> float``.

       **cp** : :class:`python:float` or :func:`python:callable`, :obj:`optional`
           Isobaric specific heat [J/(kg·K)] or ``cp(T, p) -> float``. Required when
           this object represents a coolant.

       **rho** : :class:`python:float` or :func:`python:callable`, :obj:`optional`
           Density [kg/m³] or ``rho(T, p) -> float``. Required when this object
           represents a coolant.

       **gamma_coolant** : :class:`python:float` or :func:`python:callable`, :obj:`optional`
           Ratio of specific heats ``γ = cp/cv`` [-] or ``gamma(T, p) -> float``.
           If provided, the object is assumed to represent a **compressible coolant**.

   :Attributes:

       **compressible_coolant** : :ref:`bool <python:bltin-boolean-values>`
           Whether the instance represents a compressible coolant (``gamma_coolant``
           was provided).

       **_Pr, _mu, _k, _cp, _rho, _gamma_coolant** : :obj:`Any`
           Backing values/callables as provided.









   .. seealso::

       
       :class:`~pyskyfire.common.fluids.Fluid`
           Fluid identifier utilities used elsewhere in the package.
       
       
   .. rubric:: Notes

   - All callables must accept ``(T [K], p [Pa])`` in that order.
   - For coolants, both ``cp`` and ``rho`` should be provided to allow
   consistent thermo-fluid calculations.



   ..
       !! processed by numpydoc !!

   .. py:method:: Pr(T, p)

      
      Prandtl number.


      :Parameters:

          **T** : :class:`python:float`
              Temperature [K].

          **p** : :class:`python:float`
              Pressure [Pa].



      :Returns:

          :class:`python:float`
              Prandtl number [-].











      ..
          !! processed by numpydoc !!


   .. py:method:: cp(T, p)

      
      Isobaric specific heat capacity.


      :Parameters:

          **T** : :class:`python:float`
              Temperature [K].

          **p** : :class:`python:float`
              Pressure [Pa].



      :Returns:

          :class:`python:float`
              ``c_p`` [J/(kg·K)].




      :Raises:

          :obj:`ValueError`
              If ``cp`` was not provided at construction.







      ..
          !! processed by numpydoc !!


   .. py:method:: gamma_coolant(T, p)

      
      Ratio of specific heats for a compressible coolant.


      :Parameters:

          **T** : :class:`python:float`
              Temperature [K].

          **p** : :class:`python:float`
              Pressure [Pa].



      :Returns:

          :class:`python:float`
              ``γ = c_p / c_v`` [-].




      :Raises:

          :obj:`ValueError`
              If ``gamma_coolant`` was not provided at construction.







      ..
          !! processed by numpydoc !!


   .. py:method:: k(T, p)

      
      Thermal conductivity.


      :Parameters:

          **T** : :class:`python:float`
              Temperature [K].

          **p** : :class:`python:float`
              Pressure [Pa].



      :Returns:

          :class:`python:float`
              Thermal conductivity [W/(m·K)].











      ..
          !! processed by numpydoc !!


   .. py:method:: mu(T, p)

      
      Dynamic (absolute) viscosity.


      :Parameters:

          **T** : :class:`python:float`
              Temperature [K].

          **p** : :class:`python:float`
              Pressure [Pa].



      :Returns:

          :class:`python:float`
              Viscosity [Pa·s].











      ..
          !! processed by numpydoc !!


   .. py:method:: rho(T, p)

      
      Density.


      :Parameters:

          **T** : :class:`python:float`
              Temperature [K].

          **p** : :class:`python:float`
              Pressure [Pa].



      :Returns:

          :class:`python:float`
              Density [kg/m³].




      :Raises:

          :obj:`ValueError`
              If ``rho`` was not provided at construction.







      ..
          !! processed by numpydoc !!

