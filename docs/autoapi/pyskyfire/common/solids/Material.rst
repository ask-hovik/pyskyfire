pyskyfire.common.solids.Material
================================

.. py:class:: pyskyfire.common.solids.Material

   
   Container for material name and temperature-dependent properties.


   :Parameters:

       **name** : :class:`python:str`
           Material name.

       **k** : :obj:`PropertyModel`, :obj:`optional`
           Thermal conductivity model [W/m·K].

       **E** : :obj:`PropertyModel`, :obj:`optional`
           Young’s modulus model [Pa].

       **alpha** : :obj:`PropertyModel`, :obj:`optional`
           Thermal expansion coefficient model [1/K].

       **nu** : :obj:`PropertyModel`, :obj:`optional`
           Poisson’s ratio model.

       **rho** : :obj:`PropertyModel`, :obj:`optional`
           Density model [kg/m³].


   .. rubric:: Methods



   ================  ==========
       **get_k(T)**  Evaluate thermal conductivity.  
       **get_E(T)**  Evaluate Young’s modulus.  
   **get_alpha(T)**  Evaluate thermal expansion coefficient.  
      **get_nu(T)**  Evaluate Poisson’s ratio.  
     **get_rho(T)**  Evaluate density.  
   ================  ==========









   .. rubric:: Notes

   Each getter raises :class:`AttributeError` if the property is not defined.



   ..
       !! processed by numpydoc !!

   .. py:method:: get_E(T: ArrayLike) -> numpy.ndarray

      
      Evaluate Young’s modulus.
















      ..
          !! processed by numpydoc !!


   .. py:method:: get_alpha(T: ArrayLike) -> numpy.ndarray

      
      Evaluate coefficient of thermal expansion.
















      ..
          !! processed by numpydoc !!


   .. py:method:: get_k(T: ArrayLike) -> numpy.ndarray

      
      Evaluate thermal conductivity.


      :Parameters:

          **T** : :obj:`ArrayLike`
              Temperature array [K].



      :Returns:

          :obj:`np.ndarray <numpy.ndarray>`
              Thermal conductivity [W/m·K].




      :Raises:

          :obj:`AttributeError`
              If no model is defined for :attr:`k`.







      ..
          !! processed by numpydoc !!


   .. py:method:: get_nu(T: ArrayLike) -> numpy.ndarray

      
      Evaluate Poisson’s ratio.
















      ..
          !! processed by numpydoc !!


   .. py:method:: get_rho(T: ArrayLike) -> numpy.ndarray

      
      Evaluate density.
















      ..
          !! processed by numpydoc !!

