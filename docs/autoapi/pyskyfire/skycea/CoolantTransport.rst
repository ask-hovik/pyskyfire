pyskyfire.skycea.CoolantTransport
=================================

.. py:class:: pyskyfire.skycea.CoolantTransport(fluid: pyskyfire.common.fluids.Fluid)

   
   CoolProp-backed transport/thermo properties for a single fluid.

   Thin wrapper around **CoolProp** that exposes convenient getters using a
   :class:`~pyskyfire.common.fluids.Fluid` to obtain the CoolProp fluid string.

   :Parameters:

       **fluid** : :obj:`pyskyfire.common.fluids.Fluid`
           Fluid descriptor; its ``.coolprop_string()`` is passed to CoolProp.

   :Attributes:

       **fluid** : :class:`python:str`
           CoolProp fluid identifier string.









   .. seealso::

       
       :obj:`CoolProp.CoolProp.PropsSI`
           Underlying property evaluator.
       :class:`~pyskyfire.common.fluids.Fluid`
           Source of the CoolProp-compatible fluid string.
       
       
   .. rubric:: Notes

   - All getters accept ``(T [K], p [Pa])`` in that order and return SI units.



   ..
       !! processed by numpydoc !!

   .. py:method:: get_Pr(T, p)


   .. py:method:: get_cp(T, p)


   .. py:method:: get_cv(T, p)


   .. py:method:: get_gamma(T, p)


   .. py:method:: get_k(T, p)


   .. py:method:: get_mu(T, p)


   .. py:method:: get_rho(T, p)

