pyskyfire.common.solids.Log10PolynomialModel
============================================

.. py:class:: pyskyfire.common.solids.Log10PolynomialModel

   Bases: :py:obj:`PropertyModel`

   .. autoapi-inheritance-diagram:: pyskyfire.common.solids.Log10PolynomialModel
      :parts: 1
      :private-bases:


   
   Log-polynomial model (NIST cryogenic form).

   Implements :math:`y = 10^{P(\log_{10} T)}`,
   where :math:`P(x) = \sum_i a_i x^i`.

   :Parameters:

       **coeffs** : :obj:`Iterable`\[:class:`python:float`]
           Polynomial coefficients ``[a0, a1, a2, ...]``.

       **Tmin, Tmax** : :class:`python:float`, :obj:`optional`
           Optional temperature bounds for validity.

       **enforce_bounds** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           If ``True``, out-of-range temperatures raise :class:`ValueError`.
           Default is ``False``.











   .. rubric:: Notes

   Useful for fitting data that spans multiple orders of magnitude in
   temperature or property values.



   ..
       !! processed by numpydoc !!

   .. py:method:: __call__(T: ArrayLike) -> numpy.ndarray

      
      Evaluate the log-polynomial model.


      :Parameters:

          **T** : :obj:`ArrayLike`
              Temperature array (must be positive).



      :Returns:

          :obj:`np.ndarray <numpy.ndarray>`
              Evaluated property values.




      :Raises:

          :obj:`ValueError`
              If any ``T <= 0`` or bounds are violated.







      ..
          !! processed by numpydoc !!

