pyskyfire.common.PolynomialModel
================================

.. py:class:: pyskyfire.common.PolynomialModel

   Bases: :py:obj:`PropertyModel`

   .. autoapi-inheritance-diagram:: pyskyfire.common.PolynomialModel
      :parts: 1
      :private-bases:


   
   Ordinary polynomial property model.

   Implements :math:`y = c_0 + c_1 T + c_2 T^2 + ...`.

   :Parameters:

       **coeffs** : :obj:`Iterable`\[:class:`python:float`]
           Polynomial coefficients ``[c0, c1, c2, ...]``.

       **Tmin, Tmax** : :class:`python:float`, :obj:`optional`
           Optional temperature bounds for validity.

       **enforce_bounds** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           If ``True``, out-of-range temperatures raise :class:`ValueError`.
           Default is ``False``.











   .. rubric:: Notes

   When used inside :class:`PiecewiseModel`, range policy is typically handled
   externally.



   ..
       !! processed by numpydoc !!

   .. py:method:: __call__(T: ArrayLike) -> numpy.ndarray

      
      Evaluate the polynomial.


      :Parameters:

          **T** : :obj:`ArrayLike`
              Temperature array.



      :Returns:

          :obj:`np.ndarray <numpy.ndarray>`
              Evaluated property values.











      ..
          !! processed by numpydoc !!

