pyskyfire.common.TabulatedModel
===============================

.. py:class:: pyskyfire.common.TabulatedModel

   Bases: :py:obj:`PropertyModel`

   .. autoapi-inheritance-diagram:: pyskyfire.common.TabulatedModel
      :parts: 1
      :private-bases:


   
   Property model based on tabulated data.

   Performs 1D linear interpolation between discrete (T, Y) points.

   :Parameters:

       **Ts** : :term:`numpy:array_like`
           Monotonic array of temperatures [K].

       **Ys** : :term:`numpy:array_like`
           Property values corresponding to :attr:`Ts`.

       **enforce_bounds** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           If ``True``, out-of-range temperatures raise :class:`ValueError`.
           Default is ``False``.







   :Raises:

       :obj:`ValueError`
           If ``Ts`` and ``Ys`` have mismatched shapes or are not strictly increasing.







   ..
       !! processed by numpydoc !!

   .. py:method:: __call__(T: ArrayLike) -> numpy.ndarray

      
      Interpolate the property at temperature ``T``.


      :Parameters:

          **T** : :obj:`ArrayLike`
              Temperatures [K].



      :Returns:

          :obj:`np.ndarray <numpy.ndarray>`
              Interpolated property values.











      ..
          !! processed by numpydoc !!


   .. py:method:: __post_init__()

