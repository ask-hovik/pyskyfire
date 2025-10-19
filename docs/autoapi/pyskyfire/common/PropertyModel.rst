pyskyfire.common.PropertyModel
==============================

.. py:class:: pyskyfire.common.PropertyModel

   
   Abstract base class for temperature-dependent property models.

   All subclasses must implement ``__call__(T) -> np.ndarray``.











   .. seealso::

       
       :obj:`ConstantModel`
           Fixed value, independent of temperature.
       :obj:`PolynomialModel`
           Ordinary polynomial :math:`y = \sum_i c_i T^i`.
       :obj:`Log10PolynomialModel`
           Log-polynomial of :math:`T` as used in cryogenic data.
       :obj:`TabulatedModel`
           Linear interpolation from discrete data points.
       :obj:`SumOfGaussiansModel`
           Smooth empirical fit as a sum of Gaussian peaks.
       :obj:`PiecewiseModel`
           Combines multiple sub-models over temperature segments.
       
       



   ..
       !! processed by numpydoc !!

   .. py:method:: __call__(T: ArrayLike) -> numpy.ndarray
      :abstractmethod:


