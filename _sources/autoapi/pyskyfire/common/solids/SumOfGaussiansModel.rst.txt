pyskyfire.common.solids.SumOfGaussiansModel
===========================================

.. py:class:: pyskyfire.common.solids.SumOfGaussiansModel

   Bases: :py:obj:`PropertyModel`

   .. autoapi-inheritance-diagram:: pyskyfire.common.solids.SumOfGaussiansModel
      :parts: 1
      :private-bases:


   
   Empirical model as a sum of Gaussian peaks.

   Implements :math:`y(T) = \sum_i b_i \exp(-((T - c_i)/d_i)^2)`.

   :Parameters:

       **params** : :class:`python:list`\[:class:`python:tuple`\[:class:`python:float`, :class:`python:float`, :class:`python:float`]]
           List of triplets ``(b, c, d)`` representing amplitude, center, and width.














   ..
       !! processed by numpydoc !!

   .. py:method:: __call__(T)

      
      Evaluate the Gaussian sum model.


      :Parameters:

          **T** : :obj:`ArrayLike`
              Temperature array.



      :Returns:

          :obj:`np.ndarray <numpy.ndarray>`
              Property values.











      ..
          !! processed by numpydoc !!

