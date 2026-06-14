pyskyfire.common.solids.ConstantModel
=====================================

.. py:class:: pyskyfire.common.solids.ConstantModel

   Bases: :py:obj:`PropertyModel`

   .. autoapi-inheritance-diagram:: pyskyfire.common.solids.ConstantModel
      :parts: 1
      :private-bases:


   
   Constant property model returning a fixed value.


   :Parameters:

       **value** : :class:`python:float`
           Constant property value.













   .. rubric:: Examples

   >>> ConstantModel(10.0)([1, 2, 3])
   array([10., 10., 10.])

   ..
       !! processed by numpydoc !!

   .. py:method:: __call__(T: ArrayLike) -> numpy.ndarray

      
      Return a constant array equal to :attr:`value`.


      :Parameters:

          **T** : :obj:`ArrayLike`
              Temperature array (ignored).



      :Returns:

          :obj:`np.ndarray <numpy.ndarray>`
              Array of same shape as ``T`` filled with :attr:`value`.











      ..
          !! processed by numpydoc !!

