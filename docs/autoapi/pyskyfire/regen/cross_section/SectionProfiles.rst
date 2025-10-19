pyskyfire.regen.cross_section.SectionProfiles
=============================================

.. py:class:: pyskyfire.regen.cross_section.SectionProfiles

   
   All geometry profiles along the cooling-channel centerline.

   Each array must have consistent length ``N`` (number of stations).


   :Attributes:

       **h** : :obj:`np.ndarray <numpy.ndarray>`
           Channel height array [m, shape (N,)].

       **theta** : :obj:`np.ndarray <numpy.ndarray>`
           Included wedge angle at each station [rad, shape (N,)].

       **t_wall** : :obj:`np.ndarray <numpy.ndarray>`
           Wall thickness between coolant and hot-gas side [m, shape (N,)].

       **centerline** : :obj:`np.ndarray <numpy.ndarray>`
           Centerline coordinates ``(x, r, θ)`` or equivalent system [shape (N, 3)].

       **local_coords** : :obj:`np.ndarray <numpy.ndarray>`
           Local orthonormal frames ``(t, n, b)`` at each station [shape (N, 3, 3)].

       **blockage_ratio** : :obj:`np.ndarray <numpy.ndarray>`
           Fraction [0–1] describing coolant blockage (1 = fully blocked).






   :Raises:

       :obj:`ValueError`
           If any array length or shape is inconsistent.







   ..
       !! processed by numpydoc !!

   .. py:method:: __post_init__()

