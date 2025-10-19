pyskyfire.common.solids
=======================

.. py:module:: pyskyfire.common.solids






Module Contents
---------------

.. py:data:: ArrayLike

.. py:class:: PropertyModel

.. py:class:: ConstantModel

   Bases: :py:obj:`PropertyModel`


   .. py:attribute:: value
      :type:  float


.. py:class:: PolynomialModel

   Bases: :py:obj:`PropertyModel`


   y = c0 + c1*T + c2*T^2 + ...  (no bounds unless provided)


   .. py:attribute:: coeffs
      :type:  Iterable[float]


   .. py:attribute:: Tmin
      :type:  Optional[float]
      :value: None



   .. py:attribute:: Tmax
      :type:  Optional[float]
      :value: None



   .. py:attribute:: enforce_bounds
      :type:  bool
      :value: False



.. py:class:: Log10PolynomialModel

   Bases: :py:obj:`PropertyModel`


   y = 10 ** P(log10(T)), P(x) = sum a_i x^i
   Bounds optional; leave unenforced and let Piecewise handle range policy.


   .. py:attribute:: coeffs
      :type:  Iterable[float]


   .. py:attribute:: Tmin
      :type:  Optional[float]
      :value: None



   .. py:attribute:: Tmax
      :type:  Optional[float]
      :value: None



   .. py:attribute:: enforce_bounds
      :type:  bool
      :value: False



.. py:class:: TabulatedModel

   Bases: :py:obj:`PropertyModel`


   .. py:attribute:: Ts
      :type:  numpy.ndarray


   .. py:attribute:: Ys
      :type:  numpy.ndarray


   .. py:attribute:: enforce_bounds
      :type:  bool
      :value: False



.. py:class:: SumOfGaussiansModel

   Bases: :py:obj:`PropertyModel`


   .. py:attribute:: params
      :type:  List[Tuple[float, float, float]]


.. py:class:: PiecewiseModel

   Bases: :py:obj:`PropertyModel`


   Piecewise wrapper over sub-models defined on [T_lo, T_hi] segments.

   New rules:
     • GAPS: if no segment covers T but T lies between two segments, linearly
       interpolate between the *endpoint values* of those segments.
     • OVERLAPS: if multiple segments cover T, evaluate all and return the average.
     • OUT-OF-RANGE: clip to nearest edge; if 'warn_clip' in range_policy, warn.

   .. rubric:: Notes

   • Sub-models must be callable: y = model(T).
   • Segments are tuples: (T_lo, T_hi, model) with T_lo < T_hi (strict).


   .. py:attribute:: segments
      :type:  List[Tuple[float, float, PropertyModel]]


   .. py:attribute:: range_policy
      :type:  str
      :value: 'warn_clip'



   .. py:attribute:: blend
      :type:  float
      :value: 0.0



.. py:class:: Material

   .. py:attribute:: name
      :type:  str


   .. py:attribute:: k
      :type:  Optional[PropertyModel]
      :value: None



   .. py:attribute:: E
      :type:  Optional[PropertyModel]
      :value: None



   .. py:attribute:: alpha
      :type:  Optional[PropertyModel]
      :value: None



   .. py:attribute:: nu
      :type:  Optional[PropertyModel]
      :value: None



   .. py:attribute:: rho
      :type:  Optional[PropertyModel]
      :value: None



   .. py:method:: get_k(T)


   .. py:method:: get_E(T)


   .. py:method:: get_alpha(T)


   .. py:method:: get_nu(T)


   .. py:method:: get_rho(T)


.. py:data:: log10poly_304_cryo

.. py:data:: mills_304_highT

.. py:data:: k_304_piecewise

.. py:data:: StainlessSteel304

.. py:data:: k718_cryo

.. py:data:: T_F

.. py:data:: k_BTUin

.. py:data:: T718_hi

.. py:data:: k718_hi

.. py:data:: k718_table

.. py:data:: k718

.. py:data:: Inconel718

.. py:data:: T625_C

.. py:data:: k625_W

.. py:data:: T625

.. py:data:: k625

.. py:data:: Inconel625

.. py:data:: k42_mono

.. py:data:: GRCop42

.. py:data:: ZirconiumOxide

.. py:data:: TEOS

