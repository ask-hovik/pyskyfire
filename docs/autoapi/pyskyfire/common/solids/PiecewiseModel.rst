pyskyfire.common.solids.PiecewiseModel
======================================

.. py:class:: pyskyfire.common.solids.PiecewiseModel

   Bases: :py:obj:`PropertyModel`

   .. autoapi-inheritance-diagram:: pyskyfire.common.solids.PiecewiseModel
      :parts: 1
      :private-bases:


   
   Combine multiple sub-models across temperature segments.

   Each segment is a tuple ``(T_lo, T_hi, model)`` where ``model`` is a callable
   that returns property values. The segments must have strictly increasing
   temperature bounds.

   :Parameters:

       **segments** : :class:`python:list`\[:class:`python:tuple`\[:class:`python:float`, :class:`python:float`, :obj:`PropertyModel`]]
           Segment list defining temperature intervals and their models.

       **range_policy** : :class:`python:str`, :obj:`optional`
           Range-handling mode: ``"warn_clip"`` (default), ``"clip"`` or ``"error"``.

       **blend** : :class:`python:float`, :obj:`optional`
           Retained for API compatibility (unused in the new rules).











   .. rubric:: Notes

   The evaluation rules are:

   - **Gaps:** linearly interpolate between the nearest segment endpoints.
   - **Overlaps:** average results from all covering segments.
   - **Out-of-range:** clip to the nearest valid edge (warns if ``"warn_clip"``).



   ..
       !! processed by numpydoc !!

   .. py:method:: __call__(T: ArrayLike) -> numpy.ndarray


   .. py:method:: __post_init__()


   .. py:method:: _covering_indices(T: float) -> numpy.ndarray

      
      Return indices of all segments that cover T (inclusive bounds).
















      ..
          !! processed by numpydoc !!


   .. py:method:: _edge_value(idx: int, side: str) -> float

      
      Evaluate model at the exact segment edge (lo/hi) and cache it.
















      ..
          !! processed by numpydoc !!


   .. py:method:: _find_neighbors(T: float) -> Tuple[Optional[int], Optional[int]]

      
      Return indices (i_left, i_right) of nearest segments such that:
        end[i_left] <= T <= start[i_right], with end[i_left] < start[i_right].
      If not found on one side, that index is None.
















      ..
          !! processed by numpydoc !!


   .. py:method:: _handle_oob(T: float) -> float

      
      Out-of-bounds policy: clip to nearest edge, warn as configured.
















      ..
          !! processed by numpydoc !!

