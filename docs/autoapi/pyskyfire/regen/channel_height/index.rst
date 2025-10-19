pyskyfire.regen.channel_height
==============================

.. py:module:: pyskyfire.regen.channel_height


Functions
---------

.. autoapisummary::

   pyskyfire.regen.channel_height.make_channel_height_fn


Module Contents
---------------

.. py:function:: make_channel_height_fn(contour: pyskyfire.regen.thrust_chamber.Contour, region_fractions: Sequence[float], flat_heights: Sequence[float], pinch_factors: Sequence[float], transition_widths: Union[float, Sequence[float]], logistic_k: float = 10.0) -> Callable[[float], float]

   
   Construct a smooth channel-height profile along a thrust chamber contour.

   Builds a continuous, piecewise-logistic blend of per-region channel heights
   across the chamber and nozzle. Each region is defined by a normalized axial
   coordinate fraction relative to the throat, and can optionally "pinch"
   channel height with radius.

   The returned function evaluates the local channel height at any axial
   coordinate ``x`` along the contour.

   :Parameters:

       **contour** : :obj:`Contour`
           Thrust-chamber contour defining the wall coordinates ``x`` and ``r(x)``.

       **region_fractions** : :obj:`Sequence`\[:class:`python:float`]
           Normalized axial coordinates marking region boundaries:
           
           - ``-1.0`` → chamber inlet
           - ``0.0`` → throat
           - ``+1.0`` → nozzle exit

       **flat_heights** : :obj:`Sequence`\[:class:`python:float`]
           Base (unpinched) channel heights for each region [m].

       **pinch_factors** : :obj:`Sequence`\[:class:`python:float`]
           Fraction in ``[0, 1]`` describing how strongly each region’s height scales
           with local radius. ``0`` means constant height, ``1`` means fully proportional
           to radius.

       **transition_widths** : :class:`python:float` or :obj:`Sequence`\[:class:`python:float`]
           Axial transition width(s) controlling how quickly heights blend between
           adjacent regions. If a single scalar is given, it is used for all boundaries.

       **logistic_k** : :class:`python:float`, :obj:`optional`
           Steepness of the logistic blending function. Higher values yield sharper
           transitions. Default is ``10.0``.



   :Returns:

       :obj:`Callable`\[[:class:`python:float`], :class:`python:float`]
           A callable function ``channel_height(x)`` that returns the local channel
           height [m] at axial coordinate ``x``.




   :Raises:

       :obj:`ValueError`
           If the number of transition widths does not match ``len(flat_heights) - 1``.




   .. rubric:: Notes

   The normalized coordinate mapping is:

   - Negative region fractions → chamber side
   - Positive region fractions → nozzle side

   Transitions are blended smoothly using a logistic weighting function centered
   at each region boundary. Channel heights can pinch with radius according to
   ``pinch_factors``.


   .. rubric:: Examples

   >>> channel_height_fn = psf.regen.make_channel_height_fn(
               contour=contour, 
               region_fractions=[-1.0, 0.25, 1.0], 
               flat_heights= [0.0032, 0.00134], 
               pinch_factors= [0.6, -5.0], 
               transition_widths=[0.1]
   >>> fn(0.12)  # Evaluate channel height at some axial x
   0.0034

   ..
       !! processed by numpydoc !!

