pyskyfire.regen.contour
=======================

.. py:module:: pyskyfire.regen.contour

.. autoapi-nested-parse::

   contour.py

   Functinos and methods that together create a thrust chamber contour (combustion chamber + nozzle).



   This code takes in two json files theta_n.json and theta_e.json. These two files define the angles of the nozzle
   near the throat and the exit of the nozzle. The data comes for discrete length fractions, with discrete
   data points, so the code interpolates in the length fraction, area ratio, theta surface to yield the answer. 


   References:
       - [1] - The Thrust Optimised Parabolic nozzle, AspireSpace, 
         http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf

   ..
       !! processed by numpydoc !!


Functions
---------

.. autoapisummary::

   pyskyfire.regen.contour.compute_chamber_volume
   pyskyfire.regen.contour.compute_cutoff_length
   pyskyfire.regen.contour.get_contour
   pyskyfire.regen.contour.get_contour_internal
   pyskyfire.regen.contour.get_theta_e_n
   pyskyfire.regen.contour.integrate_area


Module Contents
---------------

.. py:function:: compute_chamber_volume(xs, rs)

   
   Compute chamber volume by revolving the contour about the x-axis.

   Integrates :math:`\pi r(x)^2` from the left boundary up to the throat
   (the minimum radius).

   :Parameters:

       **xs** : :term:`numpy:array_like`
           Monotone-increasing axial coordinates [m].

       **rs** : :term:`numpy:array_like`
           Radii [m] corresponding to ``xs``.



   :Returns:

       :class:`python:float`
           Volume [m³].








   .. rubric:: Notes

   Uses :func:`numpy.trapezoid` for the integral. If the throat occurs at the
   first point, returns ``0.0``.



   ..
       !! processed by numpydoc !!

.. py:function:: compute_cutoff_length(V_goal, xs_chamber, ys_chamber)

   
   Find ``L_c`` such that the chamber volume from 0 to ``-L_c`` equals ``V_goal``.


   :Parameters:

       **V_goal** : :class:`python:float`
           Target volume [m³].

       **xs_chamber** : :term:`numpy:array_like`, :obj:`shape` (N,)
           Axial coordinates; must include non-positive values.

       **ys_chamber** : :term:`numpy:array_like`, :obj:`shape` (N,)
           Radii corresponding to ``xs_chamber``.



   :Returns:

       :class:`python:float`
           ``L_c = |x_cutoff|`` where the cumulative volume first reaches ``V_goal``.




   :Raises:

       :obj:`ValueError`
           If ``V_goal < 0``, no non-positive ``x`` are available, or the target
           volume exceeds the total available volume.







   ..
       !! processed by numpydoc !!

.. py:function:: get_contour(r_t, area_ratio, r_c=None, L_c=None, V_c=None, eps_c=None, AR_c=None, theta_conv=45, theta_div=15, nozzle='rao', R_1f=1.5, R_2f=0.5, R_3f=0.382, length_fraction=0.8, angle_input='degrees', export_tikz=False)

   
   High-level API to generate a nozzle contour using four input modes.

   Exactly one of the following input combinations must be provided:

   1. **Direct**: ``r_c`` and ``L_c``.
   2. **Volume+eps**: ``V_c`` and ``eps_c`` → solve for ``r_c`` then for ``L_c``.
   3. **Volume+AR**: ``V_c`` and ``AR_c`` → solve for ``r_c`` then for ``L_c``.
   4. **Volume+r_c**: ``V_c`` and ``r_c`` → solve for ``L_c``.

   :Parameters:

       **r_t** : :class:`python:float`
           Throat radius [m].

       **area_ratio** : :class:`python:float`
           Nozzle area ratio :math:`\epsilon = A_e/A_t`.

       **r_c, L_c** : :class:`python:float`, :obj:`optional`
           Chamber radius/length [m] (Direct mode).

       **V_c** : :class:`python:float`, :obj:`optional`
           Chamber volume target [m³].

       **eps_c** : :class:`python:float`, :obj:`optional`
           Chamber area ratio (used with ``V_c`` to infer ``r_c``).

       **AR_c** : :class:`python:float`, :obj:`optional`
           Chamber aspect ratio target (used with ``V_c`` to infer ``r_c``).

       **theta_conv, theta_div** : :class:`python:float`, :obj:`optional`
           Convergence/divergence angles (degrees if ``angle_input='degrees'``).

       **nozzle** : {'rao', 'conical'}, :obj:`optional`
           Nozzle geometry. Default ``'rao'``.

       **R_1f, R_2f, R_3f** : :class:`python:float`, :obj:`optional`
           Throat/chamber fillet/exit curvature scale factors (× ``r_t``).
           ``R_2f=0`` or ``None`` gives a hard corner.

       **length_fraction** : :class:`python:float`, :obj:`optional`
           Rao length fraction in ``[0.60, 1.00]`` (Default ``0.8``).

       **angle_input** : {'degrees', 'radians'}, :obj:`optional`
           Units for ``theta_conv``/``theta_div``. Default ``'degrees'``.

       **export_tikz** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           Placeholder flag for downstream export.



   :Returns:

       :class:`python:tuple`\[:obj:`np.ndarray <numpy.ndarray>`, :obj:`np.ndarray <numpy.ndarray>`]
           ``(xs, ys)`` arrays defining the contour.




   :Raises:

       :obj:`ValueError`
           If the input combination is invalid, or the internal minimization/root
           solve fails to produce a consistent chamber length.



   .. seealso::

       
       :obj:`get_contour_internal`
           Low-level builder used by all modes.
       :obj:`compute_chamber_volume`
           Volume integral used by modes that solve for ``L_c``.
       
       



   ..
       !! processed by numpydoc !!

.. py:function:: get_contour_internal(r_c, r_t, area_ratio, L_c, theta_conv, theta_div, nozzle, R_1f, R_2f, R_3f, length_fraction, export_tikz)

   
   Generate the full nozzle contour coordinates from geometric inputs.

   Builds chamber and nozzle segments (with optional chamber fillet), throat
   arcs (entrant/exit), and either a conical or Rao-type nozzle, then joins
   and post-processes them to ensure strictly increasing ``x`` values.

   :Parameters:

       **r_c** : :class:`python:float`
           Chamber radius [m].

       **r_t** : :class:`python:float`
           Throat radius [m].

       **area_ratio** : :class:`python:float`
           Nozzle area ratio :math:`\epsilon = A_e/A_t`.

       **L_c** : :class:`python:float`
           Chamber length [m].

       **theta_conv** : :class:`python:float`
           Convergence angle (radians).

       **theta_div** : :class:`python:float`
           Divergence angle (radians). Used for conical nozzle and fallback.

       **nozzle** : :class:`python:str`
           ``"rao"`` or ``"conical"``.

       **R_1f** : :class:`python:float`
           Scale factor for throat entrant curvature radius (multiplies ``r_t``).

       **R_2f** : :class:`python:float` or :data:`python:None`
           Scale factor for chamber fillet radius (``0``/``None`` → hard corner).

       **R_3f** : :class:`python:float`
           Scale factor for throat exit curvature radius (multiplies ``r_t``).

       **length_fraction** : :class:`python:float`
           Normalized Rao length fraction in ``[0.60, 1.00]``.

       **export_tikz** : :ref:`bool <python:bltin-boolean-values>`
           (Unused here) placeholder for downstream export.



   :Returns:

       :class:`python:tuple`\[:obj:`np.ndarray <numpy.ndarray>`, :obj:`np.ndarray <numpy.ndarray>`]
           ``(xs, ys)`` where ``xs`` are axial coordinates [m] and ``ys`` are
           radii [m], strictly increasing in ``x``.




   :Raises:

       :obj:`ValueError`
           If ``nozzle`` is not ``"rao"`` or ``"conical"``, or if post-processing
           detects non-monotonic ``x`` ordering that cannot be repaired.







   ..
       !! processed by numpydoc !!

.. py:function:: get_theta_e_n(length_fraction, epsilon_value)

   
   Interpolate exit and throat angles from tabulated JSON data.


   :Parameters:

       **length_fraction** : :class:`python:float`
           Normalized nozzle length fraction in ``[0.60, 1.00]``.

       **epsilon_value** : :class:`python:float`
           Nozzle area ratio :math:`\epsilon = A_e / A_t`.



   :Returns:

       :class:`python:tuple`\[:class:`python:float`, :class:`python:float`]
           ``(theta_e, theta_n)`` in **radians**, where ``theta_e`` is the exit
           angle and ``theta_n`` is the throat (diverging-side) angle.








   .. rubric:: Notes

   The function performs a two-stage interpolation:

   1. For each bracketing length fraction dataset, linearly interpolate
      angle vs. area ratio.
   2. Interpolate those two results across length fraction.

   JSON files are read from ``<this module>/data/{theta_e,theta_n}.json``.



   ..
       !! processed by numpydoc !!

.. py:function:: integrate_area(L_c, xs_chamber, ys_chamber)

   
   Integrate the area under ``y(x)`` from ``x=0`` down to ``x=-L_c``.


   :Parameters:

       **L_c** : :class:`python:float`
           Positive cutoff length (m).

       **xs_chamber** : :term:`numpy:array_like`, :obj:`shape` (N,)
           Axial coordinates; must include non-positive values.

       **ys_chamber** : :term:`numpy:array_like`, :obj:`shape` (N,)
           Radii corresponding to ``xs_chamber``.



   :Returns:

       :class:`python:float`
           Area :math:`\int_0^{-L_c} y(x)\,dx` (returned as a positive number).




   :Raises:

       :obj:`ValueError`
           If ``L_c < 0`` or ``-L_c`` lies outside the negative portion of
           ``xs_chamber``.







   ..
       !! processed by numpydoc !!

