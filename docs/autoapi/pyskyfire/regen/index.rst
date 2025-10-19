pyskyfire.regen
===============

.. py:module:: pyskyfire.regen


Submodules
----------

.. toctree::
   :maxdepth: 1

   /autoapi/pyskyfire/regen/channel_height/index
   /autoapi/pyskyfire/regen/contour/index
   /autoapi/pyskyfire/regen/cross_section/index
   /autoapi/pyskyfire/regen/physics/index
   /autoapi/pyskyfire/regen/solver/index
   /autoapi/pyskyfire/regen/thrust_chamber/index


Attributes
----------

.. autoapisummary::

   pyskyfire.regen.ReDh_laminar
   pyskyfire.regen.ReDh_turbulent


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/regen/BoundaryConditions
   /autoapi/pyskyfire/regen/ChannelPlacement
   /autoapi/pyskyfire/regen/ChannelSection
   /autoapi/pyskyfire/regen/Contour
   /autoapi/pyskyfire/regen/Contour
   /autoapi/pyskyfire/regen/Contour
   /autoapi/pyskyfire/regen/ContourToroidalAerospike
   /autoapi/pyskyfire/regen/CoolingCircuit
   /autoapi/pyskyfire/regen/CoolingCircuitGroup
   /autoapi/pyskyfire/regen/CrossSectionRounded
   /autoapi/pyskyfire/regen/CrossSectionSquared
   /autoapi/pyskyfire/regen/HeatExchangerPhysics
   /autoapi/pyskyfire/regen/InternalPlacement
   /autoapi/pyskyfire/regen/SectionProfiles
   /autoapi/pyskyfire/regen/SectionProfiles
   /autoapi/pyskyfire/regen/SurfacePlacement
   /autoapi/pyskyfire/regen/ThrustChamber
   /autoapi/pyskyfire/regen/Wall
   /autoapi/pyskyfire/regen/WallGroup

.. autoapisummary::

   pyskyfire.regen.BoundaryConditions
   pyskyfire.regen.ChannelPlacement
   pyskyfire.regen.ChannelSection
   pyskyfire.regen.Contour
   pyskyfire.regen.Contour
   pyskyfire.regen.Contour
   pyskyfire.regen.ContourToroidalAerospike
   pyskyfire.regen.CoolingCircuit
   pyskyfire.regen.CoolingCircuitGroup
   pyskyfire.regen.CrossSectionRounded
   pyskyfire.regen.CrossSectionSquared
   pyskyfire.regen.HeatExchangerPhysics
   pyskyfire.regen.InternalPlacement
   pyskyfire.regen.SectionProfiles
   pyskyfire.regen.SectionProfiles
   pyskyfire.regen.SurfacePlacement
   pyskyfire.regen.ThrustChamber
   pyskyfire.regen.Wall
   pyskyfire.regen.WallGroup


Functions
---------

.. autoapisummary::

   pyskyfire.regen.T_aw
   pyskyfire.regen.analyse_residuals
   pyskyfire.regen.compute_chamber_volume
   pyskyfire.regen.compute_cutoff_length
   pyskyfire.regen.f_darcy
   pyskyfire.regen.f_darcy_laminar
   pyskyfire.regen.f_darcy_turbulent
   pyskyfire.regen.get_contour
   pyskyfire.regen.get_contour_internal
   pyskyfire.regen.get_theta_e_n
   pyskyfire.regen.h_coolant_colburn
   pyskyfire.regen.h_gas_bartz
   pyskyfire.regen.h_gas_bartz_enthalpy_driven
   pyskyfire.regen.integrate_area
   pyskyfire.regen.interleaved_indices
   pyskyfire.regen.make_channel_height_fn
   pyskyfire.regen.phi_curv
   pyskyfire.regen.radius_of_curvature
   pyskyfire.regen.reynolds
   pyskyfire.regen.sigma
   pyskyfire.regen.solve_heat_exchanger_euler
   pyskyfire.regen.steady_heating_analysis
   pyskyfire.regen.u_coolant


Package Contents
----------------

.. py:function:: T_aw(gamma, M_inf, T_inf, Pr)

   
   Compute the adiabatic-wall (recovery) temperature.


   :Parameters:

       **gamma** : :class:`python:float`
           Ratio of specific heats.

       **M_inf** : :class:`python:float`
           Freestream Mach number.

       **T_inf** : :class:`python:float`
           Freestream static temperature [K].

       **Pr** : :class:`python:float`
           Prandtl number.



   :Returns:

       :class:`python:float`
           Adiabatic-wall temperature [K].








   .. rubric:: Notes

   The recovery factor ``r`` is approximated as ``Pr^(1/3)``,
   valid for turbulent flow. For laminar conditions, a smaller
   recovery factor should be used (≈ Pr⁰˙⁵).

   TODO: Implement a more robust version that could also handle laminar flow, in case this function is suddenly
   used in a function where this makes sense.



   ..
       !! processed by numpydoc !!

.. py:function:: analyse_residuals(residual_log, n_cells, p=2)

   
   Aggregate local solver residuals into global history and final per-cell vector.


   :Parameters:

       **residual_log** : :class:`python:list` or :data:`python:None`
           List of tuples ``(cell, iter, R1, R2)`` recorded during solves.
           If ``None`` or empty, returns ``(None, None)``.

       **n_cells** : :class:`python:int`
           Number of axial cells.

       **p** : :class:`python:int` or :class:`python:float`, :obj:`optional`
           Norm order for global residual history: ``2`` for RMS,
           ``np.inf`` for :math:`L_\infty`, etc. Default is 2.



   :Returns:

       **history** : :obj:`ndarray <numpy.ndarray>` or :data:`python:None`
           Global residual norm for iterations ``0..k_max``, or ``None``.

       **final_per_cell** : :obj:`ndarray <numpy.ndarray>` or :data:`python:None`
           Final residual magnitude per cell at its last local iteration, or ``None``.











   ..
       !! processed by numpydoc !!

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

.. py:function:: f_darcy(ReDh, Dh, x, roughness)

   
   Composite laminar–turbulent Darcy friction factor with smooth transition.


   :Parameters:

       **ReDh** : :class:`python:float`
           Hydraulic Reynolds number.

       **Dh** : :class:`python:float`
           Hydraulic diameter [m].

       **x** : :class:`python:float`
           Axial coordinate.

       **roughness** : :func:`python:callable` or :data:`python:None`
           Function returning surface roughness at ``x`` [m], or ``None`` for smooth wall.



   :Returns:

       :class:`python:float`
           Darcy friction factor (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: f_darcy_laminar(ReDh, Dh, x)

   
   Laminar Darcy friction factor.


   :Parameters:

       **ReDh** : :class:`python:float`
           Hydraulic Reynolds number.

       **Dh** : :class:`python:float`
           Hydraulic diameter [m].

       **x** : :class:`python:float`
           Axial coordinate (unused, kept for API consistency).



   :Returns:

       :class:`python:float`
           Darcy friction factor (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: f_darcy_turbulent(ReDh, Dh, x, roughness)

   
   Turbulent Darcy friction factor using Putukhov or Colebrook–White.


   :Parameters:

       **ReDh** : :class:`python:float`
           Hydraulic Reynolds number.

       **Dh** : :class:`python:float`
           Hydraulic diameter [m].

       **x** : :class:`python:float`
           Axial coordinate.

       **roughness** : :func:`python:callable` or :data:`python:None`
           Function returning surface roughness at ``x`` [m], or ``None`` for smooth wall.



   :Returns:

       :class:`python:float`
           Darcy friction factor (dimensionless).











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

.. py:function:: h_coolant_colburn(k_cf, D_c, Cp_cr, mu_cf, mdot_c, A_c, phi_curv=1)

   
   Colburn correlation for turbulent coolant-side heat transfer.


   :Parameters:

       **k_cf** : :class:`python:float`
           Thermal conductivity of coolant film [W m⁻¹ K⁻¹].

       **D_c** : :class:`python:float`
           Coolant hydraulic diameter [m].

       **Cp_cr** : :class:`python:float`
           Mean specific heat (reference enthalpy) [J kg⁻¹ K⁻¹].

       **mu_cf** : :class:`python:float`
           Dynamic viscosity of coolant film [Pa s].

       **mdot_c** : :class:`python:float`
           Coolant mass-flow rate [kg s⁻¹].

       **A_c** : :class:`python:float`
           Coolant-flow area [m²].

       **phi_curv** : :class:`python:float`, :obj:`optional`
           Curvature-correction factor (dimensionless). Default is 1.



   :Returns:

       :class:`python:float`
           Coolant-side heat-transfer coefficient [W m⁻² K⁻¹].








   .. rubric:: Notes

   Subscripts: 
       Cf denotes film coolant film condition
       c denotes bulk coolant conditions
       r denotes reference enthalpy condition which are averaged between the free stream and wall metal conditions



   ..
       !! processed by numpydoc !!

.. py:function:: h_gas_bartz(D_t, mu_g, cp_g, Pr_g, p_c, c_star, A_t, A_x, sigma)

   
   Compute the Bartz hot-gas correlation without curvature correction.


   :Parameters:

       **D_t** : :class:`python:float`
           Throat diameter [m].

       **mu_g** : :class:`python:float`
           Gas viscosity [Pa s].

       **cp_g** : :class:`python:float`
           Gas specific heat capacity [J kg⁻¹ K⁻¹].

       **Pr_g** : :class:`python:float`
           Gas Prandtl number.

       **p_c** : :class:`python:float`
           Chamber pressure [Pa].

       **c_star** : :class:`python:float`
           Characteristic velocity [m s⁻¹].

       **A_t** : :class:`python:float`
           Throat area [m²].

       **A_x** : :class:`python:float`
           Local flow area [m²].

       **sigma** : :class:`python:float`
           Boundary-layer property-variation correction (dimensionless).



   :Returns:

       :class:`python:float`
           Heat-transfer coefficient [W m⁻² K⁻¹].











   ..
       !! processed by numpydoc !!

.. py:function:: h_gas_bartz_enthalpy_driven(k_gr, D_hyd, Cp_gr, mu_gr, mdot_g, A_chmb, T_g, T_gr)

   
   Compute the hot-gas-side heat-transfer coefficient (Bartz correlation, enthalpy-based).


   :Parameters:

       **k_gr** : :class:`python:float`
           Thermal conductivity at the reference enthalpy condition [W m⁻¹ K⁻¹].

       **D_hyd** : :class:`python:float`
           Hydraulic diameter of the gas passage [m].

       **Cp_gr** : :class:`python:float`
           Specific heat capacity at the reference enthalpy condition [J kg⁻¹ K⁻¹].

       **mu_gr** : :class:`python:float`
           Dynamic viscosity at the reference condition [Pa s].

       **mdot_g** : :class:`python:float`
           Total hot-gas mass flow rate [kg s⁻¹].

       **A_chmb** : :class:`python:float`
           Flow cross-sectional area in the chamber [m²].

       **T_g** : :class:`python:float`
           Free-stream gas temperature [K].

       **T_gr** : :class:`python:float`
           Reference-enthalpy temperature [K].



   :Returns:

       :class:`python:float`
           Hot-gas-side heat-transfer coefficient ``h_g`` [W m⁻² K⁻¹].








   .. rubric:: Notes

   This form of the Bartz correlation uses the ratio ``T_g/T_gr`` to
   account for variable-property effects between free-stream and reference conditions.

   subscript g denotes free stream gas properties
   subscript r denotes reference enthalpy conditions which are averaged between the free stream and wall metal conditions



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

.. py:function:: interleaved_indices(circuit_counts)

   
   Given a list of circuit_counts = [n0, n1, ..., nK], produce an array
   'owners' of length sum(circuit_counts), where each index i is assigned
   to exactly one circuit in an interleaved ratio of n0 : n1 : ... : nK.

   Example: circuit_counts = [30, 60].
   Then we have total=90, ratio=1:2.  The owners array might look like
     [0,1,1, 0,1,1, 0,1,1, ...]
   So circuit #0 gets 30 slots, circuit #1 gets 60 slots, interleaved 1:2.















   ..
       !! processed by numpydoc !!

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

.. py:function:: phi_curv(Re_c, D_c, R_curv)

   
   Curvature correction factor for the coolant-side heat-transfer coefficient.


   :Parameters:

       **Re_c** : :class:`python:float`
           Coolant Reynolds number.

       **D_c** : :class:`python:float`
           Hydraulic diameter [m].

       **R_curv** : :class:`python:float`
           Radius of curvature of the coolant passage [m].
           Use ``np.inf`` for straight sections.



   :Returns:

       :class:`python:float`
           Curvature factor ``φ`` (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: radius_of_curvature(points: numpy.ndarray, axis: str = 'x', eps: float = 1e-12) -> numpy.ndarray

   
   Signed radius of curvature for a curve expressed in cylindrical coordinates
   [x, r, θ].

   • Positive  → curve bends *away* from the symmetry axis  
   • Negative  → curve bends *toward* the symmetry axis  
   • np.inf    → locally straight (|κ| below `eps`)

   :Parameters:

       **points** : (:obj:`N`, 3) :obj:`ndarray <numpy.ndarray>`
           [[x, r, theta], …] ordered along the curve.

       **axis** : {'x', 'y', 'z'}, :obj:`optional`
           Which coordinate is the symmetry axis.  Default 'x'.

       **eps** : :class:`python:float`, :obj:`optional`
           Curvature values with |κ| < eps are treated as zero (straight).



   :Returns:

       **R** : (N,) :obj:`ndarray <numpy.ndarray>`
           Signed radius of curvature at each sample.











   ..
       !! processed by numpydoc !!

.. py:function:: reynolds(rho, u, L, mu)

   
   Compute the Reynolds number.


   :Parameters:

       **rho** : :class:`python:float`
           Fluid density [kg m⁻³].

       **u** : :class:`python:float`
           Mean velocity [m s⁻¹].

       **L** : :class:`python:float`
           Characteristic length or hydraulic diameter [m].

       **mu** : :class:`python:float`
           Dynamic viscosity [Pa s].



   :Returns:

       :class:`python:float`
           Reynolds number (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: sigma(T_hw, T_c, gamma_g, M_g, omega)

   
   Dimensionless property-variation factor used in Bartz correlation.


   :Parameters:

       **T_hw** : :class:`python:float`
           Wall temperature [K].

       **T_c** : :class:`python:float`
           Core-flow (free-stream) temperature [K].

       **gamma_g** : :class:`python:float`
           Ratio of specific heats.

       **M_g** : :class:`python:float`
           Local Mach number.

       **omega** : :class:`python:float`
           Empirical property exponent (≈ 0.68 for diatomic gases).



   :Returns:

       :class:`python:float`
           ``σ`` correction factor (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: solve_heat_exchanger_euler(thrust_chamber, boundary_conditions, n_nodes, circuit_index, output, log_residuals=True)

   
   Solve 1-D steady heating with a marching Euler scheme.


   :Parameters:

       **thrust_chamber** : :obj:`Any`
           Chamber model exposing geometry and property methods.

       **boundary_conditions** : :obj:`BoundaryConditions`
           Inlet temperature/pressure/mass-flow conditions.

       **n_nodes** : :class:`python:int`
           Number of axial nodes.

       **circuit_index** : :class:`python:int`
           Which cooling circuit to simulate.

       **output** : :ref:`bool <python:bltin-boolean-values>`
           If True, print progress to stdout.

       **log_residuals** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           If True, record local residuals each Newton iteration per cell.



   :Returns:

       :class:`python:dict`
           Results with keys:
           
           - ``x`` : axial coordinates [m]
           - ``T`` : temperatures array, shape ``(n_nodes, 1 + n_walls + 1)``
             (coolant + reversed wall interfaces from cold→hot)
           - ``T_static`` : coolant static temperature [K]
           - ``T_stagnation`` : coolant stagnation temperature [K]
           - ``p_static`` : coolant static pressure [Pa]
           - ``p_stagnation`` : coolant stagnation pressure [Pa]
           - ``dQ_dA`` : local heat flux [W m⁻²]
           - ``velocity`` : coolant velocity [m s⁻¹]
           - ``residuals`` : tuple of (global history, final-per-cell) or (None, None)











   ..
       !! processed by numpydoc !!

.. py:function:: steady_heating_analysis(thrust_chamber, boundary_conditions, n_nodes=100, circuit_index=0, solver='newton', output=True)

   
   Run the steady heating analysis.


   :Parameters:

       **thrust_chamber** : :obj:`Any`
           Chamber model exposing the required geometry & property APIs.

       **boundary_conditions** : :obj:`BoundaryConditions`
           Coolant inlet boundary conditions.

       **n_nodes** : :class:`python:int`, :obj:`optional`
           Number of axial nodes. Default is 100.

       **circuit_index** : :class:`python:int`, :obj:`optional`
           Cooling-circuit index. Default is 0.

       **solver** : {'newton'}, :obj:`optional`
           Solver selector. Currently only ``'newton'`` is implemented.

       **output** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           If True, print progress. Default is True.



   :Returns:

       :class:`python:dict`
           See :func:`solve_heat_exchanger_euler` for keys.











   ..
       !! processed by numpydoc !!

.. py:function:: u_coolant(rho, mdot_c_single_channel, A_cool)

   
   Compute coolant velocity.


   :Parameters:

       **rho** : :class:`python:float`
           Coolant density [kg m⁻³].

       **mdot_c_single_channel** : :class:`python:float`
           Mass-flow rate through one channel [kg s⁻¹].

       **A_cool** : :class:`python:float`
           Coolant cross-sectional area [m²].



   :Returns:

       :class:`python:float`
           Coolant velocity [m s⁻¹].











   ..
       !! processed by numpydoc !!

.. py:data:: ReDh_laminar
   :value: 2300


.. py:data:: ReDh_turbulent
   :value: 3500


