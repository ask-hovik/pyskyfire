pyskyfire.pump.impeller_new
===========================

.. py:module:: pyskyfire.pump.impeller_new

.. autoapi-nested-parse::

   Preliminary centrifugal impeller geometry generation for pyskyfire.

   This module is intentionally a geometry/design-preparation layer, not a pump
   meanline/performance solver.  It sizes a first-pass closed radial centrifugal
   impeller from a design point, generates meridional streamlines, integrates a
   blade camber surface from a blade-angle law, and exposes arrays suitable for
   visualisation/meshing/CAD export.

   Coordinate convention
   ---------------------
   The internal streamline convention follows the rest of pyskyfire's geometry
   style: cylindrical coordinates are stored as ``(x, r, theta)`` where ``x`` is the
   axial coordinate [m], ``r`` is radius [m], and ``theta`` is angle [rad].

   The implemented correlations are preliminary engineering correlations based on
   Gülich-style pump sizing.  They are appropriate for smooth radial, closed,
   single-stage pump geometry generation at BEP-like design points; they are not a
   validated efficiency/cavitation/stress model.

   ..
       !! processed by numpydoc !!


Attributes
----------

.. autoapisummary::

   pyskyfire.pump.impeller_new.BladeAngleLaw
   pyskyfire.pump.impeller_new.G0
   pyskyfire.pump.impeller_new.TWO_PI
   pyskyfire.pump.impeller_new._RA_STAR
   pyskyfire.pump.impeller_new._RI_STAR
   pyskyfire.pump.impeller_new._ZA_STAR
   pyskyfire.pump.impeller_new._ZI_STAR


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/pump/impeller_new/Impeller
   /autoapi/pyskyfire/pump/impeller_new/ImpellerDimensions
   /autoapi/pyskyfire/pump/impeller_new/ImpellerGeometry
   /autoapi/pyskyfire/pump/impeller_new/ImpellerInputs

.. autoapisummary::

   pyskyfire.pump.impeller_new.Impeller
   pyskyfire.pump.impeller_new.ImpellerDimensions
   pyskyfire.pump.impeller_new.ImpellerGeometry
   pyskyfire.pump.impeller_new.ImpellerInputs


Functions
---------

.. autoapisummary::

   pyskyfire.pump.impeller_new.add_blade_theta
   pyskyfire.pump.impeller_new.blade_angle_distribution
   pyskyfire.pump.impeller_new.cylindrical_to_cartesian
   pyskyfire.pump.impeller_new.estimate_blade_count
   pyskyfire.pump.impeller_new.inlet_diameter_from_phi
   pyskyfire.pump.impeller_new.interpolate_curve_by_arclength
   pyskyfire.pump.impeller_new.meridional_streamlines
   pyskyfire.pump.impeller_new.outlet_diameter
   pyskyfire.pump.impeller_new.outlet_width_ratio
   pyskyfire.pump.impeller_new.pressure_coefficient
   pyskyfire.pump.impeller_new.replicate_blades
   pyskyfire.pump.impeller_new.specific_speed
   pyskyfire.pump.impeller_new.tangential_velocity


Module Contents
---------------

.. py:function:: add_blade_theta(meridional: numpy.ndarray, *, beta1_deg: float, beta2_deg: float, law: BladeAngleLaw = 'cosine') -> tuple[numpy.ndarray, numpy.ndarray]

   
   Integrate blade wrap angle on each meridional streamline.

   Uses ``dtheta = dm / (r tan(beta))``.  This is a geometric camber-surface
   construction law, not a solved velocity-triangle model.















   ..
       !! processed by numpydoc !!

.. py:function:: blade_angle_distribution(n_points: int, *, beta1_deg: float, beta2_deg: float, law: BladeAngleLaw = 'cosine') -> numpy.ndarray

   
   Blade metal angle distribution from inlet to outlet [deg].
















   ..
       !! processed by numpydoc !!

.. py:function:: cylindrical_to_cartesian(cyl: numpy.ndarray) -> numpy.ndarray

   
   Convert ``(..., 3)`` array from ``(x, r, theta)`` to ``(x, y, z)``.
















   ..
       !! processed by numpydoc !!

.. py:function:: estimate_blade_count(nq: float) -> tuple[int, str]

   
   Return a conservative preliminary main-blade count.

   Blade count is usually a design decision, not a reliable scalar correlation.
   The heuristic below exists to keep geometry generation automatic; serious
   designs should pass ``blade_count`` explicitly.















   ..
       !! processed by numpydoc !!

.. py:function:: inlet_diameter_from_phi(*, Q: float, n: float, phi1: float, hub_ratio: float, fd1: float = 1.05) -> float

   
   Estimate inlet diameter from continuity and an inlet flow coefficient.

   ``Q = c_m1 A1`` with ``c_m1 = phi1 u1`` and
   ``A1 = pi/4 * d1^2 * (1 - hub_ratio^2)``.















   ..
       !! processed by numpydoc !!

.. py:function:: interpolate_curve_by_arclength(points: numpy.ndarray, n_points: int) -> numpy.ndarray

   
   Resample a 2D/3D polyline to approximately uniform arc length.
















   ..
       !! processed by numpydoc !!

.. py:function:: meridional_streamlines(*, nq: float, b2: float, d1: float, dn: float, d2: float, n_streamlines: int = 5, n_points: int = 96) -> numpy.ndarray

   
   Generate hub-to-shroud meridional streamlines ``(x, r)``.





   :Returns:

       :obj:`np.ndarray <numpy.ndarray>`
           Shape ``(n_streamlines, n_points, 2)``.  The first streamline is the hub
           side, the last streamline is the shroud side.  Points run from inlet to
           outlet.











   ..
       !! processed by numpydoc !!

.. py:function:: outlet_diameter(psi: float, n: float, H: float) -> float

   
   Outlet diameter from ``psi = 2 g H / u2^2`` and ``u2 = pi d2 n/60``.
















   ..
       !! processed by numpydoc !!

.. py:function:: outlet_width_ratio(nq: float) -> float

   
   Approximate impeller outlet width ratio ``b2/d2``.

   This preserves the polynomial already used in the prototype code, with the
   reference value made explicit.















   ..
       !! processed by numpydoc !!

.. py:function:: pressure_coefficient(nq: float, *, fT: float = 1.1) -> float

   
   Approximate optimum pressure/head coefficient.

   The coefficient convention used here is ``psi = 2 g H / u2^2``, matching
   the outlet-diameter expression below.  ``fT`` is retained as an explicit
   tuning factor rather than buried as a magic number.















   ..
       !! processed by numpydoc !!

.. py:function:: replicate_blades(camber_surface: numpy.ndarray, blade_count: int) -> list[numpy.ndarray]

   
   Replicate one blade camber surface around the pump axis.
















   ..
       !! processed by numpydoc !!

.. py:function:: specific_speed(n: float, Q: float, H: float) -> float

   
   Gülich-style pump specific speed ``n_q`` [rpm, m^3/s, m].
















   ..
       !! processed by numpydoc !!

.. py:function:: tangential_velocity(n: float, d: float) -> float

   
   Circumferential speed [m/s] for rpm and diameter [m].
















   ..
       !! processed by numpydoc !!

.. py:data:: BladeAngleLaw

.. py:data:: G0
   :value: 9.80665


.. py:data:: TWO_PI
   :value: 6.283185307179586


.. py:data:: _RA_STAR

.. py:data:: _RI_STAR

.. py:data:: _ZA_STAR

.. py:data:: _ZI_STAR

