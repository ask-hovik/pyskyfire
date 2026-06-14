pyskyfire.pump.impeller_new.Impeller
====================================

.. py:class:: pyskyfire.pump.impeller_new.Impeller(Q: float, H: float, n: float, *, d1: float | None = None, dn: float | None = None, blade_count: int | None = None, beta1_deg: float = 18.0, beta2_deg: float = 25.0, blade_angle_law: BladeAngleLaw = 'cosine', phi1: float = 0.2, fd1: float = 1.05, hub_ratio: float = 0.3, fT: float = 1.1, n_streamlines: int = 5, n_points: int = 96)

   
   First-pass closed radial centrifugal impeller design object.

   The class is deliberately staged internally but convenient externally: on
   initialisation it computes main dimensions and geometry, then stores them in
   explicit dataclasses.  Solvers and visualisation code should consume the
   object through ``dimensions`` and ``geometry`` rather than recomputing.















   ..
       !! processed by numpydoc !!

   .. py:method:: __repr__() -> str


   .. py:method:: __str__() -> str


   .. py:method:: _compute_dimensions() -> ImpellerDimensions


   .. py:method:: _compute_geometry() -> ImpellerGeometry


   .. py:method:: export_json(path: str | pathlib.Path) -> pathlib.Path

      
      Export dimensions and construction curves for CAD/visualisation.
















      ..
          !! processed by numpydoc !!


   .. py:method:: from_pressure_rise(*, mdot: float, rho: float, dp: float, n: float, **kwargs) -> Impeller
      :classmethod:


      
      Construct from mass flow, density and pressure rise.

      This is useful for rocket-engine cycle models where the pump load is
      usually specified as ``dp`` rather than as head.















      ..
          !! processed by numpydoc !!


   .. py:method:: to_dict() -> dict

      
      Return scalar design data plus geometry arrays as JSON-ready lists.
















      ..
          !! processed by numpydoc !!

