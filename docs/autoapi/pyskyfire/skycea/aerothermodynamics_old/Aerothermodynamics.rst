pyskyfire.skycea.aerothermodynamics_old.Aerothermodynamics
==========================================================

.. py:class:: pyskyfire.skycea.aerothermodynamics_old.Aerothermodynamics(optimum: dict[str, float | str], chemrep_map: Optional[dict[str, str]] = None)

   
   Aerothermodynamic property precomputation and lookup along an engine contour.

   This class computes 2-D thermo/fluid property maps along a prescribed
   nozzle/combustor contour using **CEA_Wrap** and exposes fast interpolating
   getters. At each axial position ``x`` the equilibrium state (column 0) is
   evaluated at the local area ratio; additional columns tabulate thermodynamic
   properties over a temperature grid at **fixed local static pressure**. Lookups
   support:
   - equilibrium at a given ``x`` (no ``T``/``h`` provided),
   - TP evaluation at specified ``T`` and local pressure ``p(x)``,
   - HP evaluation at specified mixture enthalpy ``h`` and local pressure ``p(x)``
   (using an adjustable-enthalpy reactant representation).

   :Parameters:

       **optimum** : :class:`python:dict`\[:class:`python:str`, :class:`python:float` | :class:`python:str`]
           Calibrated/derived inputs (e.g., ``fu``, ``ox``, ``MR``, ``p_c``, ``F``,
           ``eps``, ``L_star``, ``c_star``, inlet temperatures, etc.). Typically
           produced by :meth:`from_F_eps_Lstar`.

       **chemrep_map** : :class:`python:dict`\[:class:`python:str`, :class:`python:str`], :obj:`optional`
           Mapping from reactant *names* to ``chemical_representation`` strings used by
           CEA for HP evaluations with enthalpy offsets (exploded formula like
           ``"C 2 H 6 O 1"``). Required only when calling HP-based getters (``h=...``).

   :Attributes:

       **# Inputs / design-point scalars**
           ..

       **fu** : :obj:`Any`
           User’s fuel specification (mixture container understood by CEA_Wrap).

       **ox** : :obj:`Any`
           User’s oxidizer specification (mixture container understood by CEA_Wrap).

       **MR** : :class:`python:float`
           Mixture ratio ``m_ox / m_fu`` [-].

       **p_c** : :class:`python:float`
           Chamber pressure [Pa].

       **F** : :class:`python:float`
           Target thrust [N].

       **eps** : :class:`python:float`
           Exit-to-throat area ratio ``A_e/A_t`` [-].

       **L_star** : :class:`python:float`
           Characteristic chamber length [m].

       **T_fu_in, T_ox_in** : :class:`python:float`
           Inlet temperatures for fuel/oxidizer [K].

       **p_amb** : :class:`python:float`
           Ambient/static back pressure used for CF/Isp_amb calculations [Pa].

       **npts** : :class:`python:int`
           Number of axial nodes used for precomputation along the contour [-].

       **# Derived performance at design point**
           ..

       **c_star** : :class:`python:float`
           Characteristic velocity [m/s].

       **Isp_vac, Isp_amb, Isp_SL, Isp_ideal_amb** : :class:`python:float`
           Vacuum / ambient / sea-level / perfectly expanded specific impulses [s].

       **CF_vac, CF_amb, CF_SL** : :class:`python:float`
           Thrust coefficients [-].

       **mdot, mdot_fu, mdot_ox** : :class:`python:float`
           Total/fuel/oxidizer mass flows [kg/s].

       **A_t, A_e** : :class:`python:float`
           Throat and exit areas [m²].

       **r_t, r_e** : :class:`python:float`
           Throat and exit radii [m].

       **V_c** : :class:`python:float`
           Chamber volume estimated from ``L_star`` [m³].

       **t_stay** : :class:`python:float`
           Mean residence time in chamber [s].

       **# Precomputed grids (after :meth:`compute_aerothermodynamics`)**
           ..

       **x_nodes** : :obj:`ndarray <numpy.ndarray>`
           Monotone axial grid spanning ``contour.xs[0]`` → ``contour.xs[-1]`` [m].

       **Nt** : :class:`python:int`
           Number of temperature samples per axial station [-].

       **T_grid** : :obj:`ndarray <numpy.ndarray>`, :obj:`shape` (:obj:`Nx`, :obj:`Nt`)
           Temperature grid per axial row [K] (col 0 equals local equilibrium T).

       **M_map, T_map, p_map, rho_map, cp_map, gamma_map, h_map, a_map, mu_map, k_map, Pr_map** : :obj:`ndarray <numpy.ndarray>`
           Property maps; column 0 is the local equilibrium; remaining columns are TP
           samples at fixed local pressure. Units: M [-], T [K], p [bar], rho [kg/m³],
           cp [kJ/kg-K], gamma [-], h [kJ/kg], a [m/s], mu [Pa·s], k [W/m-K], Pr [-].

       **X_map** : :class:`python:list`\[:class:`python:dict`\[:class:`python:str`, :class:`python:float`]]
           Equilibrium product mole fractions per axial node (unnormalized dicts).









   .. seealso::

       
       :obj:`CEA_Wrap.RocketProblem`
           Equilibrium nozzle/chamber solves used for the equilibrium column.
       :obj:`CEA_Wrap.TPProblem`
           Temperature-pressure state solves used for the TP columns and fallbacks.
       :obj:`CEA_Wrap.HPProblem`
           Enthalpy-pressure equilibrium solves used when ``h`` is specified.
       
       
   .. rubric:: Notes

   - **Units:** Inputs use SI (Pa, K, N). CEA_Wrap expects/returns some fields in
   *psi* / *bar*; all conversions are handled internally. Stored ``p_map`` is in
   **bar** (matching CEA); :meth:`get_p` returns **Pa**.
   - **Equilibrium column:** column 0 in each map corresponds to the equilibrium
   state at the local area ratio (subsonic for ``x < 0``; supersonic for
   ``x >= 0`` with your sign convention).
   - **Interpolation:** Lookups first try a bilinear map interpolation on
   ``(x, T)``. If the requested ``T`` falls outside the tabulated range, a live
   TP solve at ``p(x)`` is performed. For ``h`` queries an HP solve at ``p(x)``
   is performed and requires valid ``chemrep_map`` entries for all reactants.



   ..
       !! processed by numpydoc !!

   .. py:method:: _bilinear_map(prop_map: numpy.ndarray, x: float, T: float) -> float | None

      
      Interpolate a TP property field in x and T.

      Temperatures below T_low clamp to the value at T_low.
      Temperatures above T_hi return None so a live TP calculation can be used.















      ..
          !! processed by numpydoc !!


   .. py:method:: _evaluate_hp(*, h_target_Jkg: float, p: float)

      
      Evaluate an equilibrium HP state at local pressure p [Pa] and mixture specific enthalpy h [J/kg].
      Uses the 'chemical_representation + hf' override path.
      Strategy: assign zero enthalpy to all reactants EXCEPT one 'adjuster' fuel species,
      whose molar enthalpy is chosen so that the mixture enthalpy equals h_target_Jkg.

      Requirements:
        - self.chemrep_map must provide exploded formulas for any species we override.
        - We'll pick the FIRST fuel in self.fu as the 'adjuster'.















      ..
          !! processed by numpydoc !!


   .. py:method:: _evaluate_tp(*, T: float, p: float)

      
      Run a live TP calculation at temperature T [K] and pressure p [Pa].
















      ..
          !! processed by numpydoc !!


   .. py:method:: _get_prop(*, x: float, T: float | None, h: float | None, eq_values: numpy.ndarray, prop_map: numpy.ndarray, res_attr: str, live_scale: float = 1.0) -> float


   .. py:method:: _interp_X_dict(x: float) -> dict[str, float]

      
      Interpolate mole-fraction dict at position x from self.X_map (list[dict]).
      Uses the union of species at the two bracketing nodes, fills missing with 0,
      clips negatives, and renormalizes to sum = 1.
















      ..
          !! processed by numpydoc !!


   .. py:method:: _interp_eq_column(map2d: numpy.ndarray, x: float) -> float

      
      Interpolate along x using the equilibrium column (col 0).
















      ..
          !! processed by numpydoc !!


   .. py:method:: _interp_scalar(x: float, xs: numpy.ndarray, ys: numpy.ndarray) -> float

      
      Linear interpolation with endpoint clamping.
















      ..
          !! processed by numpydoc !!


   .. py:method:: _make_cea_reactants()

      
      Create fresh CEA_Wrap reactant objects from the stored mixtures.
















      ..
          !! processed by numpydoc !!


   .. py:method:: compute_aerothermodynamics(contour, Nt: int = 64, *, T_low: float = 500.0, T_hi: float | None = None)

      
      Precompute equilibrium states and common-grid TP property fields.

      The RocketProblem solution defines the equilibrium nozzle trajectory.
      The TP fields are evaluated at the local equilibrium pressure using one
      shared temperature axis from T_low to T_hi.












      .. rubric:: Notes

      - ``*_eq`` arrays hold RocketProblem equilibrium values along x.
      - ``*_map`` arrays hold TPProblem values on a common x-T grid.
      - Temperatures below T_low are clamped during property lookup.



      ..
          !! processed by numpydoc !!


   .. py:method:: from_F_eps_Lstar(fu, ox, MR, p_c, F, eps, L_star, T_fu_in=298.15, T_ox_in=298.15, p_amb=101300.0, npts=15)
      :classmethod:


      
      Construct from thrust, area ratio, and L* at a given chamber pressure.

      This helper solves a CEA rocket problem at the specified design point and
      assembles the ``optimum`` dict used to initialize :class:`Aerothermodynamics`.

      :Parameters:

          **fu** : :obj:`Any`
              Fuel mixture descriptor (must expose ``propellants`` and ``fractions`` and
              be consumable by ``CEA_Wrap.Fuel``).

          **ox** : :obj:`Any`
              Oxidizer mixture descriptor (same structural expectations as ``fu``).

          **MR** : :class:`python:float`
              Mixture ratio ``m_ox / m_fu`` [-].

          **p_c** : :class:`python:float`
              Chamber pressure [Pa].

          **F** : :class:`python:float`
              Target thrust [N].

          **eps** : :class:`python:float`
              Nozzle area ratio ``A_e / A_t`` [-].

          **L_star** : :class:`python:float`
              Characteristic chamber length [m].

          **T_fu_in** : :class:`python:float`, :obj:`default` 298.15
              Fuel inlet temperature [K].

          **T_ox_in** : :class:`python:float`, :obj:`default` 298.15
              Oxidizer inlet temperature [K].

          **p_amb** : :class:`python:float`, :obj:`default` 1.013e5
              Ambient/static pressure for CF/Isp_amb calculations [Pa].

          **npts** : :class:`python:int`, :obj:`default` 15
              Number of axial stations to precompute along the contour.



      :Returns:

          :obj:`Aerothermodynamics`
              Initialized instance with populated design-point performance and all fields
              necessary for :meth:`compute_aerothermodynamics`.








      .. rubric:: Notes

      - The CEA call assumes perfect expansion for ``Isp_vac`` and uses standard
      thrust-coefficient relations to compute ``Isp_amb`` and ``Isp_SL``.
      - Areas and chamber volume are derived from ``c_star``, mass flow, ``L_star``,
      and mixture density at chamber conditions.



      ..
          !! processed by numpydoc !!


   .. py:method:: from_F_pe_Lstar(fu, ox, MR, p_c, F, p_e, L_star, T_fu_in=298.15, T_ox_in=298.15, p_amb=101300.0, npts=15)
      :classmethod:


      
      Construct from thrust, area ratio, and L* at a given chamber pressure.

      This helper solves a CEA rocket problem at the specified design point and
      assembles the ``optimum`` dict used to initialize :class:`Aerothermodynamics`.

      :Parameters:

          **fu** : :obj:`Any`
              Fuel mixture descriptor (must expose ``propellants`` and ``fractions`` and
              be consumable by ``CEA_Wrap.Fuel``).

          **ox** : :obj:`Any`
              Oxidizer mixture descriptor (same structural expectations as ``fu``).

          **MR** : :class:`python:float`
              Mixture ratio ``m_ox / m_fu`` [-].

          **p_c** : :class:`python:float`
              Chamber pressure [Pa].

          **F** : :class:`python:float`
              Target thrust [N].

          **eps** : :class:`python:float`
              Nozzle area ratio ``A_e / A_t`` [-].

          **L_star** : :class:`python:float`
              Characteristic chamber length [m].

          **T_fu_in** : :class:`python:float`, :obj:`default` 298.15
              Fuel inlet temperature [K].

          **T_ox_in** : :class:`python:float`, :obj:`default` 298.15
              Oxidizer inlet temperature [K].

          **p_amb** : :class:`python:float`, :obj:`default` 1.013e5
              Ambient/static pressure for CF/Isp_amb calculations [Pa].

          **npts** : :class:`python:int`, :obj:`default` 15
              Number of axial stations to precompute along the contour.



      :Returns:

          :obj:`Aerothermodynamics`
              Initialized instance with populated design-point performance and all fields
              necessary for :meth:`compute_aerothermodynamics`.








      .. rubric:: Notes

      - The CEA call assumes perfect expansion for ``Isp_vac`` and uses standard
      thrust-coefficient relations to compute ``Isp_amb`` and ``Isp_SL``.
      - Areas and chamber volume are derived from ``c_star``, mass flow, ``L_star``,
      and mixture density at chamber conditions.



      ..
          !! processed by numpydoc !!


   .. py:method:: get_H(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_M(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_Pr(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_T(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_X(x: float) -> dict[str, float]

      
      Interpolated mole-fraction dict at position x.
















      ..
          !! processed by numpydoc !!


   .. py:method:: get_a(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_cp(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_gamma(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_h(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_k(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_molecular_weight(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_mu(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_p(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_rho(x: float, T: float | None = None, h: float | None = None) -> float

