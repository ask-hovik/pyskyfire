pyskyfire.skycea.Aerothermodynamics
===================================

.. py:class:: pyskyfire.skycea.Aerothermodynamics(optimum: dict[str, float | str], chemrep_map: Optional[dict[str, str]] = None)

   .. py:method:: _bilinear_map(prop_map: numpy.ndarray, x: float, T: float) -> float | None


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

      
      Equilibrium evaluation at (T, p) using your reactants and MR (no frozen option).
















      ..
          !! processed by numpydoc !!


   .. py:method:: _get_prop(*, x: float, T: float | None, h: float | None, map_attr: str, res_attr: str) -> float


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


   .. py:method:: compute_aerothermodynamics(contour, Nt: int = 64)

      
      Create 2-D property maps on (x, T). Column 0 = equilibrium at that x.
















      ..
          !! processed by numpydoc !!


   .. py:method:: from_F_eps_Lstar(fu, ox, MR, p_c, F, eps, L_star, T_fu_in=298.15, T_ox_in=298.15, p_amb=101300.0, npts=15)
      :classmethod:


      
      Calculate optimal values using thrust, exit pressure and L-star
















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


   .. py:method:: get_mu(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_p(x: float, T: float | None = None, h: float | None = None) -> float


   .. py:method:: get_rho(x: float, T: float | None = None, h: float | None = None) -> float

