pyskyfire.skycea
================

.. py:module:: pyskyfire.skycea


Submodules
----------

.. toctree::
   :maxdepth: 1

   /api/pyskyfire/skycea/aerothermodynamics/index
   /api/pyskyfire/skycea/coolant_transport/index
   /api/pyskyfire/skycea/nozzle_solver/index






Package Contents
----------------

.. py:class:: Fluid(type, propellants, fractions, basis='mass', precision=3)

   
   :param type: "fuel" | "oxidizer" | "coolant" etc.
   :type type: str
   :param propellants: Component names (CoolProp canonical names)
   :type propellants: list[str]
   :param fractions: Fractions (mass or mole, depending on basis)
   :type fractions: list[float]
   :param basis: "mass" or "mole"
   :type basis: str
   :param precision: number of decimal digits when exporting mole fractions
   :type precision: int


   .. py:attribute:: type


   .. py:attribute:: propellants


   .. py:attribute:: fractions


   .. py:attribute:: basis
      :value: 'mass'



   .. py:attribute:: precision
      :value: 3



   .. py:method:: molar_masses()

      Return molar masses [kg/mol] for each propellant.



   .. py:method:: as_mole_fractions()

      Convert stored fractions to mole fractions.



   .. py:method:: coolprop_string()

      Return a CoolProp HEOS string like 'HEOS::Ethanol[0.800]&Water[0.200]'.



.. py:class:: TransportProperties(Pr, mu, k, cp=None, rho=None, gamma_coolant=None)

   
   Container for specifying your transport properties. Each input can either be a function of temperature (K) and pressure (Pa) in that order, e.g. mu(T, p). Otherwise they can be constant floats.

   :param Pr: Prandtl number.
   :type Pr: float or callable
   :param mu: Absolute viscosity (Pa s).
   :type mu: float or callable
   :param k: Thermal conductivity (W/m/K).
   :type k: float or callable
   :param cp: Isobaric specific heat capacity (J/kg/K) - only required for coolants.
   :type cp: float or callable, optional
   :param rho: Density (kg/m^3) - only required for coolants.
   :type rho: float or callable, optional
   :param gamma_coolant: Ratio of specific heats (cp/cv) for a compressible coolant. If this is submitted, it is assumed that this object represents a compressible coolant.
   :type gamma_coolant: float or callable, optional

   .. attribute:: compressible_coolant

      Whether or not this TransportProperties object represents a compressible coolant.

      :type: bool


   .. py:attribute:: type


   .. py:method:: Pr(T, p)

      Prandtl number.

      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Prandtl number
      :rtype: float



   .. py:method:: mu(T, p)

      Absolute viscosity (Pa s)

      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Absolute viscosity (Pa s)
      :rtype: float



   .. py:method:: k(T, p)

      Thermal conductivity (W/m/K)

      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Thermal conductivity (W/m/K)
      :rtype: float



   .. py:method:: rho(T, p)

      Density (kg/m^3)
      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Density (kg/m^3)
      :rtype: float



   .. py:method:: cp(T, p)

      Isobaric specific heat capacity (J/kg/K)

      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Isobaric specific heat capacity (J/kg/K)
      :rtype: float



   .. py:method:: gamma_coolant(T, p)

      Ratio of specific heat capacities for a compressible coolant.

      :param T: Temperature (K)
      :type T: float
      :param p: Pressure (Pa)
      :type p: float

      :returns: Ratio of specific heat capacities (cp/cv).
      :rtype: float



.. py:class:: CoolantTransport(fluid)

   .. py:attribute:: fluid


   .. py:method:: get_Pr(T, p)


   .. py:method:: get_mu(T, p)


   .. py:method:: get_k(T, p)


   .. py:method:: get_cp(T, p)


   .. py:method:: get_rho(T, p)


   .. py:method:: get_cv(T, p)


   .. py:method:: get_gamma(T, p)


.. py:data:: script_dir

.. py:data:: trans_path

.. py:class:: Aerothermodynamics(optimum, chemrep_map = None)

   
   Add all values in the optimum dict to self


   .. py:attribute:: optimum


   .. py:attribute:: chemrep_map


   .. py:method:: from_F_eps_Lstar(fu, ox, MR, p_c, F, eps, L_star, T_fu_in=298.15, T_ox_in=298.15, p_amb=101300.0, npts=15)
      :classmethod:


      Calculate optimal values using thrust, exit pressure and L-star



   .. py:method:: compute_aerothermodynamics(contour, Nt = 64)

      Create 2-D property maps on (x, T). Column 0 = equilibrium at that x.



   .. py:method:: get_T(x, T = None, h = None)


   .. py:method:: get_p(x, T = None, h = None)


   .. py:method:: get_M(x, T = None, h = None)


   .. py:method:: get_rho(x, T = None, h = None)


   .. py:method:: get_cp(x, T = None, h = None)


   .. py:method:: get_gamma(x, T = None, h = None)


   .. py:method:: get_h(x, T = None, h = None)


   .. py:method:: get_H(x, T = None, h = None)


   .. py:method:: get_a(x, T = None, h = None)


   .. py:method:: get_mu(x, T = None, h = None)


   .. py:method:: get_k(x, T = None, h = None)


   .. py:method:: get_Pr(x, T = None, h = None)


   .. py:method:: get_X(x)

      Interpolated mole-fraction dict at position x.



.. py:class:: EquilNozzleInputs

   .. py:attribute:: mech
      :type:  str


   .. py:attribute:: chamber_solution
      :type:  Optional[cantera.Solution]
      :value: None



   .. py:attribute:: chamber_T
      :type:  Optional[float]
      :value: None



   .. py:attribute:: chamber_p
      :type:  Optional[float]
      :value: None



   .. py:attribute:: chamber_X
      :type:  Optional[Union[str, Dict[str, float], numpy.ndarray]]
      :value: None



   .. py:attribute:: A_throat
      :type:  float
      :value: 0.0



   .. py:attribute:: A_exit
      :type:  Optional[float]
      :value: None



   .. py:attribute:: p_amb
      :type:  Optional[float]
      :value: None



   .. py:attribute:: fd_eps
      :type:  float
      :value: 0.0001



.. py:class:: EquilibriumNozzle(inp)

   .. py:attribute:: inp


   .. py:attribute:: gas


   .. py:attribute:: h0


   .. py:attribute:: s0


   .. py:attribute:: throat


   .. py:attribute:: mdot


   .. py:method:: solve_throat()

      Unknowns: (T*, p*) with constraints:
        1) s(T*,p*,eq) = s0
        2) h0 - h(T*,p*,eq) - a_eq(T*,p*)^2 / 2 = 0   [since M=1 ⇒ V = a_eq]



   .. py:method:: solve_exit_at_area_ratio(Ae_over_At)


   .. py:method:: thrust(exit_state)

      1-D thrust (no loss model): F = mdot*V_e + (p_e - p_amb) A_e



