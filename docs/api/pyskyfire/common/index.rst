pyskyfire.common
================

.. py:module:: pyskyfire.common


Submodules
----------

.. toctree::
   :maxdepth: 1

   /api/pyskyfire/common/blocks/index
   /api/pyskyfire/common/engine_network/index
   /api/pyskyfire/common/fluids/index
   /api/pyskyfire/common/results/index
   /api/pyskyfire/common/solids/index








Package Contents
----------------

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



.. py:class:: Station

   .. py:attribute:: p
      :type:  float


   .. py:attribute:: T
      :type:  float


   .. py:attribute:: mdot
      :type:  float


.. py:class:: BoundaryConditions(T_coolant_in, p_coolant_in, mdot_coolant)

   Object for storing boundary conditions for the solver.

   :param T_coolant_in: Static(?) temperature of coolant at cooling channel inlet (K)
   :type T_coolant_in: float
   :param p_coolant_in: Static(?) pressure of coolant at cooling channel inlet (Pa)
   :type p_coolant_in: float
   :param mdot_coolant: mass flow rate of coolant through cooling channel inlet (kg/s)
   :type mdot_coolant: float


   .. py:attribute:: T_coolant_in


   .. py:attribute:: p_coolant_in


   .. py:attribute:: mdot_coolant


.. py:function:: steady_heating_analysis(thrust_chamber, boundary_conditions, n_nodes=100, circuit_index=0, solver='newton', output=True)

   Run the steady heating analysis.

   :param engine: Engine object with geometry, boundary conditions, and physics.
   :param n_nodes: Number of nodes (for Newton) or resolution for post-processing (for Radau).
   :param solver: String, currently only "newton" available

   :returns: A dictionary with simulation results.


.. py:data:: g0
   :value: 9.81


.. py:class:: FluidBlock(medium)

   Bases: :py:obj:`abc.ABC`


   Helper class that provides a standard way to create an ABC using
   inheritance.

   Base constructor for fluid-processing blocks.


   .. py:attribute:: station_inputs
      :type:  list[str]


   .. py:attribute:: station_outputs
      :type:  list[str]


   .. py:attribute:: signal_inputs
      :type:  list[str]


   .. py:attribute:: signal_outputs
      :type:  list[str]


   .. py:attribute:: medium


   .. py:method:: compute(stations_in, signals_in)
      :abstractmethod:


      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



.. py:class:: SignalBlock

   Bases: :py:obj:`abc.ABC`


   Helper class that provides a standard way to create an ABC using
   inheritance.


   .. py:attribute:: station_inputs
      :type:  list[str]


   .. py:attribute:: station_outputs
      :type:  list[str]


   .. py:attribute:: signal_inputs
      :type:  list[str]


   .. py:attribute:: signal_outputs
      :type:  list[str]


   .. py:method:: compute(stations_in, signals_in)
      :abstractmethod:


      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



.. py:class:: PumpBlock(name, st_in, st_out, overcome, p_base, input_p, load_fraction, n, eta, medium)

   Bases: :py:obj:`FluidBlock`


   Encapsulates pump sizing and performance logic.

   params is a dictionary containing keys like:
     - n_fu (pump rotational speed)
     - eta_pump_fu (pump efficiency)
     - rho_fu_tank (liquid density)
     ...


   .. py:attribute:: overcome


   .. py:attribute:: load_fraction


   .. py:attribute:: p_base


   .. py:attribute:: name


   .. py:attribute:: n


   .. py:attribute:: eta


   .. py:attribute:: input_p


   .. py:attribute:: st_in


   .. py:attribute:: st_out


   .. py:attribute:: station_inputs


   .. py:attribute:: station_outputs


   .. py:attribute:: signal_inputs
      :value: []



   .. py:attribute:: signal_outputs


   .. py:method:: compute(stations, signals)

      Given pump inlet/outlet conditions, compute the required power,
      outlet temperature, mass, etc.

      Returns a dict with keys:
        'P_pump' (pump power, W)
        'T_pump_out' (pump outlet temperature, K)




   .. py:attribute:: medium


   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



.. py:class:: RegenBlock(name, st_in, st_out, circuit_index, thrust_chamber, medium)

   Bases: :py:obj:`FluidBlock`


   Regenerative-cooling segment (single circuit).

   :param name: Unique tag; used to form the scalar Δp key.
   :type name: str
   :param st_in: Station key for coolant entering the circuit.
   :type st_in: str
   :param st_out: Station key written back to the network (exit of circuit).
   :type st_out: str
   :param circuit_index: Which cooling circuit of the engine model this block represents.
   :type circuit_index: int
   :param engine: Your pre-built rocket engine object (geometry + prop data).
   :type engine: Engine
   :param Base constructor for fluid-processing blocks.:


   .. py:attribute:: name


   .. py:attribute:: st_in


   .. py:attribute:: st_out


   .. py:attribute:: circuit_index


   .. py:attribute:: thrust_chamber


   .. py:attribute:: station_inputs


   .. py:attribute:: station_outputs


   .. py:attribute:: signal_inputs
      :value: []



   .. py:attribute:: dp_key
      :value: 'dp_Uninferable'



   .. py:attribute:: signal_outputs


   .. py:method:: compute(stations, signals)

      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



   .. py:attribute:: medium


.. py:class:: TurbineBlock(name, st_in, st_out, P_req_key, eta, medium)

   Bases: :py:obj:`FluidBlock`


   Single-stage impulse/expander turbine.
   Produces an outlet Station and its pressure drop; consumes the
   shaft-power demand that the TransmissionBlock summed.

   Base constructor for fluid-processing blocks.


   .. py:attribute:: name


   .. py:attribute:: st_in


   .. py:attribute:: st_out


   .. py:attribute:: P_req_key


   .. py:attribute:: eta


   .. py:attribute:: station_inputs


   .. py:attribute:: station_outputs


   .. py:attribute:: signal_inputs


   .. py:attribute:: dp_key
      :value: 'dp_Uninferable'



   .. py:attribute:: signal_outputs


   .. py:method:: compute(stations, signals)

      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:attribute:: medium


   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



.. py:class:: SimpleDuctBlock(name, st_in, st_out, eta, medium)

   Bases: :py:obj:`FluidBlock`


   Applies a fixed efficiency η to simulate a homogeneous duct loss:
       p_out = η · p_in

   Base constructor for fluid-processing blocks.


   .. py:attribute:: name


   .. py:attribute:: st_in


   .. py:attribute:: st_out


   .. py:attribute:: eta


   .. py:attribute:: station_inputs


   .. py:attribute:: station_outputs


   .. py:attribute:: signal_inputs
      :value: []



   .. py:attribute:: dp_key
      :value: 'dp_Uninferable'



   .. py:attribute:: signal_outputs


   .. py:method:: compute(stations, signals)

      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:attribute:: medium


   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



.. py:class:: FractionalMassFlowLossBlock(*, name, st_in, st_out, fraction, medium)

   Bases: :py:obj:`FluidBlock`


   Peel off a fixed fraction of mdot from st_in and inject it into st_out
   *without* showing the underlying split / merge to normal users.

   Visible metadata:
       station_inputs  = [st_in, st_out]
       station_outputs = [st_in, st_out]
   Internally:
       ┌─ MassFlowSplitterBlock ─┐
       │  st_in  ──► main + bleed│
       └─────────────────────────┘
                        │
                        ▼
       ┌─ MassFlowMergerBlock  ─┐
       │  st_out + bleed ──► out│
       └─────────────────────────┘

   Base constructor for fluid-processing blocks.


   .. py:attribute:: name


   .. py:attribute:: st_in


   .. py:attribute:: st_out


   .. py:attribute:: fraction


   .. py:attribute:: station_inputs


   .. py:attribute:: station_outputs


   .. py:attribute:: signal_inputs
      :value: []



   .. py:attribute:: signal_outputs
      :value: []



   .. py:method:: compute(stations, signals)

      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:attribute:: medium


   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



.. py:class:: OrificePlateBlock(*, name, st_in, st_out, Cd, A, medium)

   Bases: :py:obj:`FluidBlock`


   Thin-plate orifice with a *permanent* pressure loss.

   The usual incompressible relation is applied

       ṁ = C_d · A · √(2 ρ Δp)   ⟹   Δp = (ṁ / (C_d A))² / (2 ρ)

   :param name: Unique tag – shows up in the network’s scalar list as ``dp_<name>``.
   :type name: str
   :param st_in: Keys of the inlet and outlet *Station* objects.
   :type st_in: str
   :param st_out: Keys of the inlet and outlet *Station* objects.
   :type st_out: str
   :param Cd: Discharge coefficient (0 < Cₙ ≤ 1).  Typical sharp-edged orifice ≈ 0.6.
   :type Cd: float
   :param A: Flow area of the hole [m²].
   :type A: float
   :param medium: CoolProp fluid identifier – must match the other stations.
   :type medium: str
   :param Base constructor for fluid-processing blocks.:


   .. py:attribute:: name


   .. py:attribute:: st_in


   .. py:attribute:: st_out


   .. py:attribute:: Cd


   .. py:attribute:: A


   .. py:attribute:: station_inputs


   .. py:attribute:: station_outputs


   .. py:attribute:: signal_inputs
      :value: []



   .. py:attribute:: dp_key
      :value: 'dp_Uninferable'



   .. py:attribute:: signal_outputs


   .. py:method:: compute(stations, signals)

      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:attribute:: medium


   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



.. py:class:: JunctionBlock(*, name, st_in, st_out, medium, equilibrate = None, block_lookup = None, tol = 1e-06)

   Bases: :py:obj:`FluidBlock`


   A generalised mass‑flow junction.

   The block can *merge* an arbitrary number of inlet streams and *split*
   the combined flow into an arbitrary number of outlet streams **without**
   hard‑coded mass‑flow ratios.  Instead, a small one‑dimensional root‑finder
   allocates the branch mass‑flows such that user‑supplied pressure
   constraints (or pressure equalisation between the branches) are met.

   :param name: A human‑readable unique tag – used only for bookkeeping.
   :param st_in: List of input *Station* keys that feed the junction.
   :param st_out: List of output *Station* keys produced by the junction.  Each element
                  must also appear as a key in *equilibrate*.
   :param medium: CoolProp fluid name – *all* inlets are assumed to be of identical
                  composition.
   :param equilibrate: Defines how the total mass‑flow is distributed among the outlets.  It
                       is a list of one‑element dictionaries so that *order* is preserved::

                           [
                               {"fu_leg": ["duct_fu", "turbine"]},
                               {"by_leg": "diff"}
                           ]

                       • If the *value* is ``"diff"`` the branch simply receives whatever is
                         left after the other branches have been solved.
                       • Otherwise, the value must be a *list of block names* that constitute
                         the hydraulic path of that branch **in the EngineNetwork’s block
                         order**.  During *compute* the junction varies the branch mass‑flow
                         until the exit‑pressure of the final block matches the *current*
                         network pressure of its outlet station.

                       If *two* branches meet again downstream (→ equal outlet station), pass
                       two dictionaries, both with a *list* of blocks.  The solver will then
                       iterate to make the *exit pressures equal* between the two branches.
   :param block_lookup: Mapping *block‑name → block‑instance* so that the junction can call
                        ``block.compute`` internally while it searches for a pressure match.
                        The simplest way is to build ``{blk.name: blk for blk in blocks}`` once
                        when assembling the network and pass it to every *JunctionBlock*.
   :param tol: Relative tolerance used by the in‑house bisection routine.
   :param Base constructor for fluid-processing blocks.:


   .. py:attribute:: name


   .. py:attribute:: station_inputs


   .. py:attribute:: station_outputs


   .. py:attribute:: signal_inputs
      :type:  List[str]
      :value: []



   .. py:attribute:: signal_outputs
      :type:  List[str]
      :value: []



   .. py:method:: compute(stations, signals)

      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



   .. py:attribute:: medium


.. py:class:: MassFlowSplitterBlock(name, st_in, st_out, medium, fractions = None, frac_keys = None)

   Bases: :py:obj:`FluidBlock`


   Split one incoming stream into N outgoing branches.

   :param name: Unique tag.  Used only for bookkeeping.
   :type name: str
   :param st_in: Key of the inlet Station.
   :type st_in: str
   :param st_outs: Keys of the outlet Stations (length = N).
   :type st_outs: list[str]
   :param fractions: Fixed mass-flow fractions that sum to 1.0, OR
                     None if you want to supply them at run-time
                     through scalar signals (see *dynamic split* below).
   :type fractions: list[float] | None
   :param frac_keys: If you choose a dynamic split, give one signal key
                     per outlet.  The network will read those each pass.
   :type frac_keys: list[str] | None
   :param Base constructor for fluid-processing blocks.:


   .. py:attribute:: name


   .. py:attribute:: st_in


   .. py:attribute:: st_out


   .. py:attribute:: station_inputs


   .. py:attribute:: station_outputs


   .. py:attribute:: signal_inputs
      :value: []



   .. py:attribute:: signal_outputs
      :value: []



   .. py:method:: compute(stations, signals)

      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:attribute:: medium


   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



.. py:class:: MassFlowMergerBlock(name, st_in, st_out, medium)

   Bases: :py:obj:`FluidBlock`


   Combine multiple inlet streams into one outlet.

   NOTE:  Assumes identical fluid species.

   Base constructor for fluid-processing blocks.


   .. py:attribute:: name


   .. py:attribute:: st_in


   .. py:attribute:: st_out


   .. py:attribute:: station_inputs


   .. py:attribute:: station_outputs


   .. py:attribute:: signal_inputs
      :value: []



   .. py:attribute:: signal_outputs
      :value: []



   .. py:method:: compute(stations, signals)

      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:attribute:: medium


   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



.. py:class:: TransmissionBlock(name, sink_keys, out_key='P_required')

   Bases: :py:obj:`SignalBlock`


   Sums shaft-power keys → one P_required key for the turbine.


   .. py:attribute:: station_inputs
      :value: []



   .. py:attribute:: station_outputs
      :value: []



   .. py:attribute:: signal_inputs


   .. py:attribute:: signal_outputs


   .. py:attribute:: sink_keys


   .. py:attribute:: out_key
      :value: 'P_required'



   .. py:attribute:: name


   .. py:method:: compute(st, sg)

      Return two dicts:
          • updated/created Station objects
          • updated/created scalar signals



   .. py:method:: post_process(stations, signals)

      Called once after convergence; override when you have extras.



.. py:class:: Results

   Bases: :py:obj:`collections.abc.MutableMapping`


   A mapping you can also dot-access, with save/load built in.


   .. py:attribute:: metadata


   .. py:method:: save(path)


   .. py:method:: load(path)
      :classmethod:



   .. py:method:: add(name, obj)

      Explicitly register something new.



   .. py:method:: list_items()


   .. py:method:: pop(key, default=__marker)

      D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
      If key is not found, d is returned if given, otherwise KeyError is raised.



   .. py:method:: popitem()

      D.popitem() -> (k, v), remove and return some (key, value) pair
      as a 2-tuple; but raise KeyError if D is empty.



   .. py:method:: clear()

      D.clear() -> None.  Remove all items from D.



   .. py:method:: update(other=(), /, **kwds)

      D.update([E, ]**F) -> None.  Update D from mapping/iterable E and F.
      If E present and has a .keys() method, does:     for k in E.keys(): D[k] = E[k]
      If E present and lacks .keys() method, does:     for (k, v) in E: D[k] = v
      In either case, this is followed by: for k, v in F.items(): D[k] = v



   .. py:method:: setdefault(key, default=None)

      D.setdefault(k[,d]) -> D.get(k,d), also set D[k]=d if k not in D



   .. py:method:: get(key, default=None)

      D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.



   .. py:method:: keys()

      D.keys() -> a set-like object providing a view on D's keys



   .. py:method:: items()

      D.items() -> a set-like object providing a view on D's items



   .. py:method:: values()

      D.values() -> an object providing a view on D's values



.. py:class:: Station

   .. py:attribute:: p
      :type:  float


   .. py:attribute:: T
      :type:  float


   .. py:attribute:: mdot
      :type:  float


.. py:class:: EngineNetwork(stations, signals, blocks)

   Fixed-point network that stores *two* dictionaries:
       • stations : thermodynamic states        (Station objects)
       • signals  : scalar data (Δp, powers, targets, flags …)
   Every Block advertises four lists:
       station_inputs,  station_outputs,
       signal_inputs,   signal_outputs
   and its compute() returns (stations_out, signals_out).


   .. py:attribute:: stations


   .. py:attribute:: signals


   .. py:attribute:: blocks


   .. py:attribute:: residuals
      :type:  List[float]
      :value: []



   .. py:method:: update()

      Executes one full sweep over all blocks.
      Returns the maximum relative change (residual) seen.



   .. py:method:: run_fixed_point(tol = 1e-06, max_iter = 100)

      Iterate until the largest relative change among *all* advertised
      outputs is below *tol* or until *max_iter* is reached.



