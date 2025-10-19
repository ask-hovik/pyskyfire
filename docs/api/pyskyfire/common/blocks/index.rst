pyskyfire.common.blocks
=======================

.. py:module:: pyskyfire.common.blocks






Module Contents
---------------

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



