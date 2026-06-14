pyskyfire.regen.thrust_chamber.ThrustChamber
============================================

.. py:class:: pyskyfire.regen.thrust_chamber.ThrustChamber(contour, cooling_circuits, combustion_transport, optimal_values=None, K_factor: float = 0.3, n_nodes: int = 50, h_gas_corr: float = 1.0, h_cold_corr: float = 1.0, compute_gas: bool = True, enable_fin: bool = True, film_cooling: FilmCooling | None = None)

   
   Simulation-only thrust-chamber assembly.
















   ..
       !! processed by numpydoc !!

   .. py:method:: _resolve_signed_fraction_to_x(f: float) -> float


   .. py:method:: build_channel_centerlines()


   .. py:method:: build_channel_heights()

      
      Evaluate the per-circuit height law h(x) on each circuit domain.
















      ..
          !! processed by numpydoc !!


   .. py:method:: build_channel_widths()

      
      Compute wedge angle (theta) arrays used by cross-section analytics.

      Priority:
      1) If placement provides channel_width(x), use it.
      2) If placement is 'internal' (non-occluding), use height-as-width (simple radial-stack proxy).
      3) Otherwise, default to an even packing among surface-occluding channels:
         theta(x) = 2π / n_occ(x), where n_occ is provided by the circuit group.















      ..
          !! processed by numpydoc !!


   .. py:method:: build_circuit_x_domain()

      
      Convert each circuit's fractional span into an x-grid whose size
      scales with the fraction of the overall contour length covered.
      Uses at least 3 nodes per circuit.
















      ..
          !! processed by numpydoc !!


   .. py:method:: number_of_channels(x, *, occluding_only=False)

      
      Return the total number of active channels at position `x`.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].

          **occluding_only** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
              If True, count only circuits that occlude the wall surface.



      :Returns:

          :class:`python:int`
              Total number of channels currently active at `x`.











      ..
          !! processed by numpydoc !!

