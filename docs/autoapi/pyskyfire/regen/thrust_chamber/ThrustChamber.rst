pyskyfire.regen.thrust_chamber.ThrustChamber
============================================

.. py:class:: pyskyfire.regen.thrust_chamber.ThrustChamber(contour, wall_group, cooling_circuit_group, combustion_transport, optimal_values=None, roughness=1.5e-05, K_factor=0.3, n_nodes=50, h_gas_corr=1.0, h_cold_corr=1.0)

   
   Full thrust-chamber assembly combining geometry, cooling, and combustion models.

   Acts as the top-level container linking the hot-gas contour, wall stack,
   and cooling circuits into a coherent physical representation.

   :Parameters:

       **contour** : :obj:`Contour`
           Hot-gas contour defining inner geometry.

       **wall_group** : :obj:`WallGroup`
           Structural wall stack.

       **cooling_circuit_group** : :obj:`CoolingCircuitGroup`
           Collection of cooling circuits.

       **combustion_transport** : :obj:`object`
           Provides combustion-gas properties and flow variables.

       **optimal_values** : :class:`python:dict`, :obj:`optional`
           Optional dictionary of reference operating conditions.

       **roughness** : :class:`python:float` or :func:`python:callable`, :obj:`optional`
           Effective roughness height of coolant walls [m].

       **K_factor** : :class:`python:float`, :obj:`optional`
           Curvature-loss coefficient.

       **n_nodes** : :class:`python:int`, :obj:`optional`
           Number of discrete axial samples.

       **h_gas_corr, h_cold_corr** : :class:`python:float`, :obj:`optional`
           Empirical correction factors for gas- and coolant-side correlations.

   :Attributes:

       **contour** : :obj:`Contour`
           Geometric shape of the chamber/nozzle.

       **wall_group** : :obj:`WallGroup`
           Walls through which conduction occurs.

       **cooling_circuit_group** : :obj:`CoolingCircuitGroup`
           All defined cooling circuits.

       **combustion_transport** : :obj:`object`
           Hot-gas property model.

       **n_nodes** : :class:`python:int`
           Number of discretization points.

       **K_factor** : :class:`python:float`
           Curvature loss coefficient.

       **h_gas_corr, h_cold_corr** : :class:`python:float`
           Correction multipliers.

       **_roughness** : :class:`python:float` or :func:`python:callable`
           Underlying roughness definition.









   .. seealso::

       
       :obj:`CoolingCircuit`
           ..
       :obj:`WallGroup`
           ..
       :obj:`Contour`
           ..
       
   .. rubric:: Notes

   The `ThrustChamber` automatically initializes derived circuit geometry
   and may call combustion-transport property generation at construction.



   ..
       !! processed by numpydoc !!

   .. py:method:: build_channel_centerlines(mode='sim')

      
      Build centerline splines for each CoolingCircuit.
      For each circuit, use its pre-built x_domain.
      Each circuit is assigned angles in an interleaved fashion.
















      ..
          !! processed by numpydoc !!


   .. py:method:: build_channel_heights()

      
      Compute the channel heights for each cooling circuit along its pre-built x_domain.
      Evaluate the channel height function at each x in the circuit's domain.
















      ..
          !! processed by numpydoc !!


   .. py:method:: build_channel_widths()

      
      Compute the channel widths (in radians) for each cooling circuit.
      Uses each circuit's pre-built x_domain and the new number_of_channels(x)
      function to determine the total active channels at each x position.
















      ..
          !! processed by numpydoc !!


   .. py:method:: build_circuit_x_domain()

      
      Build the x-domain for each cooling circuit by converting its fractional span
      into actual x-values. The sign and ordering of the span determine the coolant flow direction.
      This function uses the overall engine x-range from the contour.
















      ..
          !! processed by numpydoc !!


   .. py:method:: build_t_wall_tot()

      
      Build an array of total wall thicknesses along each circuit's x-domain and
      assign it to the corresponding cooling circuit using set_t_wall_tot.
















      ..
          !! processed by numpydoc !!


   .. py:method:: roughness(x)

      
      Get the channel roughness, at a position, x.
















      ..
          !! processed by numpydoc !!

