pyskyfire.viz.plot_regen.PlotTemperatureProfile
===============================================

.. py:class:: pyskyfire.viz.plot_regen.PlotTemperatureProfile(results, thrust_chamber, circuit_index: int, x_query: float, n_bl: int = 1000)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`

   .. autoapi-inheritance-diagram:: pyskyfire.viz.plot_regen.PlotTemperatureProfile
      :parts: 1
      :private-bases:


   
   Temperature profile T(y) across gas boundary layer, wall, and coolant
   boundary layer at a given axial location x_query.

   Inputs:
     - results: dict with keys "x", "T" (N×n_layers), "T_static", "dQ_dA", "p_static"
     - thrust_chamber: provides combustion_transport and cooling circuits
     - circuit_index: which cooling circuit to use for coolant props
     - x_query: axial location (m) to sample















   ..
       !! processed by numpydoc !!
