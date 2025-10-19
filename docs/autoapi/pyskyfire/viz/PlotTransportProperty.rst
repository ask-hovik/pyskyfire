pyskyfire.viz.PlotTransportProperty
===================================

.. py:class:: pyskyfire.viz.PlotTransportProperty(*ats, prop: str, template: str = 'plotly_white')

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`

   .. autoapi-inheritance-diagram:: pyskyfire.viz.PlotTransportProperty
      :parts: 1
      :private-bases:


   
   Plot a single transport-property map (equilibrium column vs x) for one or more
   Aerothermodynamics objects.

   Each object must have:
     - .x_nodes (built by compute_aerothermodynamics)
     - .<prop>_map with shape (Nx, Nt); we use column 0 (equilibrium).















   ..
       !! processed by numpydoc !!
