pyskyfire.viz.plot_common.PlotEngineNetwork
===========================================

.. py:class:: pyskyfire.viz.plot_common.PlotEngineNetwork(engine_network, title: str | None = 'Engine Network', station_mode: str = 'values', mass_flow_based_arrows: bool = False, edge_length: int = 200, template: str = 'plotly_white', height: int | None = 800, width: int | None = None)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`

   .. autoapi-inheritance-diagram:: pyskyfire.viz.plot_common.PlotEngineNetwork
      :parts: 1
      :private-bases:


   
   Plot an EngineNetwork interactively with Plotly.

   Features:
     • Station label modes: "name" | "values" | "both" | "hidden"
     • Mass-flow–scaled arrow widths (optional)
     • Dashed signal edges
     • Lightweight layout: networkx.spring_layout if available, else circular















   ..
       !! processed by numpydoc !!
