pyskyfire.viz.plot_common
=========================

.. py:module:: pyskyfire.viz.plot_common






Module Contents
---------------

.. py:class:: PlotResidualHistory(residuals, name='Residual', title='Residual Convergence History', template='plotly_white')

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Common helpers + Plotly passthrough.


   .. py:attribute:: fig


   .. py:attribute:: layout


   .. py:attribute:: traces


   .. py:attribute:: xaxis


   .. py:attribute:: yaxis


   .. py:attribute:: add


   .. py:method:: template(name)


   .. py:method:: config(**cfg)


   .. py:method:: show()


   .. py:method:: save_html(path)


   .. py:method:: save_png(path, scale = 2)


.. py:class:: PlotStationProperty(station_dicts, station_list, property_name, labels=None, title=True, ylabel=None, template='plotly_white', ylim=None)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Common helpers + Plotly passthrough.


   .. py:attribute:: fig


   .. py:attribute:: layout


   .. py:attribute:: traces


   .. py:attribute:: xaxis


   .. py:attribute:: yaxis


   .. py:attribute:: add


   .. py:method:: template(name)


   .. py:method:: config(**cfg)


   .. py:method:: show()


   .. py:method:: save_html(path)


   .. py:method:: save_png(path, scale = 2)


.. py:class:: PlotPTDiagram(station_dicts, station_list, fluid_name=None, title=None, sat_points=200, labels=None, template='plotly_white', annotate_ha=None, annotate_va=None, scale='log')

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Common helpers + Plotly passthrough.


   .. py:attribute:: fig


   .. py:attribute:: layout


   .. py:attribute:: traces


   .. py:attribute:: xaxis


   .. py:attribute:: yaxis


   .. py:attribute:: add


   .. py:method:: template(name)


   .. py:method:: config(**cfg)


   .. py:method:: show()


   .. py:method:: save_html(path)


   .. py:method:: save_png(path, scale = 2)


.. py:class:: PlotEngineNetwork(engine_network, title = 'Engine Network', station_mode = 'values', mass_flow_based_arrows = False, edge_length = 200, template = 'plotly_white', height = 800, width = None)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Plot an EngineNetwork interactively with Plotly.

   Features:
     • Station label modes: "name" | "values" | "both" | "hidden"
     • Mass-flow–scaled arrow widths (optional)
     • Dashed signal edges
     • Lightweight layout: networkx.spring_layout if available, else circular


   .. py:attribute:: fig


   .. py:attribute:: layout


   .. py:attribute:: traces


   .. py:attribute:: xaxis


   .. py:attribute:: yaxis


   .. py:attribute:: add


   .. py:method:: template(name)


   .. py:method:: config(**cfg)


   .. py:method:: show()


   .. py:method:: save_html(path)


   .. py:method:: save_png(path, scale = 2)


.. py:function:: render_engine_network(engine_network, *, height = '900px', width = '100%', edge_length = 200, physics_settings = 'default', mass_flow_based_arrows = False, station_mode = 'values')

   Build a PyVis network and return an <iframe srcdoc="..."> HTML snippet
   suited for Report.add_raw_html(...). Produces a single-file embed.


