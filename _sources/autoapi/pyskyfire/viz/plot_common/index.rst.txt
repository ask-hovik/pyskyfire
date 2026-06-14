pyskyfire.viz.plot_common
=========================

.. py:module:: pyskyfire.viz.plot_common


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/viz/plot_common/PlotEngineNetwork
   /autoapi/pyskyfire/viz/plot_common/PlotPTDiagram
   /autoapi/pyskyfire/viz/plot_common/PlotResidualHistory
   /autoapi/pyskyfire/viz/plot_common/PlotStationProperty

.. autoapisummary::

   pyskyfire.viz.plot_common.PlotEngineNetwork
   pyskyfire.viz.plot_common.PlotPTDiagram
   pyskyfire.viz.plot_common.PlotResidualHistory
   pyskyfire.viz.plot_common.PlotStationProperty


Functions
---------

.. autoapisummary::

   pyskyfire.viz.plot_common.render_engine_network


Module Contents
---------------

.. py:function:: render_engine_network(engine_network, *, height: str = '900px', width: str = '100%', edge_length: int = 200, physics_settings: Dict | str | None = 'default', mass_flow_based_arrows: bool = False, station_mode: Literal['name', 'values', 'both', 'hidden'] = 'values') -> str

   
   Build a PyVis network and return an <iframe srcdoc="..."> HTML snippet
   suited for Report.add_raw_html(...). Produces a single-file embed.
















   ..
       !! processed by numpydoc !!

