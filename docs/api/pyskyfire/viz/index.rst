pyskyfire.viz
=============

.. py:module:: pyskyfire.viz


Submodules
----------

.. toctree::
   :maxdepth: 1

   /api/pyskyfire/viz/core/index
   /api/pyskyfire/viz/embed_stl/index
   /api/pyskyfire/viz/export_3d/index
   /api/pyskyfire/viz/lookup_tables/index
   /api/pyskyfire/viz/plot_common/index
   /api/pyskyfire/viz/plot_regen/index
   /api/pyskyfire/viz/plot_skycea/index
   /api/pyskyfire/viz/report/index








Package Contents
----------------

.. py:class:: PlotBase(fig = None)

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


.. py:class:: PlotContour(*contours, show_labels = True, title = 'Contour Profiles', template = 'plotly_white')

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Plot bell-nozzle or toroidal-aerospike contours, mirrored about the x-axis.

   Each contour can be either:
     • classic: has .xs and .rs
     • aerospike: has .xs_outer, .rs_outer, .xs_inner, .rs_inner

   If a contour has a .name attribute, it's used for the legend label.


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


.. py:class:: PlotWallTemperature(*cooling_data_dicts, plot_hot = True, plot_interfaces = False, plot_coolant_wall = False, template = 'plotly_white')

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Build a wall-temperature plot from one or more cooling-data dicts.

   Each input dict must contain:
     - 'x': 1D array of axial positions
     - 'T': 2D array (len(x) × n_layers)

   Optional per-dataset label:
     - 'name': str  (used as legendgroup; defaults to 'Set {i+1}')


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


.. py:class:: PlotCoolantTemperature(*cooling_data_dicts)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Expects dicts with keys: 'x', 'T_static', optional 'name'.


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


.. py:class:: PlotCoolantPressure(*cooling_data_dicts, static = True, stagnation = True)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Expects dicts with keys: 'x', 'p_static', 'p_stagnation' (either/both may exist), optional 'name'.
   Use flags to include static and/or stagnation traces.


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


.. py:class:: PlotHeatFlux(*cooling_data_dicts)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Expects dicts with keys: 'x', 'dQ_dA', optional 'name'.


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


.. py:class:: PlotVelocity(*cooling_data_dicts)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Expects dicts with keys: 'x', 'velocity', optional 'name'.


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


.. py:data:: IndexLike

.. py:class:: PlotdAdxThermalHotGas(thrust_chamber, circuit_index = None)

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


.. py:class:: PlotdAdxThermalCoolant(thrust_chamber, circuit_index = None)

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


.. py:class:: PlotCoolantArea(thrust_chamber, circuit_index = None)

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


.. py:class:: PlotdAdxCoolantArea(thrust_chamber, circuit_index = None)

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


.. py:class:: PlotHydraulicDiameter(thrust_chamber, circuit_index = None)

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


.. py:class:: PlotRadiusOfCurvature(thrust_chamber, circuit_index = None)

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


.. py:class:: PlotTemperatureProfile(results, thrust_chamber, circuit_index, x_query, n_bl = 1000)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Temperature profile T(y) across gas boundary layer, wall, and coolant
   boundary layer at a given axial location x_query.

   Inputs:
     - results: dict with keys "x", "T" (N×n_layers), "T_static", "dQ_dA", "p_static"
     - thrust_chamber: provides combustion_transport and cooling circuits
     - circuit_index: which cooling circuit to use for coolant props
     - x_query: axial location (m) to sample


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


.. py:function:: make_engine_gmsh(thrust_chamber, filename='engine', display_channels=None)

.. py:function:: format_value(value, precision = 6)

   Human-friendly string for table cells.


.. py:function:: dict_to_table_html(data, *, col_key = 'Key', col_val = 'Value', caption = None, precision = 6)

   Return HTML for a simple 2-column table of a dict.


.. py:class:: Tab

   .. py:attribute:: title
      :type:  str


   .. py:attribute:: blocks
      :type:  List[_Block]
      :value: []



   .. py:method:: add_text(text)


   .. py:method:: add_figure(fig, caption = None)


   .. py:method:: add_table(data, *, caption = None, key_title = 'Key', value_title = 'Value', precision = 6)

      Append a key/value HTML table made from a dict.



   .. py:method:: add_raw_html(raw_html)


   .. py:method:: add_image(path, *, alt = '', caption = None, style = 'max-width:100%;height:auto;')


   .. py:method:: add_svg(svg, *, caption = None)


.. py:class:: Report(title = 'pyskyfire Report')

   .. py:attribute:: title
      :value: 'pyskyfire Report'



   .. py:attribute:: tabs
      :type:  List[Tab]
      :value: []



   .. py:method:: add_tab(title)


   .. py:method:: save_html(path)


.. py:class:: PlotBase(fig = None)

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


.. py:class:: PlotBase(fig = None)

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


.. py:class:: EmbedSTL(stl_path, *, color = '#9cc4ff', opacity = 1.0, show_wireframe = False, template = 'plotly_white', process = False)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Thin wrapper that turns an STL on disk into a Plotly Mesh3d figure.
   Keeps the figure open for further Plotly commands via PlotBase.

   :param stl_path: Path to the STL file on disk.
   :type stl_path: str
   :param color: Mesh color (any Plotly color string).
   :type color: str
   :param opacity: Mesh opacity in [0, 1].
   :type opacity: float
   :param show_wireframe: If True, also overlays an approximate wireframe using triangle edges.
   :type show_wireframe: bool
   :param template: Plotly layout template.
   :type template: str
   :param process: Let trimesh 'process' the mesh (repairs, merges). Defaults off to preserve input.
   :type process: bool


   .. py:property:: vertices
      :type: numpy.ndarray



   .. py:property:: faces
      :type: numpy.ndarray



   .. py:method:: add_wireframe(*, name = 'Wireframe', width = 1.0)

      Overlay a wireframe by drawing all triangle edges as line segments.



   .. py:method:: add_section_plane(plane_point, plane_normal, *, name = 'Section', width = 3.0)

      Slice the mesh with a plane and overlay polyline(s) of the intersection.



   .. py:method:: recenter()

      Translate so the model centroid is near the origin.



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


.. py:class:: PlotBase(fig = None)

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


.. py:function:: f_darcy(ReDh, Dh, x, roughness)

.. py:class:: PlotThetaVsEpsilon(template = 'plotly_white')

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


.. py:class:: PlotMoodyDiagram(rel_rough_list = None, Re_min = 700.0, Re_max = 100000000.0, n_pts = 400, template = 'plotly_white')

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Moody-style plot of Darcy friction factor vs Reynolds number, using
   the correlations implemented in pyskyfire.physics.f_darcy.

   Each curve corresponds to a chosen relative roughness ε/D.


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


.. py:class:: PlotBase(fig = None)

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


.. py:class:: PlotTransportProperty(*ats, prop, template = 'plotly_white')

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`


   Plot a single transport-property map (equilibrium column vs x) for one or more
   Aerothermodynamics objects.

   Each object must have:
     - .x_nodes (built by compute_aerothermodynamics)
     - .<prop>_map with shape (Nx, Nt); we use column 0 (equilibrium).


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


.. py:class:: PlotBase(fig = None)

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


