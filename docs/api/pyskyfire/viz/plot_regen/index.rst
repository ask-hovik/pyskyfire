pyskyfire.viz.plot_regen
========================

.. py:module:: pyskyfire.viz.plot_regen






Module Contents
---------------

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


