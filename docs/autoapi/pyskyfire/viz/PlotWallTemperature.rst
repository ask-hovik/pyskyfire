pyskyfire.viz.PlotWallTemperature
=================================

.. py:class:: pyskyfire.viz.PlotWallTemperature(*cooling_data_dicts, plot_hot: bool = True, plot_interfaces: bool = False, plot_coolant_wall: bool = False, template: str = 'plotly_white')

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`

   .. autoapi-inheritance-diagram:: pyskyfire.viz.PlotWallTemperature
      :parts: 1
      :private-bases:


   
   Build a wall-temperature plot from one or more cooling-data dicts.

   Each input dict must contain:
     - 'x': 1D array of axial positions
     - 'T': 2D array (len(x) × n_layers)

   Optional per-dataset label:
     - 'name': str  (used as legendgroup; defaults to 'Set {i+1}')















   ..
       !! processed by numpydoc !!
