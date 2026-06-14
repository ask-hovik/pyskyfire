pyskyfire.viz.plot_regen.PlotHeatTransferCoefficient
====================================================

.. py:class:: pyskyfire.viz.plot_regen.PlotHeatTransferCoefficient(*cooling_data_dicts, hot: bool = True, cold: bool = True)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`

   .. autoapi-inheritance-diagram:: pyskyfire.viz.plot_regen.PlotHeatTransferCoefficient
      :parts: 1
      :private-bases:


   
   Plot returned regenerative-cooling heat transfer coefficients.

   Expects dicts with keys:
     - 'x'
     - 'h_hot'   : effective hot-side temperature-based h [W/m²/K]
     - 'h_cold'  : coolant-side h [W/m²/K]
   Optional:
     - 'name'















   ..
       !! processed by numpydoc !!
