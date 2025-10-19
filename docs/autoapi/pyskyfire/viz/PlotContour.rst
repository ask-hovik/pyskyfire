pyskyfire.viz.PlotContour
=========================

.. py:class:: pyskyfire.viz.PlotContour(*contours: Any, show_labels: bool = True, title: str = 'Contour Profiles', template: str = 'plotly_white')

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`

   .. autoapi-inheritance-diagram:: pyskyfire.viz.PlotContour
      :parts: 1
      :private-bases:


   
   Plot bell-nozzle or toroidal-aerospike contours, mirrored about the x-axis.

   Each contour can be either:
     • classic: has .xs and .rs
     • aerospike: has .xs_outer, .rs_outer, .xs_inner, .rs_inner

   If a contour has a .name attribute, it's used for the legend label.















   ..
       !! processed by numpydoc !!
