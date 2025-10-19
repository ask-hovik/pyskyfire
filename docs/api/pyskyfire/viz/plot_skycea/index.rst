pyskyfire.viz.plot_skycea
=========================

.. py:module:: pyskyfire.viz.plot_skycea




Module Contents
---------------

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


