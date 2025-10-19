pyskyfire.viz.lookup_tables
===========================

.. py:module:: pyskyfire.viz.lookup_tables




Module Contents
---------------

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


