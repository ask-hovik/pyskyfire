pyskyfire.viz.PlotMoodyDiagram
==============================

.. py:class:: pyskyfire.viz.PlotMoodyDiagram(rel_rough_list: Optional[Sequence[float]] = None, Re_min: float = 700.0, Re_max: float = 100000000.0, n_pts: int = 400, template: str = 'plotly_white')

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`

   .. autoapi-inheritance-diagram:: pyskyfire.viz.PlotMoodyDiagram
      :parts: 1
      :private-bases:


   
   Moody-style plot of Darcy friction factor vs Reynolds number, using
   the correlations implemented in pyskyfire.physics.f_darcy.

   Each curve corresponds to a chosen relative roughness ε/D.















   ..
       !! processed by numpydoc !!
