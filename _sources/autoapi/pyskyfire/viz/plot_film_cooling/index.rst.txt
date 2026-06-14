pyskyfire.viz.plot_film_cooling
===============================

.. py:module:: pyskyfire.viz.plot_film_cooling

.. autoapi-nested-parse::

   pyskyfire/viz/plot_film_cooling.py

   PlotBase subclasses for visualising contour-based Grisson film cooling results.

   These plots are built around the current result containers returned by:
       GrissonFilmCoolingModel.solve(...)

   Returned data
   -------------
   LiquidFilmResults:
       x
       Gamma
       m_vap
       h_conv
       Q_conv
       Q_rad
       T_film
       x_dryout
       Gamma_dryout

   GaseousFilmResults:
       x
       Mbl
       T_aw
       T_w
       h_conv
       Q_rad

   ..
       !! processed by numpydoc !!


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/viz/plot_film_cooling/PlotFilmBoundaryLayer
   /autoapi/pyskyfire/viz/plot_film_cooling/PlotFilmDryout
   /autoapi/pyskyfire/viz/plot_film_cooling/PlotFilmEffectiveWallHTC
   /autoapi/pyskyfire/viz/plot_film_cooling/PlotFilmEvaporationRate
   /autoapi/pyskyfire/viz/plot_film_cooling/PlotFilmHeatFlux
   /autoapi/pyskyfire/viz/plot_film_cooling/PlotFilmHeatTransferCoefficient
   /autoapi/pyskyfire/viz/plot_film_cooling/PlotFilmRadiativeFraction
   /autoapi/pyskyfire/viz/plot_film_cooling/PlotFilmWallTemperature
   /autoapi/pyskyfire/viz/plot_film_cooling/PlotGaseousFilmRadiation
   /autoapi/pyskyfire/viz/plot_film_cooling/PlotGaseousFilmTemperatureGap

.. autoapisummary::

   pyskyfire.viz.plot_film_cooling.PlotFilmBoundaryLayer
   pyskyfire.viz.plot_film_cooling.PlotFilmDryout
   pyskyfire.viz.plot_film_cooling.PlotFilmEffectiveWallHTC
   pyskyfire.viz.plot_film_cooling.PlotFilmEvaporationRate
   pyskyfire.viz.plot_film_cooling.PlotFilmHeatFlux
   pyskyfire.viz.plot_film_cooling.PlotFilmHeatTransferCoefficient
   pyskyfire.viz.plot_film_cooling.PlotFilmRadiativeFraction
   pyskyfire.viz.plot_film_cooling.PlotFilmWallTemperature
   pyskyfire.viz.plot_film_cooling.PlotGaseousFilmRadiation
   pyskyfire.viz.plot_film_cooling.PlotGaseousFilmTemperatureGap


Functions
---------

.. autoapisummary::

   pyskyfire.viz.plot_film_cooling._maybe_array


Module Contents
---------------

.. py:function:: _maybe_array(lst) -> numpy.ndarray

