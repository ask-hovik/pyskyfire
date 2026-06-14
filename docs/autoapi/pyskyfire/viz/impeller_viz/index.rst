pyskyfire.viz.impeller_viz
==========================

.. py:module:: pyskyfire.viz.impeller_viz

.. autoapi-nested-parse::

   PyVista visualisation helpers for pyskyfire pump impellers.

   The functions here deliberately consume an ``Impeller`` instance rather than
   recomputing geometry.  This mirrors the thrust-chamber visualisation style: the
   engineering object prepares geometry; the viz module only turns it into meshes
   and plots.

   ..
       !! processed by numpydoc !!


Functions
---------

.. autoapisummary::

   pyskyfire.viz.impeller_viz._structured_surface_from_xyz
   pyskyfire.viz.impeller_viz.export_impeller_surfaces
   pyskyfire.viz.impeller_viz.make_blade_meshes
   pyskyfire.viz.impeller_viz.make_impeller_meshes
   pyskyfire.viz.impeller_viz.make_shroud_meshes
   pyskyfire.viz.impeller_viz.plot_impeller


Module Contents
---------------

.. py:function:: _structured_surface_from_xyz(xyz: numpy.ndarray) -> pyvista.StructuredGrid

   
   Create a PyVista structured surface from ``(ni, nj, 3)`` points.
















   ..
       !! processed by numpydoc !!

.. py:function:: export_impeller_surfaces(impeller: pyskyfire.pump.impeller_new.Impeller, path: str | pathlib.Path, *, binary: bool = True) -> pathlib.Path

   
   Export current visualisation surfaces as a single PolyData mesh.

   This is a visual/CAD-reference mesh, not a watertight manufacturing model.















   ..
       !! processed by numpydoc !!

.. py:function:: make_blade_meshes(impeller: pyskyfire.pump.impeller_new.Impeller) -> list[pyvista.StructuredGrid]

   
   Return one camber-surface mesh per blade.
















   ..
       !! processed by numpydoc !!

.. py:function:: make_impeller_meshes(impeller: pyskyfire.pump.impeller_new.Impeller, *, n_theta: int = 160) -> dict[str, object]

   
   Build all currently available impeller visualisation meshes.
















   ..
       !! processed by numpydoc !!

.. py:function:: make_shroud_meshes(impeller: pyskyfire.pump.impeller_new.Impeller, *, n_theta: int = 160) -> tuple[pyvista.StructuredGrid, pyvista.StructuredGrid]

   
   Create hub and shroud surfaces by revolving their meridional curves.
















   ..
       !! processed by numpydoc !!

.. py:function:: plot_impeller(impeller: pyskyfire.pump.impeller_new.Impeller, *, show_edges: bool = True, show_shrouds: bool = True, show_blades: bool = True, blade_opacity: float = 0.75, shroud_opacity: float = 0.22, window_size: tuple[int, int] = (1300, 900), screenshot: str | pathlib.Path | None = None, return_plotter: bool = False) -> pyvista.Plotter | None

   
   Plot a preliminary 3D impeller using PyVista.


   :Parameters:

       **impeller:**
           ``pyskyfire.pump.impeller.Impeller`` instance.

       **screenshot:**
           Optional path.  If supplied, an off-screen screenshot is written there.

       **return_plotter:**
           If true, return the ``pv.Plotter`` after adding meshes.  Useful for GUI
           embedding or adding custom annotations.














   ..
       !! processed by numpydoc !!

