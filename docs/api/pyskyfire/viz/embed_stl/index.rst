pyskyfire.viz.embed_stl
=======================

.. py:module:: pyskyfire.viz.embed_stl




Module Contents
---------------

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


