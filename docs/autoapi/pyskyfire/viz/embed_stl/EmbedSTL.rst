pyskyfire.viz.embed_stl.EmbedSTL
================================

.. py:class:: pyskyfire.viz.embed_stl.EmbedSTL(stl_path: str, *, color: str = '#9cc4ff', opacity: float = 1.0, show_wireframe: bool = False, template: str = 'plotly_white', process: bool = False)

   Bases: :py:obj:`pyskyfire.viz.core.PlotBase`

   .. autoapi-inheritance-diagram:: pyskyfire.viz.embed_stl.EmbedSTL
      :parts: 1
      :private-bases:


   
   Thin wrapper that turns an STL on disk into a Plotly Mesh3d figure.
   Keeps the figure open for further Plotly commands via PlotBase.
















   ..
       !! processed by numpydoc !!

   .. py:method:: _add_mesh_trace(*, color: str, opacity: float)


   .. py:method:: add_section_plane(plane_point: Tuple[float, float, float], plane_normal: Tuple[float, float, float], *, name: str = 'Section', width: float = 3.0) -> EmbedSTL

      
      Slice the mesh with a plane and overlay polyline(s) of the intersection.
















      ..
          !! processed by numpydoc !!


   .. py:method:: add_wireframe(*, name: str = 'Wireframe', width: float = 1.0) -> EmbedSTL

      
      Overlay a wireframe by drawing all triangle edges as line segments.
















      ..
          !! processed by numpydoc !!


   .. py:method:: recenter() -> EmbedSTL

      
      Translate so the model centroid is near the origin.
















      ..
          !! processed by numpydoc !!


   .. py:property:: faces
      :type: numpy.ndarray



   .. py:property:: vertices
      :type: numpy.ndarray


