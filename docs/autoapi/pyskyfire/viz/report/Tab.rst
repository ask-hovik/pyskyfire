pyskyfire.viz.report.Tab
========================

.. py:class:: pyskyfire.viz.report.Tab

   .. py:method:: add_figure(fig: plotly.graph_objects.Figure, caption: Optional[str] = None) -> Tab


   .. py:method:: add_image(path: str, *, alt: str = '', caption: str | None = None, style: str = 'max-width:100%;height:auto;') -> Tab


   .. py:method:: add_raw_html(raw_html: str) -> Tab


   .. py:method:: add_svg(svg: str, *, caption: str | None = None) -> Tab


   .. py:method:: add_table(data: Mapping[str, Any], *, caption: Optional[str] = None, key_title: str = 'Key', value_title: str = 'Value', precision: int = 6) -> Tab

      
      Append a key/value HTML table made from a dict.
















      ..
          !! processed by numpydoc !!


   .. py:method:: add_text(text: str) -> Tab

