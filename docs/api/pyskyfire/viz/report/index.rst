pyskyfire.viz.report
====================

.. py:module:: pyskyfire.viz.report






Module Contents
---------------

.. py:function:: format_value(value, precision = 6)

   Human-friendly string for table cells.


.. py:function:: dict_to_table_html(data, *, col_key = 'Key', col_val = 'Value', caption = None, precision = 6)

   Return HTML for a simple 2-column table of a dict.


.. py:class:: Tab

   .. py:attribute:: title
      :type:  str


   .. py:attribute:: blocks
      :type:  List[_Block]
      :value: []



   .. py:method:: add_text(text)


   .. py:method:: add_figure(fig, caption = None)


   .. py:method:: add_table(data, *, caption = None, key_title = 'Key', value_title = 'Value', precision = 6)

      Append a key/value HTML table made from a dict.



   .. py:method:: add_raw_html(raw_html)


   .. py:method:: add_image(path, *, alt = '', caption = None, style = 'max-width:100%;height:auto;')


   .. py:method:: add_svg(svg, *, caption = None)


.. py:class:: Report(title = 'pyskyfire Report')

   .. py:attribute:: title
      :value: 'pyskyfire Report'



   .. py:attribute:: tabs
      :type:  List[Tab]
      :value: []



   .. py:method:: add_tab(title)


   .. py:method:: save_html(path)


