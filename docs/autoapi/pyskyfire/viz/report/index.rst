pyskyfire.viz.report
====================

.. py:module:: pyskyfire.viz.report


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/viz/report/Report
   /autoapi/pyskyfire/viz/report/Tab
   /autoapi/pyskyfire/viz/report/_Block
   /autoapi/pyskyfire/viz/report/_Block

.. autoapisummary::

   pyskyfire.viz.report.Report
   pyskyfire.viz.report.Tab
   pyskyfire.viz.report._Block
   pyskyfire.viz.report._Block


Functions
---------

.. autoapisummary::

   pyskyfire.viz.report._file_to_data_uri
   pyskyfire.viz.report._format_fluid
   pyskyfire.viz.report._format_iterable
   pyskyfire.viz.report._format_material
   pyskyfire.viz.report._format_number
   pyskyfire.viz.report._is_fluid
   pyskyfire.viz.report._is_material
   pyskyfire.viz.report.dict_to_table_html
   pyskyfire.viz.report.format_value


Module Contents
---------------

.. py:function:: _file_to_data_uri(path: str) -> str

.. py:function:: _format_fluid(fluid: Any, precision_pct: int = 2) -> str

.. py:function:: _format_iterable(it, precision: int = 6) -> str

.. py:function:: _format_material(mat: Any) -> str

   
   Render just the material's name, e.g., 'Inconel 718'.
   Keep it minimal per your request.
















   ..
       !! processed by numpydoc !!

.. py:function:: _format_number(x: float, precision: int = 6) -> str

.. py:function:: _is_fluid(obj: Any) -> bool

.. py:function:: _is_material(obj: Any) -> bool

   
   Duck-typing to detect your Material objects without importing from solids.py.
















   ..
       !! processed by numpydoc !!

.. py:function:: dict_to_table_html(data: Mapping[str, Any], *, col_key: str = 'Key', col_val: str = 'Value', caption: str | None = None, precision: int = 6) -> str

   
   Return HTML for a simple 2-column table of a dict.
















   ..
       !! processed by numpydoc !!

.. py:function:: format_value(value: Any, precision: int = 6) -> str

   
   Human-friendly string for table cells.
















   ..
       !! processed by numpydoc !!

