pyskyfire.viz
=============

.. py:module:: pyskyfire.viz


Submodules
----------

.. toctree::
   :maxdepth: 1

   /autoapi/pyskyfire/viz/core/index
   /autoapi/pyskyfire/viz/embed_stl/index
   /autoapi/pyskyfire/viz/export_3d/index
   /autoapi/pyskyfire/viz/lookup_tables/index
   /autoapi/pyskyfire/viz/plot_common/index
   /autoapi/pyskyfire/viz/plot_regen/index
   /autoapi/pyskyfire/viz/plot_skycea/index
   /autoapi/pyskyfire/viz/report/index


Attributes
----------

.. autoapisummary::

   pyskyfire.viz.IndexLike
   pyskyfire.viz._PROP_INFO


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/viz/EmbedSTL
   /autoapi/pyskyfire/viz/PlotBase
   /autoapi/pyskyfire/viz/PlotBase
   /autoapi/pyskyfire/viz/PlotBase
   /autoapi/pyskyfire/viz/PlotBase
   /autoapi/pyskyfire/viz/PlotBase
   /autoapi/pyskyfire/viz/PlotBase
   /autoapi/pyskyfire/viz/PlotContour
   /autoapi/pyskyfire/viz/PlotCoolantArea
   /autoapi/pyskyfire/viz/PlotCoolantPressure
   /autoapi/pyskyfire/viz/PlotCoolantTemperature
   /autoapi/pyskyfire/viz/PlotEngineNetwork
   /autoapi/pyskyfire/viz/PlotHeatFlux
   /autoapi/pyskyfire/viz/PlotHydraulicDiameter
   /autoapi/pyskyfire/viz/PlotMoodyDiagram
   /autoapi/pyskyfire/viz/PlotPTDiagram
   /autoapi/pyskyfire/viz/PlotRadiusOfCurvature
   /autoapi/pyskyfire/viz/PlotResidualHistory
   /autoapi/pyskyfire/viz/PlotStationProperty
   /autoapi/pyskyfire/viz/PlotTemperatureProfile
   /autoapi/pyskyfire/viz/PlotThetaVsEpsilon
   /autoapi/pyskyfire/viz/PlotTransportProperty
   /autoapi/pyskyfire/viz/PlotVelocity
   /autoapi/pyskyfire/viz/PlotWallTemperature
   /autoapi/pyskyfire/viz/PlotdAdxCoolantArea
   /autoapi/pyskyfire/viz/PlotdAdxThermalCoolant
   /autoapi/pyskyfire/viz/PlotdAdxThermalHotGas
   /autoapi/pyskyfire/viz/Report
   /autoapi/pyskyfire/viz/Tab
   /autoapi/pyskyfire/viz/_Adder
   /autoapi/pyskyfire/viz/_Block
   /autoapi/pyskyfire/viz/_Node

.. autoapisummary::

   pyskyfire.viz.EmbedSTL
   pyskyfire.viz.PlotBase
   pyskyfire.viz.PlotBase
   pyskyfire.viz.PlotBase
   pyskyfire.viz.PlotBase
   pyskyfire.viz.PlotBase
   pyskyfire.viz.PlotBase
   pyskyfire.viz.PlotContour
   pyskyfire.viz.PlotCoolantArea
   pyskyfire.viz.PlotCoolantPressure
   pyskyfire.viz.PlotCoolantTemperature
   pyskyfire.viz.PlotEngineNetwork
   pyskyfire.viz.PlotHeatFlux
   pyskyfire.viz.PlotHydraulicDiameter
   pyskyfire.viz.PlotMoodyDiagram
   pyskyfire.viz.PlotPTDiagram
   pyskyfire.viz.PlotRadiusOfCurvature
   pyskyfire.viz.PlotResidualHistory
   pyskyfire.viz.PlotStationProperty
   pyskyfire.viz.PlotTemperatureProfile
   pyskyfire.viz.PlotThetaVsEpsilon
   pyskyfire.viz.PlotTransportProperty
   pyskyfire.viz.PlotVelocity
   pyskyfire.viz.PlotWallTemperature
   pyskyfire.viz.PlotdAdxCoolantArea
   pyskyfire.viz.PlotdAdxThermalCoolant
   pyskyfire.viz.PlotdAdxThermalHotGas
   pyskyfire.viz.Report
   pyskyfire.viz.Tab
   pyskyfire.viz._Adder
   pyskyfire.viz._Block
   pyskyfire.viz._Node


Functions
---------

.. autoapisummary::

   pyskyfire.viz._file_to_data_uri
   pyskyfire.viz._format_fluid
   pyskyfire.viz._format_iterable
   pyskyfire.viz._format_material
   pyskyfire.viz._format_number
   pyskyfire.viz._is_fluid
   pyskyfire.viz._is_material
   pyskyfire.viz._normalize_indices
   pyskyfire.viz._require_trimesh
   pyskyfire.viz.dict_to_table_html
   pyskyfire.viz.f_darcy
   pyskyfire.viz.format_value
   pyskyfire.viz.make_engine_gmsh
   pyskyfire.viz.render_engine_network


Package Contents
----------------

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

.. py:function:: _normalize_indices(thrust_chamber, circuit_index: IndexLike)

.. py:function:: _require_trimesh()

.. py:function:: dict_to_table_html(data: Mapping[str, Any], *, col_key: str = 'Key', col_val: str = 'Value', caption: str | None = None, precision: int = 6) -> str

   
   Return HTML for a simple 2-column table of a dict.
















   ..
       !! processed by numpydoc !!

.. py:function:: f_darcy(ReDh, Dh, x, roughness)

   
   Composite laminar–turbulent Darcy friction factor with smooth transition.


   :Parameters:

       **ReDh** : :class:`python:float`
           Hydraulic Reynolds number.

       **Dh** : :class:`python:float`
           Hydraulic diameter [m].

       **x** : :class:`python:float`
           Axial coordinate.

       **roughness** : :func:`python:callable` or :data:`python:None`
           Function returning surface roughness at ``x`` [m], or ``None`` for smooth wall.



   :Returns:

       :class:`python:float`
           Darcy friction factor (dimensionless).











   ..
       !! processed by numpydoc !!

.. py:function:: format_value(value: Any, precision: int = 6) -> str

   
   Human-friendly string for table cells.
















   ..
       !! processed by numpydoc !!

.. py:function:: make_engine_gmsh(thrust_chamber, filename='engine', display_channels=None)

.. py:function:: render_engine_network(engine_network, *, height: str = '900px', width: str = '100%', edge_length: int = 200, physics_settings: Dict | str | None = 'default', mass_flow_based_arrows: bool = False, station_mode: Literal['name', 'values', 'both', 'hidden'] = 'values') -> str

   
   Build a PyVis network and return an <iframe srcdoc="..."> HTML snippet
   suited for Report.add_raw_html(...). Produces a single-file embed.
















   ..
       !! processed by numpydoc !!

.. py:data:: IndexLike

.. py:data:: _PROP_INFO

