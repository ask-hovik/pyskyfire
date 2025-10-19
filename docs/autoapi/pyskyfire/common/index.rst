pyskyfire.common
================

.. py:module:: pyskyfire.common


Submodules
----------

.. toctree::
   :maxdepth: 1

   /autoapi/pyskyfire/common/blocks/index
   /autoapi/pyskyfire/common/constants/index
   /autoapi/pyskyfire/common/engine_network/index
   /autoapi/pyskyfire/common/fluids/index
   /autoapi/pyskyfire/common/results/index
   /autoapi/pyskyfire/common/solids/index


Attributes
----------

.. autoapisummary::

   pyskyfire.common.ArrayLike
   pyskyfire.common.GRCop42
   pyskyfire.common.Inconel625
   pyskyfire.common.Inconel718
   pyskyfire.common.StainlessSteel304
   pyskyfire.common.T625
   pyskyfire.common.T625_C
   pyskyfire.common.T718_hi
   pyskyfire.common.TEOS
   pyskyfire.common.T_F
   pyskyfire.common.ZirconiumOxide
   pyskyfire.common.k42_mono
   pyskyfire.common.k625
   pyskyfire.common.k625_W
   pyskyfire.common.k718
   pyskyfire.common.k718_cryo
   pyskyfire.common.k718_hi
   pyskyfire.common.k718_table
   pyskyfire.common.k_304_piecewise
   pyskyfire.common.k_BTUin
   pyskyfire.common.log10poly_304_cryo
   pyskyfire.common.mills_304_highT


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/common/BoundaryConditions
   /autoapi/pyskyfire/common/ConstantModel
   /autoapi/pyskyfire/common/EngineNetwork
   /autoapi/pyskyfire/common/Fluid
   /autoapi/pyskyfire/common/FluidBlock
   /autoapi/pyskyfire/common/Log10PolynomialModel
   /autoapi/pyskyfire/common/MassFlowMergerBlock
   /autoapi/pyskyfire/common/MassFlowSplitterBlock
   /autoapi/pyskyfire/common/Material
   /autoapi/pyskyfire/common/PiecewiseModel
   /autoapi/pyskyfire/common/PolynomialModel
   /autoapi/pyskyfire/common/PropertyModel
   /autoapi/pyskyfire/common/PumpBlock
   /autoapi/pyskyfire/common/RegenBlock
   /autoapi/pyskyfire/common/Results
   /autoapi/pyskyfire/common/SignalBlock
   /autoapi/pyskyfire/common/SimpleDuctBlock
   /autoapi/pyskyfire/common/Station
   /autoapi/pyskyfire/common/Station
   /autoapi/pyskyfire/common/SumOfGaussiansModel
   /autoapi/pyskyfire/common/TabulatedModel
   /autoapi/pyskyfire/common/TransmissionBlock
   /autoapi/pyskyfire/common/TurbineBlock

.. autoapisummary::

   pyskyfire.common.BoundaryConditions
   pyskyfire.common.ConstantModel
   pyskyfire.common.EngineNetwork
   pyskyfire.common.Fluid
   pyskyfire.common.FluidBlock
   pyskyfire.common.Log10PolynomialModel
   pyskyfire.common.MassFlowMergerBlock
   pyskyfire.common.MassFlowSplitterBlock
   pyskyfire.common.Material
   pyskyfire.common.PiecewiseModel
   pyskyfire.common.PolynomialModel
   pyskyfire.common.PropertyModel
   pyskyfire.common.PumpBlock
   pyskyfire.common.RegenBlock
   pyskyfire.common.Results
   pyskyfire.common.SignalBlock
   pyskyfire.common.SimpleDuctBlock
   pyskyfire.common.Station
   pyskyfire.common.Station
   pyskyfire.common.SumOfGaussiansModel
   pyskyfire.common.TabulatedModel
   pyskyfire.common.TransmissionBlock
   pyskyfire.common.TurbineBlock


Functions
---------

.. autoapisummary::

   pyskyfire.common.steady_heating_analysis


Package Contents
----------------

.. py:function:: steady_heating_analysis(thrust_chamber, boundary_conditions, n_nodes=100, circuit_index=0, solver='newton', output=True)

   
   Run the steady heating analysis.


   :Parameters:

       **thrust_chamber** : :obj:`Any`
           Chamber model exposing the required geometry & property APIs.

       **boundary_conditions** : :obj:`BoundaryConditions`
           Coolant inlet boundary conditions.

       **n_nodes** : :class:`python:int`, :obj:`optional`
           Number of axial nodes. Default is 100.

       **circuit_index** : :class:`python:int`, :obj:`optional`
           Cooling-circuit index. Default is 0.

       **solver** : {'newton'}, :obj:`optional`
           Solver selector. Currently only ``'newton'`` is implemented.

       **output** : :ref:`bool <python:bltin-boolean-values>`, :obj:`optional`
           If True, print progress. Default is True.



   :Returns:

       :class:`python:dict`
           See :func:`solve_heat_exchanger_euler` for keys.











   ..
       !! processed by numpydoc !!

.. py:data:: ArrayLike

.. py:data:: GRCop42

.. py:data:: Inconel625

.. py:data:: Inconel718

.. py:data:: StainlessSteel304

.. py:data:: T625

.. py:data:: T625_C

.. py:data:: T718_hi

.. py:data:: TEOS

.. py:data:: T_F

.. py:data:: ZirconiumOxide

.. py:data:: k42_mono

.. py:data:: k625

.. py:data:: k625_W

.. py:data:: k718

.. py:data:: k718_cryo

.. py:data:: k718_hi

.. py:data:: k718_table

.. py:data:: k_304_piecewise

.. py:data:: k_BTUin

.. py:data:: log10poly_304_cryo

.. py:data:: mills_304_highT

