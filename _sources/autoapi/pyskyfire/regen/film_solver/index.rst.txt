pyskyfire.regen.film_solver
===========================

.. py:module:: pyskyfire.regen.film_solver

.. autoapi-nested-parse::

   film_solver.py

   Contour-based Grisson-style film cooling model for pyskyfire.

   Key design choices in this version
   ----------------------------------
   - Geometry is read directly from thrust_chamber.contour.
   - Film-cooling inputs live on thrust_chamber.film_cooling.
   - ThrustChamber is expected to resolve film_cooling.x_fraction -> film_cooling.x
     before this solver is used.
   - Coolant liquid/vapor properties are pulled from CoolProp through the existing
     coolant transport object on the selected cooling circuit.
   - H2O and CO2 mole fractions used by the current Grisson gas-radiation model
     are taken directly from FilmCooling, not from combustion_transport.

   Important assumptions
   ---------------------
   - The injected coolant is a pure fluid.
   - Gas-side mean molecular weight is still expected to be available from
     thrust_chamber.combustion_transport.

   ..
       !! processed by numpydoc !!


Attributes
----------

.. autoapisummary::

   pyskyfire.regen.film_solver.PropsSI


Classes
-------

.. toctree::
   :hidden:

   /autoapi/pyskyfire/regen/film_solver/ContourGeometryEvaluator
   /autoapi/pyskyfire/regen/film_solver/GasEmittanceCalculator
   /autoapi/pyskyfire/regen/film_solver/GasProperties
   /autoapi/pyskyfire/regen/film_solver/GaseousFilmResults
   /autoapi/pyskyfire/regen/film_solver/GaseousFilmSolver
   /autoapi/pyskyfire/regen/film_solver/GrissonFilmCoolingModel
   /autoapi/pyskyfire/regen/film_solver/LiquidCoolantProperties
   /autoapi/pyskyfire/regen/film_solver/LiquidFilmResults
   /autoapi/pyskyfire/regen/film_solver/LiquidFilmSolver

.. autoapisummary::

   pyskyfire.regen.film_solver.ContourGeometryEvaluator
   pyskyfire.regen.film_solver.GasEmittanceCalculator
   pyskyfire.regen.film_solver.GasProperties
   pyskyfire.regen.film_solver.GaseousFilmResults
   pyskyfire.regen.film_solver.GaseousFilmSolver
   pyskyfire.regen.film_solver.GrissonFilmCoolingModel
   pyskyfire.regen.film_solver.LiquidCoolantProperties
   pyskyfire.regen.film_solver.LiquidFilmResults
   pyskyfire.regen.film_solver.LiquidFilmSolver


Functions
---------

.. autoapisummary::

   pyskyfire.regen.film_solver._build_coolant_properties
   pyskyfire.regen.film_solver._build_gas_properties
   pyskyfire.regen.film_solver._coolprop_fluid_name
   pyskyfire.regen.film_solver._coolprop_scalar
   pyskyfire.regen.film_solver._extract_gas_molecular_weight
   pyskyfire.regen.film_solver._require
   pyskyfire.regen.film_solver._try_call


Module Contents
---------------

.. py:function:: _build_coolant_properties(thrust_chamber, boundary_conditions, circuit_index: int) -> LiquidCoolantProperties

.. py:function:: _build_gas_properties(thrust_chamber, film_cooling) -> GasProperties

.. py:function:: _coolprop_fluid_name(coolant_transport) -> str

.. py:function:: _coolprop_scalar(output: str, name1: str, value1: float, name2: str, value2: float, fluid: str) -> float

.. py:function:: _extract_gas_molecular_weight(combustion_transport, x: float) -> float

.. py:function:: _require(value: Any, name: str) -> Any

.. py:function:: _try_call(obj: Any, method_name: str, *args, default=None)

.. py:data:: PropsSI
   :value: None


