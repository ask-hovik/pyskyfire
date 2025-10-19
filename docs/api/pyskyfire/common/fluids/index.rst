pyskyfire.common.fluids
=======================

.. py:module:: pyskyfire.common.fluids




Module Contents
---------------

.. py:class:: Fluid(type, propellants, fractions, basis='mass', precision=3)

   
   :param type: "fuel" | "oxidizer" | "coolant" etc.
   :type type: str
   :param propellants: Component names (CoolProp canonical names)
   :type propellants: list[str]
   :param fractions: Fractions (mass or mole, depending on basis)
   :type fractions: list[float]
   :param basis: "mass" or "mole"
   :type basis: str
   :param precision: number of decimal digits when exporting mole fractions
   :type precision: int


   .. py:attribute:: type


   .. py:attribute:: propellants


   .. py:attribute:: fractions


   .. py:attribute:: basis
      :value: 'mass'



   .. py:attribute:: precision
      :value: 3



   .. py:method:: molar_masses()

      Return molar masses [kg/mol] for each propellant.



   .. py:method:: as_mole_fractions()

      Convert stored fractions to mole fractions.



   .. py:method:: coolprop_string()

      Return a CoolProp HEOS string like 'HEOS::Ethanol[0.800]&Water[0.200]'.



