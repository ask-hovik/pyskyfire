pyskyfire.regen.thrust_chamber.FilmCooling
==========================================

.. py:class:: pyskyfire.regen.thrust_chamber.FilmCooling

   
   User-facing film cooling inputs attached to a ThrustChamber.


   :Parameters:

       **x_fraction** : :class:`python:float`
           Signed axial fraction in [-1, 1], where:
               -1 -> chamber start
                0 -> throat
               +1 -> nozzle exit
           This is resolved to ``x`` by ``ThrustChamber``.

       **coolant_mass_flow_rate** : :class:`python:float`
           Total coolant mass flow injected into the film [kg/s].

       **film_injection_perimeter** : :class:`python:float`
           Wetted injection perimeter [m].

       **liquid_absorptivity** : :class:`python:float`
           Liquid-film absorptivity for radiation [0..1].

       **mole_fraction_H2O** : :class:`python:float`
           User-supplied mole fraction of H2O used by the current Grisson
           gas-radiation submodel.

       **mole_fraction_CO2** : :class:`python:float`
           User-supplied mole fraction of CO2 used by the current Grisson
           gas-radiation submodel.

       **turbulence_intensity** : :class:`python:float`, :obj:`optional`
           Free-stream turbulence intensity used by Grisson's correlations.

       **x** : :class:`python:float` | :data:`python:None`, :obj:`optional`
           Absolute axial injection location [m]. This should be populated by
           ``ThrustChamber`` after it resolves ``x_fraction``.














   ..
       !! processed by numpydoc !!
