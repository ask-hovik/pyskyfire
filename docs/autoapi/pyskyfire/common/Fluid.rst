pyskyfire.common.Fluid
======================

.. py:class:: pyskyfire.common.Fluid(type, propellants, fractions, basis='mass', precision=3)

   
   Mixture helper for fluids.

   Stores a set of component names and composition.

   :Parameters:

       **type** : :class:`python:str`
           Role of the mixture (e.g., ``"fuel"``, ``"oxidizer"``, ``"coolant"``).

       **propellants** : :class:`python:list`\[:class:`python:str`]
           Component names (CoolProp canonical names, e.g., ``"Ethanol"``, ``"Water"``).

       **fractions** : :class:`python:list`\[:class:`python:float`]
           Composition fractions corresponding to ``propellants``. Interpreted
           as mass or mole fractions depending on ``basis``.

       **basis** : :class:`python:str`, :obj:`optional`
           ``"mass"`` or ``"mole"``. Default is ``"mass"``.

       **precision** : :class:`python:int`, :obj:`optional`
           Number of decimal digits when exporting mole fractions (e.g., in
           :meth:`coolprop_string`). Default is ``3``.

   :Attributes:

       **type** : :class:`python:str`
           Role marker for the mixture.

       **propellants** : :class:`python:list`\[:class:`python:str`]
           CoolProp component names.

       **fractions** : :class:`python:list`\[:class:`python:float`]
           Stored composition in the given ``basis``.

       **basis** : :class:`python:str`
           Basis of the stored composition (``"mass"`` or ``"mole"``).

       **precision** : :class:`python:int`
           Decimal places for exported mole fractions.









   .. seealso::

       
       :obj:`CoolProp`
           External property library used for molar masses and HEOS string format.
       
       
   .. rubric:: Notes

   Fractions are re-normalized internally when converting to mole basis to
   guard against small rounding errors.



   ..
       !! processed by numpydoc !!

   .. py:method:: as_mole_fractions()

      
      Convert the stored composition to mole fractions.





      :Returns:

          :class:`python:list`\[:class:`python:float`]
              Mole fractions for each component in :attr:`propellants`, rounded to
              :attr:`precision` decimals and re-normalized to sum to 1.0.




      :Raises:

          :obj:`ValueError`
              If :attr:`basis` is not ``"mass"`` or ``"mole"``.




      .. rubric:: Notes

      For a mass-basis input, the conversion uses
      :math:`x_i = \frac{w_i/M_i}{\sum_j w_j/M_j}` with molar masses
      :math:`M_i` from :meth:`molar_masses`.



      ..
          !! processed by numpydoc !!


   .. py:method:: coolprop_string()

      
      Format a CoolProp HEOS mixture string.





      :Returns:

          :class:`python:str`
              A string such as ``'HEOS::Ethanol[0.800]&Water[0.200]'``, where the
              bracketed values are **mole fractions** with :attr:`precision` digits.








      .. rubric:: Notes

      The composition is always exported on a **mole** basis, using
      :meth:`as_mole_fractions` for conversion when needed.



      ..
          !! processed by numpydoc !!


   .. py:method:: molar_masses()

      
      Return molar masses for each component.





      :Returns:

          :class:`python:list`\[:class:`python:float`]
              Molar masses in ``kg/mol`` for each name in :attr:`propellants`.








      .. rubric:: Notes

      Values are retrieved via ``CoolProp.PropsSI("M", name)``.



      ..
          !! processed by numpydoc !!

