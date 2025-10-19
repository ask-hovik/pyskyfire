pyskyfire.common.PumpBlock
==========================

.. py:class:: pyskyfire.common.PumpBlock(name, st_in, st_out, overcome, p_base, input_p, load_fraction, n, eta, medium)

   Bases: :py:obj:`FluidBlock`

   .. autoapi-inheritance-diagram:: pyskyfire.common.PumpBlock
      :parts: 1
      :private-bases:


   
   Centrifugal/positive-displacement pump model (lumped).

   The block raises inlet pressure to meet a target set by the load
   (sum of downstream pressure drops) and computes shaft power and outlet state.

   :Parameters:

       **name** : :class:`python:str`
           Block name used to form signal keys (e.g., ``P_{name}``).

       **st_in** : :class:`python:str`
           Inlet station key.

       **st_out** : :class:`python:str`
           Outlet station key.

       **overcome** : :class:`python:list`\[:class:`python:str`]
           Names of blocks whose ``dp_*`` signals form the load to overcome.

       **p_base** : :class:`python:float`
           Baseline pressure added to the load [Pa].

       **input_p** : :class:`python:float`
           Upstream pressure offset already present [Pa].

       **load_fraction** : :class:`python:float`
           Fraction of total load to apply in this pump (0..1).

       **n** : :class:`python:float`
           Rotational speed (metadata).

       **eta** : :class:`python:float`
           Pump efficiency (0..1).

       **medium** : :class:`python:str`
           CoolProp fluid string.

   :Attributes:

       **dp_key** : :class:`python:str`
           Name of emitted power signal (``P_{name}``).

       **station_inputs, station_outputs, signal_inputs, signal_outputs** : :class:`python:list`\[:class:`python:str`]
           Network I/O metadata.









   .. seealso::

       
       :obj:`TurbineBlock`
           Consumes power to satisfy aggregate shaft load.
       
       



   ..
       !! processed by numpydoc !!

   .. py:method:: compute(stations, signals)

      
      Compute outlet state and required pump power.


      :Parameters:

          **stations** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              Network stations (must contain ``st_in``).

          **signals** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              Scalar signals (must contain all ``dp_{blk}`` listed in ``overcome``).



      :Returns:

          **stations_out** : :class:`python:dict`\[:class:`python:str`, :obj:`Station <pyskyfire.common.engine_network.Station>`]
              ``{st_out: Station}`` with updated pressure and temperature.

          **signals_out** : :class:`python:dict`\[:class:`python:str`, :class:`python:float`]
              ``{f"P_{name}": float}`` required pump power [W].




      :Raises:

          :obj:`KeyError`
              If required stations or signals are missing.




      .. rubric:: Notes

      Temperature property calls are offset by ``1e-3 K`` to avoid
      saturation-line singularities in CoolProp.



      ..
          !! processed by numpydoc !!

