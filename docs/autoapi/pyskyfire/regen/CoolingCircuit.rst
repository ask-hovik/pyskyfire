pyskyfire.regen.CoolingCircuit
==============================

.. py:class:: pyskyfire.regen.CoolingCircuit(name, contour, cross_section, span, placement, channel_height, coolant_transport, blockage_ratio=None)

   
   Representation of a cooling circuit following the chamber contour.

   Each circuit defines the geometry of one group of coolant channels,
   including their cross-sectional shape, span, and placement relative to
   the hot wall.

   :Parameters:

       **name** : :class:`python:str`
           Identifier of the circuit.

       **contour** : :obj:`Contour`
           Hot-gas wall geometry defining the outer boundary.

       **cross_section** : :obj:`ChannelSection`
           Cross-sectional geometry model providing A, Dh, and perimeters.

       **span** : :class:`python:tuple`\[:class:`python:float`, :class:`python:float`]
           Normalized start and end of circuit (-1 = chamber inlet, +1 = nozzle exit).

       **placement** : :obj:`ChannelPlacement`
           Strategy describing how channel centerlines are positioned.

       **channel_height** : :func:`python:callable`
           Function returning the local channel height [m].

       **coolant_transport** : :obj:`object`
           Object providing coolant thermophysical properties.

       **blockage_ratio** : :class:`python:float` or :term:`numpy:array_like`, :obj:`optional`
           Fraction of the sector blocked by solid wall or rib.

   :Attributes:

       **name** : :class:`python:str`
           Circuit name.

       **contour** : :obj:`Contour`
           Reference hot-gas contour.

       **cross_section** : :obj:`ChannelSection`
           Shape model used for thermal and hydraulic quantities.

       **placement** : :obj:`ChannelPlacement`
           Placement rule describing angular and radial positioning.

       **channel_height** : :func:`python:callable`
           Function returning height as a function of x.

       **coolant_transport** : :obj:`object`
           Provides coolant properties (k, μ, Cp, ρ, etc.).

       **span** : (:class:`python:float`, :class:`python:float`)
           Normalized start/end bounds.

       **direction** : :class:`python:int`
           +1 if flow is forward (increasing x), −1 if reversed.

       **A_coolant_vals, Dh_coolant_vals** : :obj:`ndarray <numpy.ndarray>`
           Precomputed local coolant area and hydraulic diameter.

       **dA_dx_thermal_exhaust_vals, dA_dx_thermal_coolant_vals** : :obj:`ndarray <numpy.ndarray>`
           Surface-area differentials on hot and cold sides.

       **radius_of_curvature_vals** : :obj:`ndarray <numpy.ndarray>`
           Local curvature radius for each axial node.









   .. seealso::

       
       :obj:`CoolingCircuitGroup`
           Groups multiple cooling circuits.
       :obj:`SectionProfiles`
           Bundles local geometry inputs for cross-section methods.
       
       
   .. rubric:: Notes

   The circuit precomputes local geometric properties for fast interpolation
   during steady-state or transient simulations.



   ..
       !! processed by numpydoc !!

   .. py:method:: A_coolant(x)


   .. py:method:: Dh_coolant(x)


   .. py:method:: _prof(centerline, local_coords)

      
      Assemble a `SectionProfiles` object for a given centerline.


      :Parameters:

          **centerline** : :obj:`ndarray <numpy.ndarray>`, :obj:`shape` (:obj:`N`, 3)
              Channel centerline coordinates (x, r, θ).

          **local_coords** : :obj:`ndarray <numpy.ndarray>`, :obj:`shape` (:obj:`N`, 3, 3)
              Local coordinate frames.



      :Returns:

          :obj:`SectionProfiles`
              Ready-to-use profile bundle for cross-section routines.











      ..
          !! processed by numpydoc !!


   .. py:method:: compute_geometry()

      
      Generate full 3D point-cloud representations for all channel centerlines.

      Each point cloud corresponds to one physical cooling channel.















      ..
          !! processed by numpydoc !!


   .. py:method:: compute_single_centerline()

      
      Generate OCC wire objects for each station along the first centerline.

      Intended for CAD or meshing visualization.















      ..
          !! processed by numpydoc !!


   .. py:method:: compute_volume()

      
      Compute total circuit volume by integrating local area along its centerline.





      :Returns:

          :class:`python:float`
              Total coolant volume [m³].











      ..
          !! processed by numpydoc !!


   .. py:method:: dA_dx_coolant(x)


   .. py:method:: dA_dx_thermal_coolant(x)


   .. py:method:: dA_dx_thermal_exhaust(x)


   .. py:method:: finalize()


   .. py:method:: precompute_thermal_properties()

      
      Precompute all cross-section-dependent thermal geometry arrays.

      Calculates effective surface-area derivatives, hydraulic diameter, and
      radius of curvature for interpolation during simulation.















      ..
          !! processed by numpydoc !!


   .. py:method:: radius_of_curvature(x)


   .. py:method:: set_blockage_ratio(blockage_ratio)

      
      blockage_ratio can be scalar or length-N array over x-domain.
















      ..
          !! processed by numpydoc !!


   .. py:method:: set_centerline(centerline_list)


   .. py:method:: set_centerline_test(centerline_list)


   .. py:method:: set_channel_height(heights)


   .. py:method:: set_channel_width(widths_rad)


   .. py:method:: set_t_wall_tot(t_wall_tot)


   .. py:method:: set_x_domain(x_domain)

