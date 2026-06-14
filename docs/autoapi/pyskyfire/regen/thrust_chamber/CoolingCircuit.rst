pyskyfire.regen.thrust_chamber.CoolingCircuit
=============================================

.. py:class:: pyskyfire.regen.thrust_chamber.CoolingCircuit(name, contour, cross_section, span: list[float], placement, channel_height, walls, coolant_transport, roughness)

   
   Simulation-only representation of a cooling circuit.













   .. rubric:: Notes

   - No OCC / geometry generation here.
   - Does NOT compute true local frames; those are a visualization concern.



   ..
       !! processed by numpydoc !!

   .. py:method:: A_coolant(x)


   .. py:method:: Dh_coolant(x)


   .. py:method:: R_coolant_per_len(x: float, h_c: float, T_wall_rep: float) -> float

      
      Effective coolant-side thermal resistance per unit *axial length* [K m / W].

      Cross-section returns resistance per unit channel length (s). Convert to per-x using:
          R_x = R_s / (ds/dx)















      ..
          !! processed by numpydoc !!


   .. py:method:: _prof(centerline: numpy.ndarray) -> pyskyfire.regen.cross_section.SectionProfiles

      
      Build a SectionProfiles with trivial frames (analytics don't use them).
















      ..
          !! processed by numpydoc !!


   .. py:method:: compute_volume()

      
      Total circuit volume (all physical lanes) via ∫ A ds, scaled by lane count.
















      ..
          !! processed by numpydoc !!


   .. py:method:: dA_dx_coolant(x)


   .. py:method:: dA_dx_thermal_coolant(x)


   .. py:method:: dA_dx_thermal_exhaust(x)


   .. py:method:: ds_dx(x)


   .. py:method:: finalize()


   .. py:method:: local_sector_angle(x)


   .. py:method:: precompute_thermal_properties()

      
      Precompute A, Dh, and thermal perimeters along a representative centerline.
















      ..
          !! processed by numpydoc !!


   .. py:method:: radius_of_curvature(x)


   .. py:method:: roughness(x)


   .. py:method:: section_profiles_at(x: float) -> pyskyfire.regen.cross_section.SectionProfiles

      
      Single-station SectionProfiles for use in local closures (rib/fin, LUTs, etc.).
      Returns arrays of length 1.
















      ..
          !! processed by numpydoc !!


   .. py:method:: set_centerline(centerline_list)

      
      Provide one or more centerlines as arrays of shape (N, 3): [x, r, theta].
      For simulation we use the first as the representative path.
















      ..
          !! processed by numpydoc !!


   .. py:method:: set_channel_height(heights)


   .. py:method:: set_channel_local_sector(local_sectors)


   .. py:method:: set_channel_width(widths_rad)


   .. py:method:: set_x_domain(x_domain)


   .. py:method:: total_thickness(x)


   .. py:method:: wedge_angle(x)

