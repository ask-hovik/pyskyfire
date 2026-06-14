pyskyfire.pump.impeller_new.ImpellerInputs
==========================================

.. py:class:: pyskyfire.pump.impeller_new.ImpellerInputs

   
   User-facing design inputs for a first-pass centrifugal impeller.


   :Parameters:

       **Q:**
           Design-point volumetric flow rate [m^3/s].

       **H:**
           Design-point pump head [m].

       **n:**
           Rotational speed [rpm].

       **d1:**
           Optional impeller inlet / suction-eye outer diameter [m].  If omitted,
           it is estimated from the selected inlet flow coefficient and hub ratio.

       **dn:**
           Optional hub diameter at impeller inlet [m].  If omitted, ``hub_ratio``
           times ``d1`` is used.

       **blade_count:**
           Optional number of main blades.  If omitted, a conservative heuristic is
           used.  For real hardware this should normally be selected deliberately.

       **beta1_deg, beta2_deg:**
           Inlet and outlet blade metal angles [deg], measured relative to the
           circumferential direction in the usual pump convention.  These are used
           only to generate a camber surface here.

       **phi1:**
           Inlet meridional-flow coefficient used when estimating ``d1``:
           ``phi1 = c_m1/u_1``.

       **fd1:**
           Safety/enlargement factor applied to the calculated minimum inlet
           diameter to account for blockage and non-uniformity.

       **hub_ratio:**
           ``dn/d1`` used if ``dn`` is not provided.

       **n_streamlines, n_points:**
           Resolution of the generated meridional/blade camber surface.














   ..
       !! processed by numpydoc !!

   .. py:method:: validate() -> None

