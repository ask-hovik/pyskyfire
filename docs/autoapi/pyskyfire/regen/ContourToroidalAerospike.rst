pyskyfire.regen.ContourToroidalAerospike
========================================

.. py:class:: pyskyfire.regen.ContourToroidalAerospike(xs_outer, rs_outer, xs_inner, rs_inner, *, name=None)

   
   Axisymmetric contour describing a toroidal-aerospike geometry.

   Defines both the inner and outer wall surfaces, allowing evaluation of
   annular areas and local geometric derivatives.

   :Parameters:

       **xs_outer, rs_outer** : :term:`numpy:array_like`
           Axial and radial coordinates of the outer wall [m].

       **xs_inner, rs_inner** : :term:`numpy:array_like`
           Axial and radial coordinates of the inner wall [m].

       **name** : :class:`python:str`, :obj:`optional`
           Identifier for this contour.

   :Attributes:

       **xs_outer, rs_outer** : :obj:`ndarray <numpy.ndarray>`
           Outer wall geometry arrays [m].

       **xs_inner, rs_inner** : :obj:`ndarray <numpy.ndarray>`
           Inner wall geometry arrays [m].

       **_dr_dx_outer, _dr_dx_inner** : :obj:`ndarray <numpy.ndarray>`
           Local wall slope arrays for each surface.

       **name** : :class:`python:str` or :data:`python:None`
           Descriptive label.

       **A_t, A_c, A_e** : :class:`python:float`
           Annular areas at throat, chamber, and exit [m²].

       **eps, eps_c** : :class:`python:float`
           Expansion and contraction ratios (dimensionless).









   .. seealso::

       
       :obj:`Contour`
           Single-wall version used in standard bell-nozzle engines.
       
       
   .. rubric:: Notes

   The inner radius is validated to always remain below the outer radius.



   ..
       !! processed by numpydoc !!

   .. py:method:: A(x)

      
      Return the local annular flow area.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].



      :Returns:

          :class:`python:float`
              Annular cross-sectional area [m²].











      ..
          !! processed by numpydoc !!


   .. py:method:: __setattr__(key, value)


   .. py:method:: _interp_inner(x)

      
      Linear interpolation of inner radius at *x*.
















      ..
          !! processed by numpydoc !!


   .. py:method:: _interp_outer(x)

      
      Linear interpolation of outer radius at *x*.
















      ..
          !! processed by numpydoc !!


   .. py:method:: dr_dx(x, which='outer')

      
      Return the wall slope :math:`dr/dx` at `x` for either wall.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].

          **which** : {'outer', 'inner'}, :obj:`optional`
              Which wall to evaluate. Default is 'outer'.



      :Returns:

          :class:`python:float`
              Local slope of the selected wall.











      ..
          !! processed by numpydoc !!


   .. py:method:: normal_angle(x, which='outer')

      
      Return the normal-plane angle for either wall.


      :Parameters:

          **x** : :class:`python:float`
              Axial coordinate [m].

          **which** : {'outer', 'inner'}, :obj:`optional`
              Wall selector. Default 'outer'.



      :Returns:

          :class:`python:float`
              Angle between outward normal and vertical plane [rad].











      ..
          !! processed by numpydoc !!


   .. py:method:: r(x)

      
      Return **outer** radius at *x* (m). Provided for API compatibility.
















      ..
          !! processed by numpydoc !!


   .. py:property:: A_c


   .. py:property:: A_e


   .. py:property:: A_t

      
      Annular area at the throat (m²).
















      ..
          !! processed by numpydoc !!


   .. py:property:: eps

      
      Area ratio exit / throat.
















      ..
          !! processed by numpydoc !!


   .. py:property:: eps_c

      
      Chamber contraction ratio.
















      ..
          !! processed by numpydoc !!


   .. py:property:: r_c

      
      Outer radius at chamber start (m).
















      ..
          !! processed by numpydoc !!


   .. py:property:: r_e

      
      Outer radius at exit (m).
















      ..
          !! processed by numpydoc !!


   .. py:property:: r_t

      
      Throat radius of the **outer** wall (m).
















      ..
          !! processed by numpydoc !!


   .. py:property:: x_t

      
      Axial location of the **outer** throat (m).
















      ..
          !! processed by numpydoc !!

