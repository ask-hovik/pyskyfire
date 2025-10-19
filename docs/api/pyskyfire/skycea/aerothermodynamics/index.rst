pyskyfire.skycea.aerothermodynamics
===================================

.. py:module:: pyskyfire.skycea.aerothermodynamics






Module Contents
---------------

.. py:data:: script_dir

.. py:data:: trans_path

.. py:class:: Aerothermodynamics(optimum, chemrep_map = None)

   
   Add all values in the optimum dict to self


   .. py:attribute:: optimum


   .. py:attribute:: chemrep_map


   .. py:method:: from_F_eps_Lstar(fu, ox, MR, p_c, F, eps, L_star, T_fu_in=298.15, T_ox_in=298.15, p_amb=101300.0, npts=15)
      :classmethod:


      Calculate optimal values using thrust, exit pressure and L-star



   .. py:method:: compute_aerothermodynamics(contour, Nt = 64)

      Create 2-D property maps on (x, T). Column 0 = equilibrium at that x.



   .. py:method:: get_T(x, T = None, h = None)


   .. py:method:: get_p(x, T = None, h = None)


   .. py:method:: get_M(x, T = None, h = None)


   .. py:method:: get_rho(x, T = None, h = None)


   .. py:method:: get_cp(x, T = None, h = None)


   .. py:method:: get_gamma(x, T = None, h = None)


   .. py:method:: get_h(x, T = None, h = None)


   .. py:method:: get_H(x, T = None, h = None)


   .. py:method:: get_a(x, T = None, h = None)


   .. py:method:: get_mu(x, T = None, h = None)


   .. py:method:: get_k(x, T = None, h = None)


   .. py:method:: get_Pr(x, T = None, h = None)


   .. py:method:: get_X(x)

      Interpolated mole-fraction dict at position x.



