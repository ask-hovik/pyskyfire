pyskyfire.pump
==============

.. py:module:: pyskyfire.pump


Submodules
----------

.. toctree::
   :maxdepth: 1

   /api/pyskyfire/pump/constants/index
   /api/pyskyfire/pump/impeller/index
   /api/pyskyfire/pump/plot/index
   /api/pyskyfire/pump/pump_cli/index
   /api/pyskyfire/pump/utils/index




Package Contents
----------------

.. py:class:: Impeller(Q, H, n, material=None)

   Bases: :py:obj:`object`


   Execute the project of a centrifugal pump.

   Take input variables and execute the project.

   :param Q (float): flow rate [m^3/s]
   :param H (float) head [m]
   :param n (int): rotational speed [1/min]


   .. py:attribute:: Q


   .. py:attribute:: H


   .. py:attribute:: n


   .. py:attribute:: z
      :value: 6



   .. py:attribute:: geometry


   .. py:method:: main_dimensions()

      Calculate the main dimensions of the impeller



   .. py:method:: compute_geometry()

      Compute meridional & and blade geometry



   .. py:method:: compute_mass()


   .. py:method:: plot_3d(a=1, b=0, c=0, beta_1B=45, beta_2B=50, num_blades=6)

      Plot a 3D view of the impeller geometry.



