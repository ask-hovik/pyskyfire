pyskyfire.regen.thrust_chamber
==============================

.. py:module:: pyskyfire.regen.thrust_chamber






Module Contents
---------------

.. py:class:: Contour(xs, rs, name=None)

   
   Class for representing the inner contour of a rocket engine, from the beginning of the combustion chamber to the nozzle exit.

   :param xs: Array of x-positions, that the 'y' list corresponds to (m). Must be increasing values of x.
   :type xs: list
   :param rs: Array, containing local engine radius (m).
   :type rs: list

   .. attribute:: x_t

      x-position of the throat (m)

      :type: float

   .. attribute:: r_t

      Throat radius (m)

      :type: float

   .. attribute:: A_t

      Throat area (m2)

      :type: float

   .. attribute:: r_e

      Exit radius (m)

      :type: float

   .. attribute:: A_e

      Exit area (m2)

      :type: float

   .. attribute:: r_curvature_t

      Radius of curvature at the throat (m)

      :type: float


   .. py:attribute:: xs


   .. py:attribute:: rs


   .. py:attribute:: name
      :value: None



   .. py:property:: x_t


   .. py:property:: r_t


   .. py:property:: A_t


   .. py:property:: r_e


   .. py:property:: r_c


   .. py:property:: A_e


   .. py:property:: A_c


   .. py:property:: eps


   .. py:property:: eps_c

      Get the contraction ratio for the chamber
      :returns: contraction ratio
      :rtype: (float)


   .. py:method:: r(x)

      Get the distance from the centreline to the inner wall of the engine.
      :param x: x position (m)
      :type x: float

      :returns: Distance from engine centreline to edge of inner wall (m)
      :rtype: float



   .. py:method:: dr_dx(x)

      Get the slope of the engine wall, dr/dx.
      :param x: Axial position (m).
      :type x: float

      :returns: Rate of change of contour radius with respect to position, dr/dx
      :rtype: float



   .. py:method:: A(x)

      Get the flow area for the exhaust gas
      :param x: x position (m)
      :type x: float

      :returns: Flow area (m2)
      :rtype: float



   .. py:method:: normal_angle(x)

      Return the smaller angle [0..pi/2] between the outward normal to the
      contour and the plane perpendicular to the x-axis (i.e. the 'vertical'
      direction in this 2D cross-section).

      Geometrically:
        - Tangent: T = (1, dr/dx)
        - Outward normal: N = (-dr/dx, 1)
        - Plane perpendicular to x-axis ~ vertical direction, V = (0, 1)
        - cos(theta) = (N dot V) / (|N| * |V|) = 1 / sqrt((dr/dx)^2 + 1).



.. py:class:: ContourToroidalAerospike(xs_outer, rs_outer, xs_inner, rs_inner, *, name=None)

   Minimum-viable toroidal-aerospike contour class.

   * All **single-radius** properties (`r_t`, `r_c`, etc.) refer to the **outer** wall so
     existing bell-nozzle code keeps working.
   * All **areas** are **annular**:  :math:`A = \pi (r_{o}^{2} - r_{i}^{2})`.

   :param xs_outer: Axial sample points of the **outer** contour (m).
   :type xs_outer: array-like
   :param rs_outer: Outer radius at the given ``xs_outer`` (m).
   :type rs_outer: array-like
   :param xs_inner: Axial sample points of the **inner** contour (m).
   :type xs_inner: array-like
   :param rs_inner: Inner radius at the given ``xs_inner`` (m).
   :type rs_inner: array-like
   :param name: Human-readable identifier.
   :type name: str, optional


   .. py:attribute:: xs_outer


   .. py:attribute:: rs_outer


   .. py:attribute:: xs_inner


   .. py:attribute:: rs_inner


   .. py:attribute:: name
      :value: None



   .. py:property:: x_t

      Axial location of the **outer** throat (m).


   .. py:property:: r_t

      Throat radius of the **outer** wall (m).


   .. py:property:: r_e

      Outer radius at exit (m).


   .. py:property:: r_c

      Outer radius at chamber start (m).


   .. py:method:: A(x)

      Annular flow area at *x* (m²).



   .. py:property:: A_t

      Annular area at the throat (m²).


   .. py:property:: A_e


   .. py:property:: A_c


   .. py:property:: eps

      Area ratio exit / throat.


   .. py:property:: eps_c

      Chamber contraction ratio.


   .. py:method:: dr_dx(x, which='outer')

      Radial slope ``dr/dx`` at *x* on chosen wall.



   .. py:method:: normal_angle(x, which='outer')

      Angle between outward normal and vertical plane (rad).



   .. py:method:: r(x)

      Return **outer** radius at *x* (m). Provided for API compatibility.



.. py:class:: Wall(material, thickness, name=None)

   
   Object for representing an engine wall.

   :param name: Name of layer, for example "Chrome Coating" or "Main Wall"
   :type name: str
   :param material: Material object to define the material the wall is made of.
   :type material: Material
   :param thickness: Thickness of the wall (m). Can be a constant float, or a function of position, i.e. t(x).
   :type thickness: float or callable


   .. py:attribute:: name
      :value: None



   .. py:attribute:: material


   .. py:method:: thickness(x)

      Get the thickness of the wall at a position x.

      :param x: Axial position along the engine (m)
      :type x: float

      :returns: Wall thickness (m)
      :rtype: float



.. py:class:: WallGroup(walls=None)

   
   Class that manages multiple Wall objects.

   :param walls: A list of Wall objects. If None, defaults to an empty list.
   :type walls: list


   .. py:attribute:: walls
      :value: None



   .. py:method:: total_thickness(x)

      Returns the sum of the thicknesses of all walls at the given position x.

      :param x: Axial position (m)
      :type x: float

      :returns: The total thickness of all walls (m)
      :rtype: float



.. py:class:: CoolingCircuit(name, contour, cross_section, span, placement, channel_height, coolant_transport, blockage_ratio=None)

   .. py:attribute:: name


   .. py:attribute:: contour


   .. py:attribute:: cross_section


   .. py:attribute:: placement


   .. py:attribute:: channel_height


   .. py:attribute:: coolant_transport


   .. py:attribute:: blockage_ratio
      :value: None



   .. py:method:: precompute_thermal_properties()


   .. py:method:: compute_volume()


   .. py:method:: compute_single_centerline()


   .. py:method:: compute_geometry()


   .. py:method:: dA_dx_thermal_exhaust(x)


   .. py:method:: dA_dx_thermal_coolant(x)


   .. py:method:: A_coolant(x)


   .. py:method:: dA_dx_coolant(x)


   .. py:method:: Dh_coolant(x)


   .. py:method:: radius_of_curvature(x)


   .. py:method:: set_centerline_test(centerline_list)


   .. py:method:: set_centerline(centerline_list)


   .. py:method:: set_channel_width(widths_rad)


   .. py:method:: set_channel_height(heights)


   .. py:method:: set_t_wall_tot(t_wall_tot)


   .. py:method:: set_blockage_ratio(blockage_ratio)

      blockage_ratio can be scalar or length-N array over x-domain.



   .. py:method:: set_x_domain(x_domain)


   .. py:method:: finalize()


.. py:function:: old_radius_of_curvature(points)

   Calculate the radius of curvature along a line given by an array of points.
   Points: numpy array of shape (N, 3) representing [x, r, theta]
   :returns: numpy array of radius of curvature at each point (length N)
   :rtype: radii


.. py:function:: radius_of_curvature(points, axis = 'x', eps = 1e-12)

   Signed radius of curvature for a curve expressed in cylindrical coordinates
   [x, r, θ].

   • Positive  → curve bends *away* from the symmetry axis
   • Negative  → curve bends *toward* the symmetry axis
   • np.inf    → locally straight (|κ| below `eps`)

   :param points: [[x, r, theta], …] ordered along the curve.
   :type points: (N, 3) ndarray
   :param axis: Which coordinate is the symmetry axis.  Default 'x'.
   :type axis: {'x', 'y', 'z'}, optional
   :param eps: Curvature values with |κ| < eps are treated as zero (straight).
   :type eps: float, optional

   :returns: **R** -- Signed radius of curvature at each sample.
   :rtype: (N,) ndarray


.. py:class:: CoolingCircuitGroup(circuit_list, configuration=None)

   .. py:attribute:: circuits


   .. py:method:: number_of_channels(x, *, occluding_only=False)

      Return the total number of cooling channels active at a given x position.

      :param x: The axial x position along the engine.
      :type x: float

      :returns: Total number of cooling channels active at the given x position.
      :rtype: int



.. py:class:: ChannelPlacement(n_channel_positions, channel_width=None, occludes = True)

   Bases: :py:obj:`abc.ABC`


   Helper class that provides a standard way to create an ABC using
   inheritance.


   .. py:attribute:: n_channel_positions


   .. py:attribute:: channel_width
      :value: None



   .. py:attribute:: occludes
      :value: True



   .. py:method:: compute_centerline_radius(x, contour, wall_group)
      :abstractmethod:


      Given axial coordinate x, the hot-gas contour and the wall stack,
      return the r-coordinate of the coolant channel centerline.



   .. py:method:: channel_count()


.. py:class:: SurfacePlacement(n_channel_positions)

   Bases: :py:obj:`ChannelPlacement`


   Helper class that provides a standard way to create an ABC using
   inheritance.


   .. py:attribute:: n_channels_per_leaf
      :value: 1



   .. py:method:: compute_centerline_radius(x, contour, wall_group)

      Given axial coordinate x, the hot-gas contour and the wall stack,
      return the r-coordinate of the coolant channel centerline.



   .. py:attribute:: n_channel_positions


   .. py:attribute:: channel_width
      :value: None



   .. py:attribute:: occludes
      :value: True



   .. py:method:: channel_count()


.. py:class:: InternalPlacement(n_channel_positions, n_channels_per_leaf, *, channel_width, occludes = False)

   Bases: :py:obj:`ChannelPlacement`


   In-chamber heat-exchanger channels.  May have their own width law.


   .. py:attribute:: n_channels_per_leaf


   .. py:method:: compute_centerline_radius(x, contour, wall_group)

      Given axial coordinate x, the hot-gas contour and the wall stack,
      return the r-coordinate of the coolant channel centerline.



   .. py:attribute:: n_channel_positions


   .. py:attribute:: channel_width
      :value: None



   .. py:attribute:: occludes
      :value: True



   .. py:method:: channel_count()


.. py:class:: ThrustChamber(contour, wall_group, cooling_circuit_group, combustion_transport, optimal_values=None, roughness=1.5e-05, K_factor=0.3, n_nodes=50, h_gas_corr=1.0, h_cold_corr=1.0)

   
   :param contour: The hot-gas contour of the engine
   :type contour: Contour
   :param walls: Collection of walls, must have walls.total_thickness(x)
   :type walls: WallCollection
   :param cooling_circuits: Master container that holds individual CoolingCircuit objects
   :type cooling_circuits: CircuitMaster
   :param channel_height: A function returning the channel height h(x)
   :type channel_height: callable
   :param n_nodes: Number of axial subdivisions to use for constructing centerlines
   :type n_nodes: int


   .. py:attribute:: contour


   .. py:attribute:: wall_group


   .. py:attribute:: cooling_circuit_group


   .. py:attribute:: combustion_transport


   .. py:attribute:: n_nodes
      :value: 50



   .. py:attribute:: optimal_values
      :value: None



   .. py:attribute:: h_gas_corr
      :value: 1.0



   .. py:attribute:: h_cold_corr
      :value: 1.0



   .. py:attribute:: K_factor
      :value: 0.3



   .. py:method:: build_circuit_x_domain()

      Build the x-domain for each cooling circuit by converting its fractional span
      into actual x-values. The sign and ordering of the span determine the coolant flow direction.
      This function uses the overall engine x-range from the contour.



   .. py:method:: build_channel_centerlines(mode='sim')

      Build centerline splines for each CoolingCircuit.
      For each circuit, use its pre-built x_domain.
      Each circuit is assigned angles in an interleaved fashion.



   .. py:method:: build_channel_widths()

      Compute the channel widths (in radians) for each cooling circuit.
      Uses each circuit's pre-built x_domain and the new number_of_channels(x)
      function to determine the total active channels at each x position.



   .. py:method:: build_channel_heights()

      Compute the channel heights for each cooling circuit along its pre-built x_domain.
      Evaluate the channel height function at each x in the circuit's domain.



   .. py:method:: build_t_wall_tot()

      Build an array of total wall thicknesses along each circuit's x-domain and
      assign it to the corresponding cooling circuit using set_t_wall_tot.



   .. py:method:: roughness(x)

      Get the channel roughness, at a position, x.



.. py:function:: interleaved_indices(circuit_counts)

   Given a list of circuit_counts = [n0, n1, ..., nK], produce an array
   'owners' of length sum(circuit_counts), where each index i is assigned
   to exactly one circuit in an interleaved ratio of n0 : n1 : ... : nK.

   Example: circuit_counts = [30, 60].
   Then we have total=90, ratio=1:2.  The owners array might look like
     [0,1,1, 0,1,1, 0,1,1, ...]
   So circuit #0 gets 30 slots, circuit #1 gets 60 slots, interleaved 1:2.


