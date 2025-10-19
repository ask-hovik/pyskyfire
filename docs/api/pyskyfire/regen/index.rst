pyskyfire.regen
===============

.. py:module:: pyskyfire.regen


Submodules
----------

.. toctree::
   :maxdepth: 1

   /api/pyskyfire/regen/channel_height/index
   /api/pyskyfire/regen/contour/index
   /api/pyskyfire/regen/cross_section/index
   /api/pyskyfire/regen/physics/index
   /api/pyskyfire/regen/solver/index
   /api/pyskyfire/regen/thrust_chamber/index








Package Contents
----------------

.. py:function:: h_gas_bartz_enthalpy_driven(k_gr, D_hyd, Cp_gr, mu_gr, mdot_g, A_chmb, T_g, T_gr)

   Bartz correlation for hot gas side heat transfer coefficient

   :param k_gr is the:

   Subscripts:
       g denotes free stream gas properties
       r denotes reference enthalpy conditions which are averaged between the free stream and wall metal conditions


.. py:function:: h_gas_bartz(D_t, mu_g, cp_g, Pr_g, p_c, c_star, A_t, A_x, sigma)

   Bartz correlation from wikipedia without the curvature correction


.. py:function:: sigma(T_hw, T_c, gamma_g, M_g, omega)

   Dimensionless parameter accounting for variation of gas properties across the boundary layer


.. py:function:: h_coolant_colburn(k_cf, D_c, Cp_cr, mu_cf, mdot_c, A_c, phi_curv=1)

   Colburn correlation for coolant side heat transfer coefficient

   :param k_Cf is the conductivity in the coolant film:
   :param D_C is the coolant tube hydraulic diameter:
   :param cp_cr is then maybe the averaged cp between the film condition and the wall metal condition:
   :param mu_Cf is the viscosity at the film conditions:
   :param mdot_c is the coolant mass flow rate:
   :param A_c is the cross sectional area of coolant flow:

   Subscripts:
       Cf denotes film coolant film condition
       c denotes bulk coolant conditions
       r denotes reference enthalpy condition which are averaged between the free stream and wall metal conditions


.. py:function:: u_coolant(rho, mdot_c_single_channel, A_cool)

   coolant velocity as a function of x


.. py:function:: reynolds(rho, u, L, mu)

.. py:function:: phi_curv(Re_c, D_c, R_curv)

   Curvature effect on the cold side heat transfer coefficient


.. py:data:: ReDh_laminar
   :value: 2300


.. py:data:: ReDh_turbulent
   :value: 3500


.. py:function:: f_darcy_laminar(ReDh, Dh, x)

.. py:function:: f_darcy_turbulent(ReDh, Dh, x, roughness)

.. py:function:: f_darcy(ReDh, Dh, x, roughness)

.. py:function:: T_aw(gamma, M_inf, T_inf, Pr)

   Adiabatic wall temperature, sometimes called recovery temperature T_r in some sources
   TODO: Implement a more robust version that could also handle laminar flow, in case this function is suddenly
   used in a function where this makes sense.


.. py:class:: BoundaryConditions(T_coolant_in, p_coolant_in, mdot_coolant)

   Object for storing boundary conditions for the solver.

   :param T_coolant_in: Static(?) temperature of coolant at cooling channel inlet (K)
   :type T_coolant_in: float
   :param p_coolant_in: Static(?) pressure of coolant at cooling channel inlet (Pa)
   :type p_coolant_in: float
   :param mdot_coolant: mass flow rate of coolant through cooling channel inlet (kg/s)
   :type mdot_coolant: float


   .. py:attribute:: T_coolant_in


   .. py:attribute:: p_coolant_in


   .. py:attribute:: mdot_coolant


.. py:class:: HeatExchangerPhysics(thrust_chamber, circuit_index)

   Encapsulates the physical calculations for the heat exchanger.
   This class is responsible for computing the heat fluxes and rates based
   on the given engine properties and operating conditions.


   .. py:attribute:: thrust_chamber


   .. py:attribute:: circuit_index


   .. py:attribute:: counter
      :value: 0



   .. py:method:: dQ_hot_dx(x, T_hw)

      Computes the heat transfer rate from the hot gas to the wall. Bartz equation, for example



   .. py:method:: dQ_cond_dx(x, T_hw, T_cw)


   .. py:method:: dQ_cold_dx(x, T_cw, T_cool)

      Computes the heat transfer rate from the wall to the coolant.



   .. py:method:: coolant_temperature_rate(T_cool, p_cool, dQ_cold_dx)

      Computes the rate of change of the coolant temperature.



   .. py:method:: coolant_pressure_rate(x, T_cool, p_cool)

      Computes the rate of change of the coolant pressure due to frictional losses.



   .. py:method:: interface_temperatures(x, T_hw, T_cw)

      Returns a list [T0, T1, ..., Tn] of wall-stack temperatures,
      where T0=T_hw, Tn=T_cw, and in between are each wall interface.



.. py:function:: solve_heat_exchanger_euler(thrust_chamber, boundary_conditions, n_nodes, circuit_index, output, log_residuals=True)

   Solve the 1D steady-state heat exchanger from x=0 to x=x_domain[-1].

   :param engine: your engine object containing geometry & property methods
   :param n_nodes: number of axial nodes along the thrust chamber

   :returns:    "x"          -> 1D array of axial positions
                "T"          -> 2D array of temperatures, shape (n_nodes, 3):
                                columns = [T_coolant, T_wall_cold_side, T_wall_hot_side]
                "T_coolant"  -> 1D array of coolant temperatures
                "p_coolant"  -> 1D array of coolant pressures
                "dQ_dA"      -> 1D array of local heat fluxes (W/m^2)
                "velocity"   -> 1D array of coolant velocities (m/s)
   :rtype: A dictionary with keys


.. py:function:: analyse_residuals(residual_log, n_cells, p=2)

   :param residual_log: The list returned by `solve_channel`.  If None, nothing happens.
   :type residual_log: list | None
   :param n_cells: Number of axial nodes in the simulation.
   :type n_cells: int
   :param p: Order of the global norm: 2 for RMS, np.inf for L∞, etc.
   :type p: int | float

   :returns: * **history** (*(n_iter,) ndarray | None*) -- Global residual norm per iteration 0..k_max.
             * **final_per_cell** (*(n_cells,) ndarray | None*) -- Residual magnitude in each cell at the last local iteration.


.. py:function:: steady_heating_analysis(thrust_chamber, boundary_conditions, n_nodes=100, circuit_index=0, solver='newton', output=True)

   Run the steady heating analysis.

   :param engine: Engine object with geometry, boundary conditions, and physics.
   :param n_nodes: Number of nodes (for Newton) or resolution for post-processing (for Radau).
   :param solver: String, currently only "newton" available

   :returns: A dictionary with simulation results.


.. py:class:: SectionProfiles

   All geometry profiles along the centerline (length N).


   .. py:attribute:: h
      :type:  numpy.ndarray


   .. py:attribute:: theta
      :type:  numpy.ndarray


   .. py:attribute:: t_wall
      :type:  numpy.ndarray


   .. py:attribute:: centerline
      :type:  numpy.ndarray


   .. py:attribute:: local_coords
      :type:  numpy.ndarray


   .. py:attribute:: blockage_ratio
      :type:  numpy.ndarray


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


.. py:class:: SectionProfiles

   All geometry profiles along the centerline (length N).


   .. py:attribute:: h
      :type:  numpy.ndarray


   .. py:attribute:: theta
      :type:  numpy.ndarray


   .. py:attribute:: t_wall
      :type:  numpy.ndarray


   .. py:attribute:: centerline
      :type:  numpy.ndarray


   .. py:attribute:: local_coords
      :type:  numpy.ndarray


   .. py:attribute:: blockage_ratio
      :type:  numpy.ndarray


.. py:class:: ChannelSection(n_points = 16)

   Bases: :py:obj:`abc.ABC`


   Common interface all cross-section shapes must implement.


   .. py:attribute:: n_points
      :value: 16



   .. py:method:: A_coolant(prof)
      :abstractmethod:



   .. py:method:: Dh_coolant(prof)
      :abstractmethod:



   .. py:method:: P_thermal(prof)
      :abstractmethod:



   .. py:method:: P_coolant(prof)
      :abstractmethod:



   .. py:method:: compute_cross_section(prof, i)
      :abstractmethod:



.. py:class:: CrossSectionSquared(n_points = 8)

   Bases: :py:obj:`ChannelSection`


   Common interface all cross-section shapes must implement.


   .. py:method:: A_coolant(prof)


   .. py:method:: Dh_coolant(prof)


   .. py:method:: P_thermal(prof)


   .. py:method:: P_coolant(prof)


   .. py:method:: compute_cross_section(prof, i)

      Build a closed OCC wire (int tag) for station i by:
      1) creating a circle EDGE in XY@origin,
      2) applying an affine transform to place it at (P_i, t_i, n_i, b_i),
      3) wrapping the transformed edge into a wire.
      NOTE: gmsh must be initialized and a model added by the caller.



   .. py:attribute:: n_points
      :value: 16



.. py:class:: CrossSectionRounded(n_points = 16)

   Bases: :py:obj:`ChannelSection`


   Common interface all cross-section shapes must implement.


   .. py:method:: A_coolant(prof)


   .. py:method:: Dh_coolant(prof)


   .. py:method:: P_thermal(prof)


   .. py:method:: P_coolant(prof)


   .. py:method:: compute_cross_section(prof, i)

      Build a closed OCC wire (int tag) for station i by:
      1) creating a circle EDGE in XY@origin,
      2) applying an affine transform to place it at (P_i, t_i, n_i, b_i),
      3) wrapping the transformed edge into a wire.
      NOTE: gmsh must be initialized and a model added by the caller.



   .. py:attribute:: n_points
      :value: 16



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



.. py:function:: make_channel_height_fn(contour, region_fractions, flat_heights, pinch_factors, transition_widths, logistic_k = 10.0)

   Like before, but now:
     - f = -1.0 maps to the chamber inlet,
     - f =  0.0 maps to the throat (min radius),
     - f = +1.0 maps to the nozzle exit.
   Values in-between interpolate linearly within the chamber or nozzle.


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



.. py:function:: get_theta_e_n(length_fraction, epsilon_value)

   This function takes in length fraction and area ratio, and two internally defined json files,
   theta_n.json and theta_e.json. These two files define the angles of the nozzle near the throat
   and the exit of the nozzle. The code interpolates twice for each angle in the length fraction,
   area ratio, theta space to yield the angles.

   :param length_fraction: The fractional length of the nozzle compared to a conical nozzle. A number between 0.60 and 1.00),
   :type length_fraction: float
   :param epsilon_value: Area ratio of the nozzle
   :type epsilon_value: float

   :returns:

             A tuple (theta_e, theta_n), both in radians.
                 theta_e (float): Interpolated exit angle in radians.
                 theta_n (float): Interpolated throat angle in radians.
   :rtype: tuple

   :returns:  1) For each bounding fraction dataset, interpolate across epsilon.
              2) Interpolate those results across the bounding fractions for the final answer.
   :rtype: (theta_e, theta_n) after performing two-stage 1D interpolation


.. py:function:: get_contour_internal(r_c, r_t, area_ratio, L_c, theta_conv, theta_div, nozzle, R_1f, R_2f, R_3f, length_fraction, export_tikz)

   Generate the full nozzle contour coordinates based on geometric parameters.

   This function calculates the x-coordinates and radii (ys) forming the nozzle contour
   for a thrust chamber. It performs several operations:
     - Computes the entrant and exit throat curves using the provided curvature factors.
     - Constructs the chamber contour, with or without a fillet depending on R_2f.
     - Constructs the nozzle contour based on the specified nozzle type ("rao" or "conical").
     - Concatenates all segments and then processes them to ensure the x-values are strictly increasing.

   :param r_c: Chamber radius.
   :type r_c: float
   :param r_t: Throat radius.
   :type r_t: float
   :param area_ratio: Nozzle area ratio.
   :type area_ratio: float
   :param L_c: Chamber length.
   :type L_c: float
   :param theta_conv: Convergence angle (in radians).
   :type theta_conv: float
   :param theta_div: Divergence angle (in radians).
   :type theta_div: float
   :param nozzle: Specifies the nozzle type, either "rao" or "conical".
   :type nozzle: str
   :param R_1f: Scaling factor for the throat entrant curvature radius.
   :type R_1f: float
   :param R_2f: Scaling factor for the chamber fillet radius. Use 0 or None for a hard corner.
   :type R_2f: float or None
   :param R_3f: Scaling factor for the throat exit curvature radius.
   :type R_3f: float
   :param length_fraction: A value between 0.60 and 1.00 used for interpolation.
   :type length_fraction: float

   :returns:

             A tuple (xs, ys) where:
                 xs (numpy.ndarray): Array of x-coordinates for the nozzle contour.
                 ys (numpy.ndarray): Array of corresponding radii for the nozzle contour.
   :rtype: tuple

   :raises ValueError: If the nozzle type is not 'rao' or 'conical', or if contour processing fails due to
       non-monotonic (non-increasing) x-values.


.. py:function:: compute_chamber_volume(xs, rs)

   Compute the chamber volume by revolving the contour around the x-axis.

   This function calculates the volume of the chamber by integrating the square of the radii
   (representing a circular cross-section) from the left boundary of the contour up to the throat,
   which is defined as the point with the minimum radius.

   :param xs: Sorted array of x-coordinates defining the contour (must be in ascending order).
   :type xs: array-like
   :param rs: Array of radii corresponding to the x-coordinates.
   :type rs: array-like

   :returns: The computed chamber volume, calculated as π times the integral of r² with respect to x.
   :rtype: float


.. py:function:: get_contour(r_t, area_ratio, r_c=None, L_c=None, V_c=None, eps_c=None, AR_c=None, theta_conv=45, theta_div=15, nozzle='rao', R_1f=1.5, R_2f=0.5, R_3f=0.382, length_fraction=0.8, angle_input='degrees', export_tikz=False)

   Generate the nozzle contour (xs, ys) using one of four valid input combinations.

   This function computes the nozzle contour for a thrust chamber using one of the following input methods:
     1. **Direct inputs**: Provide both chamber radius (r_c) and chamber length (L_c).
     2. **Volume & eps**: Provide chamber volume (V_c) and epsilon (eps_c). The chamber radius is computed from eps_c.
     3. **Volume & AR**: Provide chamber volume (V_c) and area ratio (AR_c). The chamber radius is computed
        from the relation r_c = (L_c * AR_c) / 2. # TODO: I need to update the aspect ratio definition to something more sensible
     4. **Volume & chamber radius**: Provide chamber volume (V_c) and chamber radius (r_c); L_c is determined by minimization.

   The input angles (theta_conv and theta_div) are expected in degrees if `angle_input` is "degrees"
   and are converted to radians internally.

   :param r_t: Throat radius.
   :type r_t: float
   :param area_ratio: Nozzle area ratio (epsilon).
   :type area_ratio: float
   :param r_c: Chamber radius.
   :type r_c: float, optional
   :param L_c: Chamber length.
   :type L_c: float, optional
   :param V_c: Chamber volume.
   :type V_c: float, optional
   :param eps_c: Epsilon value used to compute the chamber radius.
   :type eps_c: float, optional
   :param AR_c: Area ratio used with V_c to compute dimensions.
   :type AR_c: float, optional
   :param theta_conv: Convergence angle in degrees (default is 45).
   :type theta_conv: float, optional
   :param theta_div: Divergence angle in degrees (default is 15).
   :type theta_div: float, optional
   :param nozzle: Nozzle type; should be either "rao" or "conical" (default is "rao").
   :type nozzle: str, optional
   :param R_1f: Scaling factor for throat entrant curvature (default is 1.5).
   :type R_1f: float, optional
   :param R_2f: Scaling factor for chamber fillet curvature. Defaults to 0 if None.
   :type R_2f: float or None, optional
   :param R_3f: Scaling factor for throat exit curvature (default is 0.382).
   :type R_3f: float, optional
   :param length_fraction: A value between 0.60 and 1.00 used for interpolation (default is 0.8).
   :type length_fraction: float, optional
   :param angle_input: Unit for theta_conv and theta_div ("degrees" or "radians"). Default is "degrees".
   :type angle_input: str, optional

   :returns:

             A tuple (xs, ys) where:
                 xs (numpy.ndarray): Array of x-coordinates for the nozzle contour.
                 ys (numpy.ndarray): Array of corresponding radii for the nozzle contour.
   :rtype: tuple

   :raises ValueError: If the provided input combination is invalid or if minimization fails
       for calculating L_c.


.. py:function:: compute_cutoff_length(V_goal, xs_chamber, ys_chamber)

   Given a target volume V_goal, and arrays xs_chamber, ys_chamber
   (where xs_chamber runs from some negatives up through positives),
   return L_c = |x_cutoff| such that the volume of revolution about
   the x-axis from x=0 out to x=x_cutoff just reaches V_goal.

   :param V_goal: Desired volume (same units as π * ∫ y^2 dx).
   :type V_goal: float
   :param xs_chamber: x-coordinates, must cover the range from negative up to positive.
   :type xs_chamber: array_like, shape (N,)
   :param ys_chamber: y-values (assumed ≥0) corresponding to xs_chamber.
   :type ys_chamber: array_like, shape (N,)

   :returns: **L_c** -- The absolute distance |x_cutoff| from the origin where the
             cumulative volume first reaches V_goal.
   :rtype: float

   :raises ValueError: If V_goal is negative, or larger than the total volume available.


.. py:function:: integrate_area(L_c, xs_chamber, ys_chamber)

   Integrate the area under y(x) from x=0 down to x=-L_c.

   :param L_c: Positive cutoff length. The integration runs from x=0 to x=-L_c.
   :type L_c: float
   :param xs_chamber: x-coordinates, must cover the range from negative up through positive.
   :type xs_chamber: array_like, shape (N,)
   :param ys_chamber: y-values corresponding to xs_chamber.
   :type ys_chamber: array_like, shape (N,)

   :returns: **area** -- The area ∫_{0}^{-L_c} y(x) dx, returned as a positive number.
   :rtype: float

   :raises ValueError: If L_c is negative, or if -L_c lies outside the negative portion of xs_chamber.


