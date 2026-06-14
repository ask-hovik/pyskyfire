import numpy as np
import math
import pyskyfire.regen.physics as physics
from scipy.optimize import least_squares
import warnings
def short_warning(message, category, filename, lineno, file=None, line=None):
    return f"{category.__name__}: {message}\n"
warnings.formatwarning = short_warning

class BoundaryConditions:
    """Boundary conditions for the coolant-side inlet.

    Parameters
    ----------
    T_coolant_in : float
        Coolant static temperature at the channel inlet [K].
    p_coolant_in : float
        Coolant static pressure at the channel inlet [Pa].
    mdot_coolant : float
        Coolant mass-flow rate through the cooling channel [kg s⁻¹].
    """
    def __init__(self, T_coolant_in, p_coolant_in, mdot_coolant):
        self.T_coolant_in = T_coolant_in
        self.p_coolant_in = p_coolant_in
        self.mdot_coolant = mdot_coolant

# ================================================================
# Physics Class: Encapsulate the heat exchanger calculations
# ================================================================
class HeatExchangerPhysics:
    """Encapsulate hot-side, wall, and coolant heat-transfer/pressure models.

    This helper evaluates local heat transfer and pressure-loss terms for a
    given thrust-chamber and cooling circuit.

    Parameters
    ----------
    thrust_chamber : Any
        Object exposing geometry and property models used by the solver
        (e.g., ``contour``, ``combustion_transport``,
        ``wall_group``).
    circuit_index : int
        Index of the cooling circuit to use.
    """
    def __init__(self, thrust_chamber, boundary_conditions, circuit_index):
        self.thrust_chamber = thrust_chamber
        self.boundary_conditions = boundary_conditions
        self.circuit_index = circuit_index
        self.counter = 0

    def hot_side_coefficients(self, x, T_hw):
        """Hot-side heat input per unit length using Bartz-style correlation.

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        T_hw : float
            Hot-side wall temperature [K].

        Returns
        -------
        float
            ``dQ_hw/dx`` [W m⁻¹], positive when heating the wall.

        Notes
        -----
        Implements an **enthalpy-driven** Bartz form. The gas-side coefficient
        is evaluated with property data drawn from the chamber model at ``x``.
        """

        h_gas_corr = self.thrust_chamber.h_gas_corr # a user supplied correction factor
        D_hyd = 2*self.thrust_chamber.contour.r(x) # TODO: Need to update the hydraulic diameter to take the cooling channel shape into account?
        A_chmb = self.thrust_chamber.contour.A(x)
        mdot_g = self.thrust_chamber.combustion_transport.mdot
        T_g = self.thrust_chamber.combustion_transport.get_T(x)
        T_gr = (T_hw + T_g)/2 # film temperature

        # sometimes a nes simulation will reach this point without a guess for H_hw, in which case it's better to use 
        # H_g for an iteration rather than crashing the simulation. 
        H_g = self.thrust_chamber.combustion_transport.get_h(x)
        try:
            H_hw = self.thrust_chamber.combustion_transport.get_h(x, T=T_hw)
        except:
            H_hw = H_g

        M_g = self.thrust_chamber.combustion_transport.get_M(x)
        a_g = self.thrust_chamber.combustion_transport.get_a(x)
        H_gr = 0.5*(H_hw + H_g) + 0.18*(0.5*M_g**2*a_g**2) # original eq: H_gr = 0.5*(H_hw + H_g) + 0.18*(H_g0 - H_g)
        
        # the following units should be implemented using the H_gr enthalpy, and not just the x station
        # But the current aerothermodynamics module can't handle this. TODO: Implemene new aerothermodynamics
        # module to handle these queries properly
        Cp_gr = self.thrust_chamber.combustion_transport.get_cp(x)
        mu_gr = self.thrust_chamber.combustion_transport.get_mu(x)
        k_gr = self.thrust_chamber.combustion_transport.get_k(x)
        Pr_gr = self.thrust_chamber.combustion_transport.get_Pr(x)

        # The glorious Bartz equation
        h_gr = physics.h_gas_bartz_enthalpy_driven(k_gr, D_hyd, Cp_gr, mu_gr, mdot_g, A_chmb, T_g, T_gr)*h_gas_corr

        dA_dx_hot = self.thrust_chamber.cooling_circuits[self.circuit_index].dA_dx_thermal_exhaust(x)
        gamma = self.thrust_chamber.combustion_transport.get_gamma(x)
        T_aw = physics.T_aw(gamma=gamma, M_inf=M_g, T_inf=T_g, Pr=Pr_gr) # this equation was supposed to be used here for something, but I can't remember what atm. 
        H_aw = H_g + 0.5*Pr_gr**(1/3)*(M_g**2*a_g**2)

        h_g = h_gr/Cp_gr # enthalpy driven heat transfer definition

                # Existing enthalpy-based coefficient used by the solver

        # Corresponding heat flux
        qpp_hot = h_g * (H_aw - H_hw)   # [W m^-2]

        # Standard effective temperature-based coefficient for plotting/reporting
        h_hot = qpp_hot / (T_aw - T_hw)

        return {
            "h_hot": h_hot,
            "h_g": h_g,
            "h_gr": h_gr,
            "T_aw": T_aw,
            "qpp_hot": qpp_hot,
        }

    def cold_side_coefficients(self, x, T_cw, T_cool):
        """Coolant-side heat removal per unit length.

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        T_cw : float
            Coolant-side wall temperature [K].
        T_cool : float
            Bulk coolant temperature [K].

        Returns
        -------
        float
            ``dQ_cw/dx`` [W m⁻¹].
        """
        n_chan = self.thrust_chamber.cooling_circuits[self.circuit_index].placement.n_channel_positions*self.thrust_chamber.cooling_circuits[self.circuit_index].placement.n_channels_per_leaf

        p = self.thrust_chamber.combustion_transport.get_p(x)
        mdot_c_single_channel = self.boundary_conditions.mdot_coolant/n_chan
        T_coolant_film = (T_cool + T_cw)/2
        k_cf = self.thrust_chamber.cooling_circuits[self.circuit_index].coolant_transport.get_k(T_coolant_film, p)
        Cp_cr = self.thrust_chamber.cooling_circuits[self.circuit_index].coolant_transport.get_cp(T_coolant_film, p)
        mu_cf = self.thrust_chamber.cooling_circuits[self.circuit_index].coolant_transport.get_mu(T_coolant_film, p)

        D_c = self.thrust_chamber.cooling_circuits[self.circuit_index].Dh_coolant(x)
        A_channel = self.thrust_chamber.cooling_circuits[self.circuit_index].A_coolant(x)

        rho_bulk = self.thrust_chamber.cooling_circuits[self.circuit_index].coolant_transport.get_rho(T_cool, p)
        mu_bulk = self.thrust_chamber.cooling_circuits[self.circuit_index].coolant_transport.get_mu(T_cool, p)
        u_c = physics.u_coolant(rho_bulk, mdot_c_single_channel, A_channel)
        Re_c = physics.reynolds(rho_bulk, u_c, D_c, mu_bulk)

        R_curv = self.thrust_chamber.cooling_circuits[self.circuit_index].radius_of_curvature(x)
        phi_curv = physics.phi_curv(Re_c, D_c, R_curv)

        h_cold_corr = self.thrust_chamber.h_cold_corr
        h_c = physics.h_coolant_colburn(k_cf, D_c, Cp_cr, mu_cf, mdot_c_single_channel, A_channel, phi_curv=1)*h_cold_corr 
        

        return {
            "h_cold": h_c,
            "phi_curv": phi_curv,
            "Re_c": Re_c,
        }

    def dQ_hot_dx(self, x, T_hw):

        coeffs = self.hot_side_coefficients(x, T_hw)
        dA_dx_hot = self.thrust_chamber.cooling_circuits[self.circuit_index].dA_dx_thermal_exhaust(x)
        dQ_hw_dx = coeffs["qpp_hot"] * dA_dx_hot
        return dQ_hw_dx
    
    def dQ_cond_dx(self, x, T_hw, T_cw):
        """Conduction heat flow through the wall stack per unit length.

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        T_hw : float
            Hot-side wall temperature [K].
        T_cw : float
            Coolant-side wall temperature [K].

        Returns
        -------
        float
            ``dQ_cond/dx`` [W m⁻¹].

        Notes
        -----
        Treats each wall as a 1-D resistor in series:
        :math:`R_j = L_j / (k_j A)` with ``A = dA_dx_hot`` per unit length.
        """
                
        # 1) get the local hot‐side area per unit length
        dA_dx_hot = self.thrust_chamber.cooling_circuits[self.circuit_index].dA_dx_thermal_exhaust(x)

        # 2) sum each wall’s thermal resistance (R = L/(k·A)) per unit length
        walls = self.thrust_chamber.cooling_circuits[self.circuit_index].walls
        R_tot = sum(wall.thickness(x) / (wall.material.get_k((T_hw + T_cw)/2) * dA_dx_hot) for wall in walls)
        dQ_cond_dx = (T_hw - T_cw) / R_tot

        # 3) conduction per unit length
        return dQ_cond_dx

    def dQ_cold_dx(self, x, T_cw, T_cool):
        coeffs = self.cold_side_coefficients(x, T_cw, T_cool)
        h_c = coeffs["h_cold"]
        
        T_rep = 0.5*(T_cw + T_cool)  # or use T_cw, depends on preference
        R_cool_per_len = self.thrust_chamber.cooling_circuits[self.circuit_index].R_coolant_per_len(x, h_c=h_c, T_wall_rep=T_rep)

        dQ_cw_dx = (T_cw - T_cool) / R_cool_per_len
        
        return dQ_cw_dx

    def coolant_temperature_rate(self, T_cool, p_cool, dQ_cold_dx):
        """Axial rate of change of coolant temperature.

        Parameters
        ----------
        T_cool : float
            Coolant static temperature [K].
        p_cool : float
            Coolant static pressure [Pa].
        dQ_cold_dx : float
            Heat removed by coolant per unit length [W m⁻¹].

        Returns
        -------
        float
            ``dT_cool/dx`` [K m⁻¹].
        """
        n_chan = self.thrust_chamber.cooling_circuits[self.circuit_index].placement.n_channel_positions*self.thrust_chamber.cooling_circuits[self.circuit_index].placement.n_channels_per_leaf
        mdot_c = self.boundary_conditions.mdot_coolant/n_chan 
        cp = self.thrust_chamber.cooling_circuits[self.circuit_index].coolant_transport.get_cp(T_cool, p_cool)
        dT_dx = dQ_cold_dx/(mdot_c * cp) 
        
        return dT_dx # temperature change per unit length
    

    def coolant_pressure_rate(self, x, T_cool, p_cool):
        """Axial rates of static and stagnation pressure (friction + area change).

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        T_cool : float
            Coolant static temperature [K].
        p_cool : float
            Coolant static pressure [Pa].

        Returns
        -------
        tuple[float, float]
            ``(dp_static/dx, dp_stagnation/dx)`` in [Pa m⁻¹].
        """
        # get friction factor
        n_chan = self.thrust_chamber.cooling_circuits[self.circuit_index].placement.n_channel_positions*self.thrust_chamber.cooling_circuits[self.circuit_index].placement.n_channels_per_leaf
        rho_cool = self.thrust_chamber.cooling_circuits[self.circuit_index].coolant_transport.get_rho(T_cool, p_cool)
        mdot_c_single_channel = self.boundary_conditions.mdot_coolant/n_chan

        A_cool = self.thrust_chamber.cooling_circuits[self.circuit_index].A_coolant(x)
        u_cool = physics.u_coolant(rho_cool, mdot_c_single_channel, A_cool) 

        D_h = self.thrust_chamber.cooling_circuits[self.circuit_index].Dh_coolant(x)
        mu = self.thrust_chamber.cooling_circuits[self.circuit_index].coolant_transport.get_mu(T_cool, p_cool)
        roughness = self.thrust_chamber.cooling_circuits[self.circuit_index].roughness
        ReDh = physics.reynolds(rho_cool, u_cool, D_h, mu)
        f = physics.f_darcy(ReDh, D_h, x, roughness)

        # equivalent length due to curvature: 
        R_curv = self.thrust_chamber.cooling_circuits[self.circuit_index].radius_of_curvature(x)
        K = self.thrust_chamber.K_factor  
        dL_eq_dx = (2 * K * D_h) / (np.pi * R_curv * f) # Differential equivalent length per unit dx
        curvature_factor = 1.0 + dL_eq_dx

        circuit = self.thrust_chamber.cooling_circuits[self.circuit_index]
        dp_stagnation_dx = - f/D_h * rho_cool*u_cool**2/2 * circuit.ds_dx(x) # replaced the slope factor with this

        dA_dx = self.thrust_chamber.cooling_circuits[self.circuit_index].dA_dx_coolant(x) 

        dp_static_dx = dp_stagnation_dx - rho_cool*u_cool**2/A_cool*dA_dx

        return dp_static_dx, dp_stagnation_dx
    
    def interface_temperatures(self, x, T_hw, T_cw):
        """Wall-interface temperatures across the stack at position ``x``.

        Parameters
        ----------
        x : float
            Axial coordinate [m].
        T_hw : float
            Hot-side wall temperature [K].
        T_cw : float
            Coolant-side wall temperature [K].

        Returns
        -------
        list[float]
            Temperatures ``[T_hot, T_1, ..., T_cold]`` across interfaces.
        """
        # same area-per-length as above
        dA_dx_hot = self.thrust_chamber.cooling_circuits[self.circuit_index].dA_dx_thermal_exhaust(x)

        # heat transfer per length
        qdx = self.dQ_cond_dx(x, T_hw, T_cw)
        T_rep = 0.5 * (T_hw + T_cw)
        Ts = [T_hw]
        for wall in self.thrust_chamber.cooling_circuits[self.circuit_index].walls:
            Rj = wall.thickness(x) / (wall.material.get_k(T_rep) * dA_dx_hot)
            # drop in temperature across this layer:
            T_next = Ts[-1] - qdx * Rj
            Ts.append(T_next)
        # Ts[-1] should equal T_cw (within numerical roundoff)
        return Ts

    

def solve_heat_exchanger(thrust_chamber, boundary_conditions, n_nodes, circuit_index, output, log_residuals=True):
    """Solve 1-D steady heating with a marching scheme.

    Parameters
    ----------
    thrust_chamber : Any
        Chamber model exposing geometry and property methods.
    boundary_conditions : BoundaryConditions
        Inlet temperature/pressure/mass-flow conditions.
    n_nodes : int
        Number of axial nodes.
    circuit_index : int
        Which cooling circuit to simulate.
    output : bool
        If True, print progress to stdout.
    log_residuals : bool, optional
        If True, record local residuals each Newton iteration per cell.

    Returns
    -------
    dict
        Results with keys:

        - ``x`` : axial coordinates [m]
        - ``T`` : temperatures array, shape ``(n_nodes, 1 + n_walls + 1)``
          (coolant + reversed wall interfaces from cold→hot)
        - ``T_static`` : coolant static temperature [K]
        - ``T_stagnation`` : coolant stagnation temperature [K]
        - ``p_static`` : coolant static pressure [Pa]
        - ``p_stagnation`` : coolant stagnation pressure [Pa]
        - ``dQ_dA`` : local heat flux [W m⁻²]
        - ``velocity`` : coolant velocity [m s⁻¹]
        - ``residuals`` : tuple of (global history, final-per-cell) or (None, None)
    """

    # 1) Build the axial grid in [x_min, x_max].
    residual_log = [] if log_residuals else None
    iter_counter = np.zeros(n_nodes, dtype=int)
    
    circuit = thrust_chamber.cooling_circuits[circuit_index]
    orig_x_domain = circuit.x_domain
    # Re-interpolate the x_domain to have exactly n_nodes points.
    if circuit.direction == 1:
        x_domain = np.linspace(orig_x_domain[0], orig_x_domain[-1], n_nodes)
    else: 
        x_domain = np.linspace(orig_x_domain[-1], orig_x_domain[0], n_nodes)
    # Compute the spacing based on the x_domain.
    dx = abs((x_domain[-1] - x_domain[0]) / (n_nodes - 1))
    
    # 2) Prepare arrays to hold results
    T_hw_arr    = np.zeros(n_nodes)  # wall on hot side
    T_cw_arr    = np.zeros(n_nodes)  # wall on coolant side
    T_cool_arr  = np.zeros(n_nodes)  # coolant temperature
    p_static_arr  = np.zeros(n_nodes)  # coolant pressure
    p_stagnation_arr = np.zeros(n_nodes)
    dQ_dA_arr   = np.zeros(n_nodes)  # local heat flux (W/m²)

    # Inlet boundary conditions
    T_cool_in = boundary_conditions.T_coolant_in
    p_cool_in = boundary_conditions.p_coolant_in
    T_cool_arr[0] = T_cool_in
    p_static_arr[0] = p_cool_in
    p_stagnation_arr[0] = p_cool_in # assuming here no velocity at inlet.....
    # Physics helper
    physics_helper = HeatExchangerPhysics(thrust_chamber, boundary_conditions, circuit_index)

    # Initial guesses for the wall temperatures
    T_hw_guess = 0.5 * (thrust_chamber.combustion_transport.get_T(x_domain[0]) + T_cool_in)
    T_cw_guess = 0.5 * (thrust_chamber.combustion_transport.get_T(x_domain[0]) + T_cool_in)

    # March along x (simulation loop)
    if output is True:
        print(f"Started heat exchanger simulation with {n_nodes} nodes")


    # ==============
    # --- SOLVER ---
    # ==============
    for i in range(n_nodes):
        if output is True:
            print(f'\rSimulating: {math.ceil(i/n_nodes*100)}%', end='', flush=True)

        x_i = x_domain[i]
        
        #all_T_interfaces = physics_helper.interface_temperatures(x_i, T_hw_sol, T_cw_sol)
        # Solve for (T_cw, dT) using bounded, robust least-squares.
        # Enforces T_hw = T_cw + dT with dT >= 0  =>  T_hw >= T_cw.

        T_gas_i = physics_helper.thrust_chamber.combustion_transport.get_T(x_i)

        # Conservative upper bound for wall temperatures to avoid sign flips in (T_gas - T_hw)
        # and keep the solver in a physical basin.
        T_hw_max = max(T_cool_arr[i] + 1.0, T_gas_i - 1e-3)

        def residuals_scaled(vars_, cell=i):
            T_cw_trial, dT_trial = vars_
            T_hw_trial = T_cw_trial + dT_trial

            # Use per-cell energy flows
            Q_hot_val  = physics_helper.dQ_hot_dx(x_i, T_hw_trial) * dx
            Q_cond_val = physics_helper.dQ_cond_dx(x_i, T_hw_trial, T_cw_trial) * dx
            Q_cold_val = physics_helper.dQ_cold_dx(x_i, T_cw_trial, T_cool_arr[i]) * dx

            R1 = Q_hot_val - Q_cond_val
            R2 = Q_cond_val - Q_cold_val

            # Optional logging of *raw* residuals (not scaled)
            if residual_log is not None:
                k = iter_counter[cell]
                residual_log.append((cell, k, R1, R2))
                iter_counter[cell] += 1

            # Scale residuals to reduce conditioning issues at extreme heat loads
            Q_ref = max(abs(Q_hot_val), abs(Q_cond_val), abs(Q_cold_val), 1.0)
            return np.array([R1 / Q_ref, R2 / Q_ref], dtype=float)

        # Initial guess in (T_cw, dT) space
        dT_guess = max(T_hw_guess - T_cw_guess, 1.0)
        x0 = np.array([T_cw_guess, dT_guess], dtype=float)

        # Bounds:
        # - T_cw cannot go below bulk coolant temp at this station.
        # - T_hw cannot exceed T_hw_max, so dT <= T_hw_max - T_cw_lo.
        T_cw_lo = T_cool_arr[i]
        T_cw_hi = T_hw_max
        dT_lo = 0.0
        dT_hi = max(1.0, T_hw_max - T_cw_lo)

        res = least_squares(
            residuals_scaled,
            x0=x0,
            bounds=([T_cw_lo, dT_lo], [T_cw_hi, dT_hi]),
            method="trf",
            loss="soft_l1",
            f_scale=1.0,
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
            max_nfev=200,
        )

        if not res.success:
            raise RuntimeError(
                f"least_squares failed at node {i} (x={x_i:.6g}): {res.message}; "
                f"x={res.x}"
            )

        T_cw_sol, dT_sol = res.x
        T_hw_sol = T_cw_sol + dT_sol

        # Store in arrays
        T_hw_arr[i] = T_hw_sol
        T_cw_arr[i] = T_cw_sol

        # Update guesses for the next node
        T_hw_guess, T_cw_guess = T_hw_sol, T_cw_sol

        # Update coolant temperature and pressure for the next node
        if i < n_nodes - 1:
            Q_cold_val = physics_helper.dQ_cold_dx(x_i, T_cw_sol, T_cool_arr[i])*dx
            dT = physics_helper.coolant_temperature_rate(T_cool_arr[i], p_stagnation_arr[i], Q_cold_val)
            dp_static, dp_stagnation = physics_helper.coolant_pressure_rate(x_i, T_cool_arr[i], p_stagnation_arr[i])
            dp_static = dp_static*dx
            dp_stagnation = dp_stagnation*dx
            T_cool_arr[i+1] = T_cool_arr[i] + dT
            p_static_arr[i+1] = p_static_arr[i] + dp_static
            p_stagnation_arr[i+1] = p_stagnation_arr[i] + dp_stagnation

    if output is True:
        print(f'\rSimulating: {100}%\n', end='', flush=True)

    # ===============
    # --- RESULTS ---
    # ===============

    # 1. Compute the local heat flux (dQ/dA) for each node 
    for i, x_i in enumerate(x_domain):
        Q_hot_val = physics_helper.dQ_hot_dx(x_i, T_hw_arr[i]) * dx
        A_hot = thrust_chamber.cooling_circuits[circuit_index].dA_dx_thermal_exhaust(x_i) * dx
        dQ_dA_arr[i] = Q_hot_val / A_hot if A_hot != 0 else 0.0

    # 2. Compute the coolant velocity at each node
    velocity_arr = np.zeros(n_nodes)
    n_chan = thrust_chamber.cooling_circuits[circuit_index].placement.n_channel_positions*thrust_chamber.cooling_circuits[circuit_index].placement.n_channels_per_leaf
    for i, x_i in enumerate(x_domain):
        A_channel = thrust_chamber.cooling_circuits[circuit_index].A_coolant(x_i)
        mdot_c_single = boundary_conditions.mdot_coolant / n_chan
        rho_cool = thrust_chamber.cooling_circuits[circuit_index].coolant_transport.get_rho(T_cool_arr[i], p_static_arr[i])
        # Use the physics helper's u_coolant function 
        u_cool = physics.u_coolant(rho_cool, mdot_c_single, A_channel)
        velocity_arr[i] = u_cool
    
    T_stagnation_arr = np.zeros_like(T_cool_arr)
    for i, x_i in enumerate(x_domain):
        # Calculate cp for the current node using static temperature and pressure
        cp_cool = thrust_chamber.cooling_circuits[circuit_index].coolant_transport.get_cp(T_cool_arr[i], p_static_arr[i])
        # Compute stagnation temperature by adding the kinetic energy term
        T_stagnation_arr[i] = T_cool_arr[i] + (velocity_arr[i]**2) / (2.0 * cp_cool)

    n_walls = len(thrust_chamber.cooling_circuits[circuit_index].walls)
    T_full = np.zeros((n_nodes, 1 + n_walls + 1))  # coolant + (n_walls+1) interfaces

    for i, x_i in enumerate(x_domain):
        # Ts = [T_hot, ..., T_cold], length = n_walls+1
        Ts = physics_helper.interface_temperatures(x_i,
                                                   T_hw_arr[i],
                                                   T_cw_arr[i]
                                                   )
        # reverse so Ts[0]=T_cold, Ts[-1]=T_hot
        Ts_rev = Ts[::-1]
        T_full[i, 0]    = T_cool_arr[i]
        T_full[i, 1:]   = Ts_rev

    p_static_corrected = np.zeros_like(p_stagnation_arr)
    for i in range(n_nodes):
        # use density at node i (assumes rho(T,p) weakly depends on p)
        rho_i = thrust_chamber.cooling_circuits[circuit_index] \
                   .coolant_transport.get_rho(T_cool_arr[i], p_stagnation_arr[i])
        q_dyn = 0.5 * rho_i * velocity_arr[i]**2
        p_static_corrected[i] = p_stagnation_arr[i] - q_dyn

    # overwrite the old (bad) static-pressure array
    p_static_arr = p_static_corrected

    # 1. Compute hot- and cold-side heat transfer coefficients at each node
    h_hot_arr = np.zeros(n_nodes)            # effective hot-side h [W/m²/K]
    h_hot_enthalpy_arr = np.zeros(n_nodes)   # enthalpy-based hot-side coeff [kg/m²/s]
    h_cold_arr = np.zeros(n_nodes)           # coolant-side h [W/m²/K]
    T_aw_hot_arr = np.zeros(n_nodes)         # adiabatic wall temperature [K]
    for i, x_i in enumerate(x_domain):
        hot = physics_helper.hot_side_coefficients(x_i, T_hw_arr[i])
        cold = physics_helper.cold_side_coefficients(x_i, T_cw_arr[i], T_cool_arr[i])

        h_hot_arr[i] = hot["h_hot"]
        h_hot_enthalpy_arr[i] = hot["h_g"]
        h_cold_arr[i] = cold["h_cold"]
        T_aw_hot_arr[i] = hot["T_aw"]

    global_R, final_R = analyse_residuals(residual_log, n_nodes)
    #print(f"TP was acessed {physics_helper.counter} times")
    cooling_data = {
        "x"             : x_domain,
        "T"             : T_full,
        "T_static"      : T_cool_arr,
        "T_stagnation"  : T_stagnation_arr,
        "p_static"      : p_static_arr,
        "p_stagnation"  : p_stagnation_arr,
        "dQ_dA"         : dQ_dA_arr,
        "velocity"      : velocity_arr,
        "h_hot"         : h_hot_arr,
        "h_hot_enthalpy": h_hot_enthalpy_arr,
        "h_cold"        : h_cold_arr,
        "T_aw_hot"      : T_aw_hot_arr,
        "residuals"     : (global_R, final_R)
    }

    return cooling_data

def analyse_residuals(residual_log, n_cells, p=2):
    """Aggregate local solver residuals into global history and final per-cell vector.

    Parameters
    ----------
    residual_log : list or None
        List of tuples ``(cell, iter, R1, R2)`` recorded during solves.
        If ``None`` or empty, returns ``(None, None)``.
    n_cells : int
        Number of axial cells.
    p : int or float, optional
        Norm order for global residual history: ``2`` for RMS,
        ``np.inf`` for :math:`L_\\infty`, etc. Default is 2.

    Returns
    -------
    history : ndarray or None
        Global residual norm for iterations ``0..k_max``, or ``None``.
    final_per_cell : ndarray or None
        Final residual magnitude per cell at its last local iteration, or ``None``.
    """
    if not residual_log:                          # catches [] and None
        return None, None

    log = np.asarray(residual_log)                # shape (m,4)
    cells = log[:, 0].astype(int)
    iters = log[:, 1].astype(int)
    rmag  = np.hypot(log[:, 2], log[:, 3])        # L2 of (R1,R2)

    # ---- global norm history  -----------------------------
    k_max = iters.max()
    history = np.empty(k_max + 1)

    for k in range(k_max + 1):
        mask = iters == k
        if p == np.inf:
            history[k] = rmag[mask].max()
        else:
            history[k] = (rmag[mask] ** p).mean() ** (1.0 / p)

    # ---- per‑cell residual at final local iteration -------
    final_per_cell = np.full(n_cells, np.nan)
    for c in range(n_cells):
        mask = cells == c
        if mask.any():
            final_per_cell[c] = rmag[mask][-1]    # last entry for cell c

    return history, final_per_cell



def steady_heating_analysis(thrust_chamber, boundary_conditions, n_nodes=100, circuit_index=0, solver="newton", output=True):
    """Run the steady heating analysis.

    Parameters
    ----------
    thrust_chamber : Any
        Chamber model exposing the required geometry & property APIs.
    boundary_conditions : BoundaryConditions
        Coolant inlet boundary conditions.
    n_nodes : int, optional
        Number of axial nodes. Default is 100.
    circuit_index : int, optional
        Cooling-circuit index. Default is 0.
    solver : {'newton'}, optional
        Solver selector. Currently only ``'newton'`` is implemented.
    output : bool, optional
        If True, print progress. Default is True.

    Returns
    -------
    dict
        See :func:`solve_heat_exchanger` for keys.
    """
    if solver.lower() == "newton":
        return solve_heat_exchanger(thrust_chamber, boundary_conditions, n_nodes, circuit_index, output)
    else: 
        print("solver name not recognized")
    # possibility to implement other solvers here
