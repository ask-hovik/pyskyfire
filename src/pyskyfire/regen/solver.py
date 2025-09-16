import numpy as np
import math
import pyskyfire.regen.physics as physics
from scipy.optimize import fsolve

class BoundaryConditions:
    """Object for storing boundary conditions for the solver.
     
    Args: 
        T_coolant_in (float): Static(?) temperature of coolant at cooling channel inlet (K)
        p_coolant_in (float): Static(?) pressure of coolant at cooling channel inlet (Pa)
        mdot_coolant (float): mass flow rate of coolant through cooling channel inlet (kg/s)
    """
    def __init__(self, T_coolant_in, p_coolant_in, mdot_coolant):
        self.T_coolant_in = T_coolant_in
        self.p_coolant_in = p_coolant_in
        self.mdot_coolant = mdot_coolant

# ================================================================
# Physics Class: Encapsulate the heat exchanger calculations
# ================================================================
class HeatExchangerPhysics:
    """
    Encapsulates the physical calculations for the heat exchanger.
    This class is responsible for computing the heat fluxes and rates based
    on the given engine properties and operating conditions.
    """
    def __init__(self, thrust_chamber, circuit_index):
        self.thrust_chamber = thrust_chamber
        self.circuit_index = circuit_index
        self.counter = 0

    def dQ_hot_dx(self, x, T_hw):
        """
        Computes the heat transfer rate from the hot gas to the wall. Bartz equation, for example
        """
        approach = "new"
        if approach == "new":
            h_gas_corr = self.thrust_chamber.h_gas_corr # a user supplied correction factor
            D_hyd = 2*self.thrust_chamber.contour.r(x) # TODO: Need to update the hydraulic diameter to take the cooling channel shape into account?
            A_chmb = self.thrust_chamber.contour.A(x)
            mdot_g = self.thrust_chamber.combustion_transport.mdot
            T_g = self.thrust_chamber.combustion_transport.get_T(x)
            T_gr = (T_hw + T_g)/2 # film temperature

            H_hw = self.thrust_chamber.combustion_transport.get_h(x, T=T_hw)
            self.counter +=1
            
            H_g = self.thrust_chamber.combustion_transport.get_h(x)

            M_g = self.thrust_chamber.combustion_transport.get_M(x)
            a_g = self.thrust_chamber.combustion_transport.get_a(x)
            H_gr = 0.5*(H_hw + H_g) + 0.18*(0.5*M_g**2*a_g**2) # original eq: H_gr = 0.5*(H_hw + H_g) + 0.18*(H_g0 - H_g)
            
            #print(f"H_gr:{H_gr}")
            #input()
            #Cp_gr = self.thrust_chamber.combustion_transport.get_cp(x, T=T_gr)
            #mu_gr = self.thrust_chamber.combustion_transport.get_mu(x, T=T_gr)
            #k_gr = self.thrust_chamber.combustion_transport.get_k(x, T=T_gr)
            #Pr_gr = Cp_gr*mu_gr/k_gr

            #props = self.thrust_chamber.combustion_transport.get_properties(x, h=H_gr, props=("cp","mu","k","Pr"))
            #Cp_gr, mu_gr, k_gr, Pr_gr = props["cp"], props["mu"], props["k"], props["Pr"]

            # could potentially evaluate all of these in a single run (much faster)
            #Cp_gr = self.thrust_chamber.combustion_transport.get_cp(x, h=H_gr)
            #mu_gr = self.thrust_chamber.combustion_transport.get_mu(x, h=H_gr)
            #k_gr = self.thrust_chamber.combustion_transport.get_k(x, h=H_gr)
            #Pr_gr = Cp_gr*mu_gr/k_gr

            Cp_gr = self.thrust_chamber.combustion_transport.get_cp(x)
            mu_gr = self.thrust_chamber.combustion_transport.get_mu(x)
            k_gr = self.thrust_chamber.combustion_transport.get_k(x)
            Pr_gr = self.thrust_chamber.combustion_transport.get_Pr(x)


            
            # The glorious Bartz equation
            h_gr = physics.h_gas_bartz_enthalpy_driven(k_gr, D_hyd, Cp_gr, mu_gr, mdot_g, A_chmb, T_g, T_gr)*h_gas_corr

            dA_dx_hot = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].dA_dx_thermal_exhaust(x)
            gamma = self.thrust_chamber.combustion_transport.get_gamma(x)
            T_aw = physics.T_aw(gamma=gamma, M_inf=M_g, T_inf=T_g, Pr=Pr_gr)
            H_aw = H_g + 0.5*Pr_gr**(1/3)*(M_g**2*a_g**2)
            #print(H_aw)
            #input()
            #H_aw = -1399116.7949471772

            #print(f"H_aw: {H_aw}")

            h_g = h_gr/Cp_gr # enthalpy driven heat transfer definition
            dQ_hw_dx = h_g*dA_dx_hot*(H_aw - H_hw)

        elif approach == "old":
            # Old Approach:
            D_hyd = 2*self.thrust_chamber.contour.r(x) # TODO: Need to update the hydraulic diameter to take the cooling channel shape into account?
            A_chmb = self.thrust_chamber.contour.A(x)
            mdot_g = self.thrust_chamber.combustion_transport.mdot
            T_g = self.thrust_chamber.combustion_transport.get_T(x)
            T_gr = (T_hw + T_g)/2 # film temperature
            
            Cp_gr = self.thrust_chamber.combustion_transport.get_cp(x)
            mu_gr = self.thrust_chamber.combustion_transport.get_mu(x)
            k_gr = self.thrust_chamber.combustion_transport.get_k(x)
            h_gas_corr = self.thrust_chamber.h_gas_corr # a user supplied correction factor

            h_gr = physics.h_gas_bartz_enthalpy_driven(k_gr, D_hyd, Cp_gr, mu_gr, mdot_g, A_chmb, T_g, T_gr)*h_gas_corr
            dA_dx_hot = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].dA_dx_thermal_exhaust(x)
            M_inf = self.thrust_chamber.combustion_transport.get_M(x)
            gamma = self.thrust_chamber.combustion_transport.get_gamma(x)
            T_aw = self.thrust_chamber.combustion_transport.get_T_aw(x, gamma, M_inf, T_inf=T_g)
            #H_aw = self.thrust_chamber.combustion_transport.get_h(x, T=T_aw)
            H_aw = -1399116.7949471772
            #H_hw = self.thrust_chamber.combustion_transport.get_h(x, T=T_hw) # problem child
            H_hw = -11247569.350990318
            h_g = h_gr/Cp_gr # enthalpy driven heat transfer definition
            dQ_hw_dx = h_g*dA_dx_hot*(H_aw - H_hw)

        elif approach == "idea":
            
            D_t = self.thrust_chamber.contour.r_t*2
            cp_g = self.thrust_chamber.combustion_transport.get_cp(x)
            mu_g = self.thrust_chamber.combustion_transport.get_mu(x)
            Pr_g = self.thrust_chamber.combustion_transport.get_Pr(x)
            p_c = self.thrust_chamber.combustion_transport.p_c
            c_star = self.thrust_chamber.combustion_transport.c_star
            A_t = self.thrust_chamber.contour.A_t
            A_x = self.thrust_chamber.contour.A(x)
            T_c = self.thrust_chamber.combustion_transport.T_c
            gamma_g = self.thrust_chamber.combustion_transport.get_gamma(x)
            M_g = self.thrust_chamber.combustion_transport.get_M(x)
            
            sigma = physics.sigma(T_hw, T_c, gamma_g, M_g, omega=0.6)
            h_g = physics.h_gas_bartz(D_t, mu_g, cp_g, Pr_g, p_c, c_star, A_t, A_x, sigma)

            dA_dx_hot = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].dA_dx_thermal_exhaust(x)
            T_g = self.thrust_chamber.combustion_transport.get_T(x)
            dQ_hw_dx = h_g*dA_dx_hot*(T_g - T_hw) 


        printout = False
        if printout: 
            print(f"D_hyd in Q_hot: {D_hyd}")
            print(f"A_chmb in Q_hot: {A_chmb}")
            print(f"mdot in Q_hot: {mdot_g}")
            print(f"T_g in Q_hot: {T_g}")
            print(f"T_gr in Q_hot: {T_gr}")
            print(f"Cp_gr in Q_hot: {Cp_gr}")
            print(f"mu_gr in Q_hot: {mu_gr}")
            print(f"k_gr in Q_hot: {k_gr}")
            print(f"h_gr in Q_hot: {h_gr}")
            print(f"dA_dx_hot in Q_hot: {dA_dx_hot}")
            #print(f"p in Q_hot: {p}")
            print(f"H_hw in Q_hot: {H_hw}")
            #print(f"M_inf in Q_hot: {M_g}")
            print(f"gamma in Q_hot: {gamma}")
            print(f"T_aw in Q_hot: {T_aw}")
            print(f"H_aw in Q_hot: {H_aw}")
            print(f"dQ_hw_dx in Q_hot: {dQ_hw_dx}\n")
            x = input()

        return dQ_hw_dx

    """def dQ_cond_dx(self, x, T_hw, T_cw):
        
        Computes the conductive heat transfer through the wall.
        
        k = self.engine.thrust_chamber.wall_group.walls[0].material.k # TODO: update to handle more walls. 

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        dA_dx_hot = self.engine.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].dA_dx_thermal_exhaust(x)

        L = self.engine.thrust_chamber.wall_group.walls[0].thickness(x)
        # k * A * (T_hw - T_cw)
        dQ_cond_dx = k*dA_dx_hot*(T_hw - T_cw)/L 
        
        printout = False
        if printout: 
            print(f"k in Q_cond: {k}")
            print(f"dA_dx_hot in Q_cond: {dA_dx_hot}")
            print(f"L in Q_cond: {L}")
            print(f"Q_cond in Q_cond: {dQ_cond_dx}\n")
            print(f"T_hw in Q_cond: {T_hw}")
            print(f"T_cw in Q_cond: {T_cw}")

        return dQ_cond_dx"""
    
    def dQ_cond_dx(self, x, T_hw, T_cw):
        # 1) get the local hot‐side area per unit length
        dA_dx_hot = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].dA_dx_thermal_exhaust(x)

        # 2) sum each wall’s thermal resistance (R = L/(k·A)) per unit length
        walls = self.thrust_chamber.wall_group.walls
        R_tot = sum(wall.thickness(x) / (wall.material.get_k((T_hw + T_cw)/2) * dA_dx_hot) for wall in walls) # TODO: this looks a little iffy
        dQ_cond_dx = (T_hw - T_cw) / R_tot

        # 3) conduction per unit length
        return dQ_cond_dx

    def dQ_cold_dx(self, x, T_cw, T_cool):
        """
        Computes the heat transfer rate from the wall to the coolant.
        """
        n_chan = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].placement.n_channel_positions*self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].placement.n_channels_per_leaf

        p = self.thrust_chamber.combustion_transport.get_p(x)
        mdot_c_single_channel = self.thrust_chamber.combustion_transport.mdot_fu/n_chan # TODO: uh, this is kinda absolutely wrong should pull mdot from the cooling circuit
        T_coolant_film = (T_cool + T_cw)/2
        k_cf = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].coolant_transport.get_k(T_coolant_film, p)
        Cp_cr = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].coolant_transport.get_cp(T_coolant_film, p)
        mu_cf = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].coolant_transport.get_mu(T_coolant_film, p)

        D_c = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].Dh_coolant(x)
        dA_dx_cool = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].dA_dx_thermal_coolant(x) # TODO: reimplement fin area
        A_channel = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].A_coolant(x)

        rho_bulk = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].coolant_transport.get_rho(T_cool, p)
        mu_bulk = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].coolant_transport.get_mu(T_cool, p)
        u_c = physics.u_coolant(rho_bulk, mdot_c_single_channel, A_channel)
        Re_c = physics.reynolds(rho_bulk, u_c, D_c, mu_bulk)

        R_curv = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].radius_of_curvature(x)
        phi_curv = physics.phi_curv(Re_c, D_c, R_curv)

        h_cold_corr = self.thrust_chamber.h_cold_corr
        h_c = physics.h_coolant_colburn(k_cf, D_c, Cp_cr, mu_cf, mdot_c_single_channel, A_channel, phi_curv=1)*h_cold_corr 
        dQ_cw_dx = h_c*dA_dx_cool*(T_cw - T_cool) 
        
        printout = False
        if printout: 
            print(f"T_cw in Q_cold: {T_cw}")
            print(f"T_cool in Q_cold: {T_cool}")
            print(f"p in Q_cold: {p}")
            print(f"mdot_single_chan in Q_cond: {mdot_c_single_channel}")
            print(f"T_coolant_film in Q_cond: {T_coolant_film}")
            print(f"k_cf in Q_cold: {k_cf}")
            print(f"Cp_cr in Q_cold: {Cp_cr}")
            print(f"mu_cf in Q_cold: {mu_cf}")
            print(f"D_c in Q_cold: {D_c}")
            print(f"A_c in Q_cold: {dA_dx_cool}")
            print(f"h_c in Q_cold: {h_c}")
            print(f"Q_cw in Q_cold: {dQ_cw_dx}\n")
        return dQ_cw_dx

    def coolant_temperature_rate(self, T_cool, p_cool, dQ_cold_dx):
        """
        Computes the rate of change of the coolant temperature.
        """
        n_chan = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].placement.n_channel_positions*self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].placement.n_channels_per_leaf
        mdot_c = self.thrust_chamber.combustion_transport.mdot_fu/n_chan # TODO: this is not actually the coolant flow, its the combustion reaction coolant flow!
        cp = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].coolant_transport.get_cp(T_cool, p_cool)
        dT_dx = dQ_cold_dx/(mdot_c * cp) 
        
        printout = False
        if printout: 
            print(f"\nQ_cold in dTdx: {dQ_cold_dx}")
            print(f"T_cool in dTdx: {T_cool}")
            print(f"p_cool in dTdx: {p_cool}")
            print(f"mdot_coolant in dTdx: {mdot_c}")
            print(f"cp in dTdx: {cp}")
            print(f"dTdx in dTdx: {dT_dx}\n")

        return dT_dx # temperature change per unit length
    

    def coolant_pressure_rate(self, x, T_cool, p_cool):
        """
        Computes the rate of change of the coolant pressure due to frictional losses.
        """
        # get friction factor
        n_chan = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].placement.n_channel_positions*self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].placement.n_channels_per_leaf
        rho_cool = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].coolant_transport.get_rho(T_cool, p_cool)
        mdot_c_single_channel = self.thrust_chamber.combustion_transport.mdot_fu/n_chan

        A_cool = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].A_coolant(x)
        u_cool = physics.u_coolant(rho_cool, mdot_c_single_channel, A_cool) 

        D_h = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].Dh_coolant(x)
        mu = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].coolant_transport.get_mu(T_cool, p_cool)
        roughness = self.thrust_chamber.roughness
        ReDh = physics.reynolds(rho_cool, u_cool, D_h, mu)
        f = physics.f_darcy(ReDh, D_h, x, roughness)

        # equivalent length due to curvature: 
        R_curv = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].radius_of_curvature(x)
        K = self.thrust_chamber.K_factor  
        dL_eq_dx = (2 * K * D_h) / (np.pi * R_curv * f) # Differential equivalent length per unit dx
        curvature_factor = 1.0 + dL_eq_dx
        
        # actual length since the thrust chamber is sloping
        dr_dx = self.thrust_chamber.contour.dr_dx(x)
        slope_factor = np.sqrt(1 + dr_dx**2) # TODO: this expression only works for vertical cooling channels, need to update for helical.
        #slope_factor = 1
        # pressure drop per unit length
        dp_stagnation_dx = - f/D_h * rho_cool*u_cool**2/2*slope_factor#*curvature_factor
        """print(f"\nrho_cool: {rho_cool}")
        print(f"u_cool: {u_cool}")
        print(f"f: {f}")
        print(f"D_h: {D_h}")
        print(f"slope_factor: {slope_factor} \n")"""

        # Remove dynamic pressure from stagnation pressure

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        dA_dx = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].dA_dx_coolant(x) 
        #if self.engine.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].direction != 1:
        #    dA_dx = -dA_dx

        #dp_dynamic_dx = - 0.5*rho_cool*u_cool**2/A_cool*dA_dx

        dp_static_dx = dp_stagnation_dx - rho_cool*u_cool**2/A_cool*dA_dx

        """print(f"\nx: {x}")
        print(f"rho_cool: {rho_cool}")
        print(f"u_cool: {u_cool}")
        print(f"A_cool: {A_cool}")
        print(f"dA_dx: {dA_dx} \n")"""

        printout = False
        if printout: 
            print(f"n_chan: {n_chan}")
            print(f"rho_cool in dpdx: {rho_cool}")
            print(f"mdot_coolant_single_cannel in dpdx: {mdot_c_single_channel}")
            print(f"A_cool in dpdx: {A_cool}")
            print(f"u_cool in dpdx: {u_cool}")
            print(f"darcy_friction in dpdx: {f}")
            print(f"Reynolds in dpdx: {ReDh}")
            print(f"roughness in dpdx: {roughness}")
            print(f"curvature_factor: {curvature_factor}")
            print(f"slope_factor: {slope_factor}")
            print(f"mu in dpdx: {mu}")
            print(f"D_h in dpdx: {D_h}")
            print(f"dA_dx: {dA_dx}")
            print(f"dp_stagnation_dx {dp_stagnation_dx}")
            print(f"dp_static_dx {dp_static_dx}\n")

        return dp_static_dx, dp_stagnation_dx
    
    def interface_temperatures(self, x, T_hw, T_cw):
        """
        Returns a list [T0, T1, ..., Tn] of wall-stack temperatures,
        where T0=T_hw, Tn=T_cw, and in between are each wall interface.
        """
        # same area-per-length as above
        dA_dx_hot = self.thrust_chamber.cooling_circuit_group.circuits[self.circuit_index].dA_dx_thermal_exhaust(x)

        # heat transfer per length
        qdx = self.dQ_cond_dx(x, T_hw, T_cw)
        T_rep = 0.5 * (T_hw + T_cw)
        Ts = [T_hw]
        for wall in self.thrust_chamber.wall_group.walls:
            Rj = wall.thickness(x) / (wall.material.get_k(T_rep) * dA_dx_hot)
            # drop in temperature across this layer:
            T_next = Ts[-1] - qdx * Rj
            Ts.append(T_next)
        # Ts[-1] should equal T_cw (within numerical roundoff)
        return Ts

    

def solve_heat_exchanger_euler(thrust_chamber, boundary_conditions, n_nodes, circuit_index, output, log_residuals=True):
    """
    Solve the 1D steady-state heat exchanger from x=0 to x=x_domain[-1].

    Arguments:
      engine       : your engine object containing geometry & property methods
      n_nodes      : number of axial nodes along the thrust chamber

    Returns:
      A dictionary with keys:
         "x"          -> 1D array of axial positions
         "T"          -> 2D array of temperatures, shape (n_nodes, 3):
                         columns = [T_coolant, T_wall_cold_side, T_wall_hot_side]
         "T_coolant"  -> 1D array of coolant temperatures
         "p_coolant"  -> 1D array of coolant pressures
         "dQ_dA"      -> 1D array of local heat fluxes (W/m^2)
         "velocity"   -> 1D array of coolant velocities (m/s)
    """

    # 1) Build the axial grid in [x_min, x_max].
    residual_log = [] if log_residuals else None
    iter_counter = np.zeros(n_nodes, dtype=int)
    
    circuit = thrust_chamber.cooling_circuit_group.circuits[circuit_index]
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
    physics_helper = HeatExchangerPhysics(thrust_chamber, circuit_index)

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
        #print(x_i)
        # Solve for T_hw, T_cw via algebraic system:
        def residuals(T_vars, cell=i):
            T_hw_trial, T_cw_trial = T_vars
            Q_hot_val  = physics_helper.dQ_hot_dx(x_i, T_hw_trial)*dx
            Q_cond_val = physics_helper.dQ_cond_dx(x_i, T_hw_trial, T_cw_trial)*dx
            Q_cold_val = physics_helper.dQ_cold_dx(x_i, T_cw_trial, T_cool_arr[i])*dx
            R1 = Q_hot_val - Q_cond_val
            R2 = Q_cond_val - Q_cold_val

            if residual_log is not None:
                k = iter_counter[cell]          # ← current Newton iteration for *this* cell
                residual_log.append((cell, k, R1, R2))
                iter_counter[cell] += 1
            
            return (R1, R2)

        sol = fsolve(residuals, x0=[T_hw_guess, T_cw_guess])
        
        T_hw_sol, T_cw_sol = sol
        #all_T_interfaces = physics_helper.interface_temperatures(x_i, T_hw_sol, T_cw_sol)

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
            #print(f"dp: {dp_static}")

            #print(dp)

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
        A_hot = thrust_chamber.cooling_circuit_group.circuits[circuit_index].dA_dx_thermal_exhaust(x_i) * dx
        dQ_dA_arr[i] = Q_hot_val / A_hot if A_hot != 0 else 0.0

    # 2. Compute the coolant velocity at each node
    velocity_arr = np.zeros(n_nodes)
    n_chan = thrust_chamber.cooling_circuit_group.circuits[circuit_index].placement.n_channel_positions*thrust_chamber.cooling_circuit_group.circuits[circuit_index].placement.n_channels_per_leaf
    for i, x_i in enumerate(x_domain):
        A_channel = thrust_chamber.cooling_circuit_group.circuits[circuit_index].A_coolant(x_i)
        mdot_c_single = thrust_chamber.combustion_transport.mdot_fu / n_chan
        rho_cool = thrust_chamber.cooling_circuit_group.circuits[circuit_index].coolant_transport.get_rho(T_cool_arr[i], p_static_arr[i])
        # Use the physics helper's u_coolant function 
        u_cool = physics.u_coolant(rho_cool, mdot_c_single, A_channel)
        velocity_arr[i] = u_cool
    
    T_stagnation_arr = np.zeros_like(T_cool_arr)
    for i, x_i in enumerate(x_domain):
        # Calculate cp for the current node using static temperature and pressure
        cp_cool = thrust_chamber.cooling_circuit_group.circuits[circuit_index].coolant_transport.get_cp(T_cool_arr[i], p_static_arr[i])
        # Compute stagnation temperature by adding the kinetic energy term
        T_stagnation_arr[i] = T_cool_arr[i] + (velocity_arr[i]**2) / (2.0 * cp_cool)

    n_walls = len(thrust_chamber.wall_group.walls)
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
        rho_i = thrust_chamber.cooling_circuit_group \
                   .circuits[circuit_index] \
                   .coolant_transport.get_rho(T_cool_arr[i], p_stagnation_arr[i])
        q_dyn = 0.5 * rho_i * velocity_arr[i]**2
        p_static_corrected[i] = p_stagnation_arr[i] - q_dyn

    # overwrite the old (bad) static-pressure array
    p_static_arr = p_static_corrected

    global_R, final_R = analyse_residuals(residual_log, n_nodes)
    #print(f"TP was acessed {physics_helper.counter} times")
    cooling_data = {
        "x"            : x_domain,
        "T"            : T_full,
        "T_static"     : T_cool_arr,
        "T_stagnation" : T_stagnation_arr,
        "p_static"     : p_static_arr,
        "p_stagnation" : p_stagnation_arr,
        "dQ_dA"        : dQ_dA_arr,
        "velocity"     : velocity_arr,
        "residuals"    : (global_R, final_R)
    }

    return cooling_data

def analyse_residuals(residual_log, n_cells, p=2):
    """
    Parameters
    ----------
    residual_log : list | None
        The list returned by `solve_channel`.  If None, nothing happens.
    n_cells : int
        Number of axial nodes in the simulation.
    p : int | float
        Order of the global norm: 2 for RMS, np.inf for L∞, etc.

    Returns
    -------
    history : (n_iter,) ndarray | None
        Global residual norm per iteration 0..k_max.
    final_per_cell : (n_cells,) ndarray | None
        Residual magnitude in each cell at the last local iteration.
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
    """
    Run the steady heating analysis.
    
    Parameters:
      engine  : Engine object with geometry, boundary conditions, and physics.
      n_nodes : Number of nodes (for Newton) or resolution for post-processing (for Radau).
      solver  : String, currently only "newton" available
    
    Returns:
      A dictionary with simulation results.
    """
    if solver.lower() == "newton":
        return solve_heat_exchanger_euler(thrust_chamber, boundary_conditions, n_nodes, circuit_index, output)
    else: 
        print("solver name not recognized")
    # possibility to implement other solvers here
