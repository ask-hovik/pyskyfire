import numpy as np
from scipy.optimize import fsolve

def h_gas_bartz_enthalpy_driven(k_gr, D_hyd, Cp_gr, mu_gr, mdot_g, A_chmb, T_g, T_gr): 
    """Bartz correlation for hot gas side heat transfer coefficient
    
    Args: 
        k_gr is the 
        
    Subscripts: 
        g denotes free stream gas properties
        r denotes reference enthalpy conditions which are averaged between the free stream and wall metal conditions"""
    
    return 0.026*k_gr/D_hyd*(Cp_gr/(k_gr*mu_gr))**0.4*(mdot_g*D_hyd/A_chmb)**0.8*(T_g/T_gr)**0.8

def h_gas_bartz(D_t, mu_g, cp_g, Pr_g, p_c, c_star, A_t, A_x, sigma):
    """Bartz correlation from wikipedia without the curvature correction"""
    h_g = 0.026/D_t**0.2*(mu_g**0.2*cp_g/Pr_g**0.6)*(p_c/c_star)**0.8*(A_t/A_x)**0.9*sigma
    return h_g

def sigma(T_hw, T_c, gamma_g, M_g, omega): 
    """Dimensionless parameter accounting for variation of gas properties across the boundary layer"""
    sig = 1/((0.5*T_hw/T_c*(1 + (gamma_g - 1)/2*M_g**2) + 0.5)**(0.8 - omega/5)*(1 + (gamma_g -1)/2*M_g**2)**(omega/5))
    return sig

def h_coolant_colburn(k_cf, D_c, Cp_cr, mu_cf, mdot_c, A_c, phi_curv=1): 
    """Colburn correlation for coolant side heat transfer coefficient
    
    Args: 
        k_Cf is the conductivity in the coolant film
        D_C is the coolant tube hydraulic diameter
        cp_cr is then maybe the averaged cp between the film condition and the wall metal condition
        mu_Cf is the viscosity at the film conditions
        mdot_c is the coolant mass flow rate
        A_c is the cross sectional area of coolant flow
    Subscripts: 
        Cf denotes film coolant film condition
        c denotes bulk coolant conditions
        r denotes reference enthalpy condition which are averaged between the free stream and wall metal conditions"""
    
    return 0.023*k_cf/D_c*(Cp_cr/(k_cf*mu_cf))**0.4*(mdot_c*D_c/A_c)**0.8*phi_curv

def u_coolant(rho, mdot_c_single_channel, A_cool): 
    """ coolant velocity as a function of x"""
    u = mdot_c_single_channel / (rho * A_cool)
    return u

def reynolds(rho, u, L, mu): 
    return rho*u*L/mu

def phi_curv(Re_c, D_c, R_curv): 
    """ Curvature effect on the cold side heat transfer coefficient """
    if R_curv == float("inf"):
        return 1
    else:
        phi = (Re_c*(0.5*D_c/R_curv)**2)**0.05
        return phi

ReDh_laminar = 2300         # Maximum Reynolds number for laminar flow in a pipe
ReDh_turbulent = 3500 

def f_darcy_laminar(ReDh, Dh, x):
    return 64.0 / ReDh      # Reference [3]

def f_darcy_turbulent(ReDh, Dh, x, roughness):
    
    if roughness == None:
        # Putukhov equation [1]
        return (0.79 * np.log(ReDh) - 1.64)**(-2)   

    else:
        # Colebrook-White by iteration
        roughness_x = roughness(x)
        def func_to_solve(f):
            return 1/(f**0.5) + 2 * np.log10( roughness_x / (3.71 * Dh) + 2.51 / (ReDh * f**0.5) )    # Reference [2]

        return fsolve(func = func_to_solve, x0 = 0.02)[0]

def f_darcy(ReDh, Dh, x, roughness):
    # Check for laminar flow
    if ReDh < ReDh_laminar:
        return f_darcy_laminar(ReDh = ReDh, Dh = Dh, x = x)     

    elif ReDh < ReDh_turbulent:
        f_lam = f_darcy_laminar(ReDh = ReDh, Dh = Dh, x = x)     
        f_turb = f_darcy_turbulent(ReDh = ReDh, Dh = Dh, x = x, roughness=roughness)     

        # "Blend" between the laminar and turbulent region
        return np.interp(ReDh, [ReDh_laminar, ReDh_turbulent], [f_lam, f_turb])

    # Turbulent flow
    else:
        return f_darcy_turbulent(ReDh = ReDh, Dh = Dh, x = x, roughness=roughness)     

def T_aw(gamma, M_inf, T_inf, Pr):
        """ Adiabatic wall temperature, sometimes called recovery temperature T_r in some sources
        TODO: Implement a more robust version that could also handle laminar flow, in case this function is suddenly
        used in a function where this makes sense."""
        
        r = Pr**(1/3) # Recovery factor. Differs between turbulent and laminar
        return T_inf*(1 + r*((gamma -1)/2)*M_inf**2)
    