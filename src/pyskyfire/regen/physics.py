import numpy as np
from scipy.optimize import fsolve

def h_gas_bartz_enthalpy_driven(k_gr, D_hyd, Cp_gr, mu_gr, mdot_g, A_chmb, T_g, T_gr): 
    """Compute the hot-gas-side heat-transfer coefficient (Bartz correlation, enthalpy-based).

    Parameters
    ----------
    k_gr : float
        Thermal conductivity at the reference enthalpy condition [W m⁻¹ K⁻¹].
    D_hyd : float
        Hydraulic diameter of the gas passage [m].
    Cp_gr : float
        Specific heat capacity at the reference enthalpy condition [J kg⁻¹ K⁻¹].
    mu_gr : float
        Dynamic viscosity at the reference condition [Pa s].
    mdot_g : float
        Total hot-gas mass flow rate [kg s⁻¹].
    A_chmb : float
        Flow cross-sectional area in the chamber [m²].
    T_g : float
        Free-stream gas temperature [K].
    T_gr : float
        Reference-enthalpy temperature [K].

    Returns
    -------
    float
        Hot-gas-side heat-transfer coefficient ``h_g`` [W m⁻² K⁻¹].

    Notes
    -----
    This form of the Bartz correlation uses the ratio ``T_g/T_gr`` to
    account for variable-property effects between free-stream and reference conditions.

    subscript g denotes free stream gas properties
    subscript r denotes reference enthalpy conditions which are averaged between the free stream and wall metal conditions
    """

    return 0.026*k_gr/D_hyd*(Cp_gr/(k_gr*mu_gr))**0.4*(mdot_g*D_hyd/A_chmb)**0.8*(T_g/T_gr)**0.8

def h_gas_bartz(D_t, mu_g, cp_g, Pr_g, p_c, c_star, A_t, A_x, sigma):
    """Compute the Bartz hot-gas correlation without curvature correction.

    Parameters
    ----------
    D_t : float
        Throat diameter [m].
    mu_g : float
        Gas viscosity [Pa s].
    cp_g : float
        Gas specific heat capacity [J kg⁻¹ K⁻¹].
    Pr_g : float
        Gas Prandtl number.
    p_c : float
        Chamber pressure [Pa].
    c_star : float
        Characteristic velocity [m s⁻¹].
    A_t : float
        Throat area [m²].
    A_x : float
        Local flow area [m²].
    sigma : float
        Boundary-layer property-variation correction (dimensionless).

    Returns
    -------
    float
        Heat-transfer coefficient [W m⁻² K⁻¹].
    """
    h_g = 0.026/D_t**0.2*(mu_g**0.2*cp_g/Pr_g**0.6)*(p_c/c_star)**0.8*(A_t/A_x)**0.9*sigma
    return h_g

def sigma(T_hw, T_c, gamma_g, M_g, omega): 
    """Dimensionless property-variation factor used in Bartz correlation.

    Parameters
    ----------
    T_hw : float
        Wall temperature [K].
    T_c : float
        Core-flow (free-stream) temperature [K].
    gamma_g : float
        Ratio of specific heats.
    M_g : float
        Local Mach number.
    omega : float
        Empirical property exponent (≈ 0.68 for diatomic gases).

    Returns
    -------
    float
        ``σ`` correction factor (dimensionless).
    """
    sig = 1/((0.5*T_hw/T_c*(1 + (gamma_g - 1)/2*M_g**2) + 0.5)**(0.8 - omega/5)*(1 + (gamma_g -1)/2*M_g**2)**(omega/5))
    return sig

def h_coolant_colburn(k_cf, D_c, Cp_cr, mu_cf, mdot_c, A_c, phi_curv=1): 
    """Colburn correlation for turbulent coolant-side heat transfer.

    Parameters
    ----------
    k_cf : float
        Thermal conductivity of coolant film [W m⁻¹ K⁻¹].
    D_c : float
        Coolant hydraulic diameter [m].
    Cp_cr : float
        Mean specific heat (reference enthalpy) [J kg⁻¹ K⁻¹].
    mu_cf : float
        Dynamic viscosity of coolant film [Pa s].
    mdot_c : float
        Coolant mass-flow rate [kg s⁻¹].
    A_c : float
        Coolant-flow area [m²].
    phi_curv : float, optional
        Curvature-correction factor (dimensionless). Default is 1.

    Returns
    -------
    float
        Coolant-side heat-transfer coefficient [W m⁻² K⁻¹].
    
    Notes
    -----
    Subscripts: 
        Cf denotes film coolant film condition
        c denotes bulk coolant conditions
        r denotes reference enthalpy condition which are averaged between the free stream and wall metal conditions"""
    
    return 0.023*k_cf/D_c*(Cp_cr/(k_cf*mu_cf))**0.4*(mdot_c*D_c/A_c)**0.8*phi_curv

def u_coolant(rho, mdot_c_single_channel, A_cool): 
    """Compute coolant velocity.

    Parameters
    ----------
    rho : float
        Coolant density [kg m⁻³].
    mdot_c_single_channel : float
        Mass-flow rate through one channel [kg s⁻¹].
    A_cool : float
        Coolant cross-sectional area [m²].

    Returns
    -------
    float
        Coolant velocity [m s⁻¹].
    """
    u = mdot_c_single_channel / (rho * A_cool)
    return u

def reynolds(rho, u, L, mu): 
    """Compute the Reynolds number.

    Parameters
    ----------
    rho : float
        Fluid density [kg m⁻³].
    u : float
        Mean velocity [m s⁻¹].
    L : float
        Characteristic length or hydraulic diameter [m].
    mu : float
        Dynamic viscosity [Pa s].

    Returns
    -------
    float
        Reynolds number (dimensionless).
    """
    return rho*u*L/mu

def phi_curv(Re_c, D_c, R_curv): 
    """Curvature correction factor for the coolant-side heat-transfer coefficient.

    Parameters
    ----------
    Re_c : float
        Coolant Reynolds number.
    D_c : float
        Hydraulic diameter [m].
    R_curv : float
        Radius of curvature of the coolant passage [m].
        Use ``np.inf`` for straight sections.

    Returns
    -------
    float
        Curvature factor ``φ`` (dimensionless).
    """
    if R_curv == float("inf"):
        return 1
    else:
        phi = (Re_c*(0.5*D_c/R_curv)**2)**0.05
        return phi

ReDh_laminar = 2300         # Maximum Reynolds number for laminar flow in a pipe
ReDh_turbulent = 3500  # TODO: move to common constants page

def f_darcy_laminar(ReDh, Dh, x):
    """Laminar Darcy friction factor.

    Parameters
    ----------
    ReDh : float
        Hydraulic Reynolds number.
    Dh : float
        Hydraulic diameter [m].
    x : float
        Axial coordinate (unused, kept for API consistency).

    Returns
    -------
    float
        Darcy friction factor (dimensionless).
    """
    return 64.0 / ReDh      # Reference [3]

def f_darcy_turbulent(ReDh, Dh, x, roughness):
    """Turbulent Darcy friction factor using Putukhov or Colebrook–White.

    Parameters
    ----------
    ReDh : float
        Hydraulic Reynolds number.
    Dh : float
        Hydraulic diameter [m].
    x : float
        Axial coordinate.
    roughness : callable or None
        Function returning surface roughness at ``x`` [m], or ``None`` for smooth wall.

    Returns
    -------
    float
        Darcy friction factor (dimensionless).
    """
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
    """Composite laminar–turbulent Darcy friction factor with smooth transition.

    Parameters
    ----------
    ReDh : float
        Hydraulic Reynolds number.
    Dh : float
        Hydraulic diameter [m].
    x : float
        Axial coordinate.
    roughness : callable or None
        Function returning surface roughness at ``x`` [m], or ``None`` for smooth wall.

    Returns
    -------
    float
        Darcy friction factor (dimensionless).
    """
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
    """Compute the adiabatic-wall (recovery) temperature.

    Parameters
    ----------
    gamma : float
        Ratio of specific heats.
    M_inf : float
        Freestream Mach number.
    T_inf : float
        Freestream static temperature [K].
    Pr : float
        Prandtl number.

    Returns
    -------
    float
        Adiabatic-wall temperature [K].

    Notes
    -----
    The recovery factor ``r`` is approximated as ``Pr^(1/3)``,
    valid for turbulent flow. For laminar conditions, a smaller
    recovery factor should be used (≈ Pr⁰˙⁵).
    
    TODO: Implement a more robust version that could also handle laminar flow, in case this function is suddenly
    used in a function where this makes sense."""
    
    r = Pr**(1/3) # Recovery factor. Differs between turbulent and laminar
    return T_inf*(1 + r*((gamma -1)/2)*M_inf**2)
