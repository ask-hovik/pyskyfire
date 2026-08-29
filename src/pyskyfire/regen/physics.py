import numpy as np
from scipy.optimize import brentq, fsolve


def kader_temperature_plus(y_plus, prandtl, eta):
    """Return Kader's dimensionless turbulent temperature profile.

    Parameters
    ----------
    y_plus : float or array-like
        Wall-normal distance in viscous units.
    prandtl : float
        Molecular Prandtl number at the representative fluid state.
    eta : float or array-like
        Outer coordinate ``y / delta``. Values must lie between zero at the
        wall and one at the thermal-boundary-layer edge.

    Returns
    -------
    float or numpy.ndarray
        Dimensionless temperature difference
        ``abs(T - T_wall) / T_tau``.

    Notes
    -----
    This is the smooth-wall, fully turbulent interpolation of Kader (1981),
    including its outer-layer correction. It blends the conductive sublayer
    and logarithmic layer without assuming a molecular Prandtl number near
    unity. Fluid properties are nevertheless treated as locally constant.
    """
    y_plus_array = np.asarray(y_plus, dtype=float)
    eta_array = np.asarray(eta, dtype=float)
    scalar = y_plus_array.ndim == 0 and eta_array.ndim == 0

    if not np.isfinite(prandtl) or prandtl <= 0.0:
        raise ValueError("prandtl must be finite and positive")
    if np.any(~np.isfinite(y_plus_array)) or np.any(y_plus_array < 0.0):
        raise ValueError("y_plus must be finite and non-negative")
    if np.any(~np.isfinite(eta_array)) or np.any(
        (eta_array < 0.0) | (eta_array > 1.0)
    ):
        raise ValueError("eta must be finite and lie in [0, 1]")

    y_plus_array, eta_array = np.broadcast_arrays(y_plus_array, eta_array)
    with np.errstate(over="ignore", invalid="ignore"):
        gamma = (
            0.01 * (prandtl * y_plus_array) ** 4
            / (1.0 + 5.0 * prandtl**3 * y_plus_array)
        )
    beta = (
        (3.85 * prandtl ** (1.0 / 3.0) - 1.3) ** 2
        + 2.12 * np.log(prandtl)
    )
    outer_factor = (
        1.5 * (2.0 - eta_array)
        / (1.0 + 2.0 * (1.0 - eta_array) ** 2)
    )

    inverse_gamma = np.full_like(gamma, np.inf)
    np.divide(1.0, gamma, out=inverse_gamma, where=gamma > 0.0)
    theta_plus = (
        prandtl * y_plus_array * np.exp(-gamma)
        + (
            2.12 * np.log((1.0 + y_plus_array) * outer_factor)
            + beta
        )
        * np.exp(-inverse_gamma)
    )
    theta_plus = np.where(y_plus_array == 0.0, 0.0, theta_plus)
    if scalar:
        return float(theta_plus)
    return theta_plus


def friction_velocity(velocity, darcy_factor):
    """Return friction velocity from a bulk velocity and Darcy factor.

    The Darcy-Weisbach convention gives
    ``u_tau = abs(velocity) * sqrt(f_D / 8)``.
    """
    velocity = float(velocity)
    darcy_factor = float(darcy_factor)
    if not np.isfinite(velocity):
        raise ValueError("velocity must be finite")
    if not np.isfinite(darcy_factor) or darcy_factor <= 0.0:
        raise ValueError("darcy_factor must be finite and positive")
    return abs(velocity) * np.sqrt(darcy_factor / 8.0)


def kader_temperature_profile(
    T_wall,
    T_edge,
    h,
    rho,
    cp,
    mu,
    prandtl,
    u_tau,
    n_points=200,
):
    """Reconstruct one local turbulent thermal boundary layer with Kader.

    The supplied heat-transfer coefficient fixes the convective wall flux as
    ``h * (T_edge - T_wall)``. The Kader edge condition is inverted for the
    friction Reynolds number, so this is a local profile reconstruction, not
    an axial boundary-layer-growth calculation.

    Returns
    -------
    distance : numpy.ndarray
        Distance from the wall, from zero to ``delta`` [m].
    temperature : numpy.ndarray
        Temperature at each returned distance [K].
    delta : float
        Reconstructed thermal-boundary-layer thickness [m].
    """
    values = {
        "T_wall": T_wall,
        "T_edge": T_edge,
        "h": h,
        "rho": rho,
        "cp": cp,
        "mu": mu,
        "prandtl": prandtl,
        "u_tau": u_tau,
    }
    values = {name: float(value) for name, value in values.items()}
    if any(not np.isfinite(value) for value in values.values()):
        raise ValueError("Kader profile inputs must be finite")
    for name in ("h", "rho", "cp", "mu", "prandtl", "u_tau"):
        if values[name] <= 0.0:
            raise ValueError(f"{name} must be positive")
    if int(n_points) != n_points or n_points < 2:
        raise ValueError("n_points must be an integer of at least two")

    T_wall = values["T_wall"]
    T_edge = values["T_edge"]
    if T_edge == T_wall:
        distance = np.zeros(int(n_points), dtype=float)
        return distance, np.full_like(distance, T_wall), 0.0

    temperature_span = abs(T_edge - T_wall)
    heat_flux = values["h"] * temperature_span
    friction_temperature = heat_flux / (
        values["rho"] * values["cp"] * values["u_tau"]
    )
    theta_edge = temperature_span / friction_temperature

    def log_edge_residual(log_re_tau):
        return (
            kader_temperature_plus(
                np.exp(log_re_tau), values["prandtl"], 1.0
            )
            - theta_edge
        )

    lower_log = -50.0
    upper_log = np.log(100.0)
    while log_edge_residual(upper_log) < 0.0 and upper_log < 700.0:
        upper_log += np.log(10.0)
    if log_edge_residual(upper_log) < 0.0:
        raise ValueError("Kader edge condition could not be bracketed")

    re_tau = np.exp(brentq(log_edge_residual, lower_log, upper_log))
    delta = re_tau * values["mu"] / (
        values["rho"] * values["u_tau"]
    )
    eta = np.linspace(0.0, 1.0, int(n_points))
    theta_plus = kader_temperature_plus(
        eta * re_tau,
        values["prandtl"],
        eta,
    )
    direction = np.sign(T_edge - T_wall)
    temperature = T_wall + direction * friction_temperature * theta_plus
    temperature[0] = T_wall
    temperature[-1] = T_edge
    distance = eta * delta
    return distance, temperature, float(delta)

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

def _curvature_group(Re_c, length_scale, R_curv):
    """Return ``Re * (length_scale / |R|)^2`` or zero if straight."""
    values = (Re_c, length_scale, R_curv)
    if not all(np.isfinite(value) for value in values):
        return 0.0
    if Re_c <= 0.0 or length_scale <= 0.0 or R_curv == 0.0:
        return 0.0
    return float(Re_c * (length_scale / abs(R_curv)) ** 2)


def phi_curv(Re_c, D_c, R_curv):
    """RL10 curvature factor for coolant-side Colburn heat transfer.

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
        Curvature factor ``φ`` (dimensionless). Positive radii are concave
        and enhance heat transfer; negative radii are convex and reduce it.

    Notes
    -----
    The RL10 model uses
    ``[Re * (0.5 D_h / R)^2]**(+/- 0.05)``.  Its expression has no sensible
    zero-curvature limit if extrapolated below a group of one, so this
    implementation returns one there.  That makes the straight-channel limit
    finite and joins the correlation continuously at its neutral value.
    """
    group = _curvature_group(Re_c, 0.5 * D_c, R_curv)
    if group <= 1.0:
        return 1.0
    exponent = 0.05 if R_curv > 0.0 else -0.05
    return group**exponent


def phi_curv_friction(Re_c, D_c, R_curv):
    """RTE/Ito curvature multiplier for the Darcy friction factor.

    The correlation uses hydraulic radius ``r_h = D_h/4`` and applies only
    for ``Re * (r_h/R)^2 > 6``.  It is independent of bend direction.
    """
    group = _curvature_group(Re_c, 0.25 * D_c, R_curv)
    if group <= 6.0:
        return 1.0
    return group**0.05

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
