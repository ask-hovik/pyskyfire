import cantera as ct
from pyskyfire.common.fluids import Fluid
from pyskyfire.skycea.nozzle_solver import EquilNozzleInputs, EquilibriumNozzle
from typing import Optional
import numpy as np


# ================
# Helper Functions
# ================

def _mix_streams_to_chamber(mech, fuel_str, Tf, ox_str, To, p_c):
    gas_f = ct.Solution(mech); gas_f.TPX = Tf, p_c, fuel_str
    gas_o = ct.Solution(mech); gas_o.TPX = To, p_c, ox_str
    qf = ct.Quantity(gas_f, mass=1.0)   # relative basis; actual ratio is applied outside
    qo = ct.Quantity(gas_o, mass=1.0)
    qm = qf + qo                        # just to get a robust initial composition guess
    gas = ct.Solution(mech)
    gas.TPX = qm.T, p_c, qm.X
    gas.equilibrate('HP')               # chamber equilibrium at p_c
    return gas                           # has Tc, Xc, hc, sc

def _equil_tp_at_entropy(gas_ref, p_target, s0):
    """
    Return a new Cantera Solution equilibrated at (T,p_target) with s=s0.
    We enforce s=s0 by root-finding T at fixed p, using TP-equilibrium.
    """
    gas = gas_ref
    # Bracket T by a wide, safe range:
    T_lo, T_hi = 200.0, max(12000.0, gas_ref.T*2.5)

    def s_minus_s0(T):
        gas.TP = T, p_target
        gas.equilibrate('TP')  # eq at this TP
        return gas.entropy_mass - s0

    # Simple bisection
    f_lo, f_hi = s_minus_s0(T_lo), s_minus_s0(T_hi)
    if f_lo * f_hi > 0:
        # Fallback: nudge bounds if not bracketed (rare)
        for k in range(12):
            T_hi *= 1.5
            f_hi = s_minus_s0(T_hi)
            if f_lo * f_hi <= 0:
                break

    for _ in range(60):
        T_mid = 0.5*(T_lo + T_hi)
        f_mid = s_minus_s0(T_mid)
        if f_lo * f_mid <= 0:
            T_hi, f_hi = T_mid, f_mid
        else:
            T_lo, f_lo = T_mid, f_mid
        if abs(f_mid) < 1e-6 or (T_hi-T_lo) < 1e-4:
            break
    # gas currently at (T_mid, p_target) & TP-equilibrium
    return gas

def _mach_at_p(gas_chamber, p, h0, s0):
    g = _equil_tp_at_entropy(gas_chamber, p, s0)
    h = g.enthalpy_mass
    a = g.sound_speed
    u2 = max(0.0, 2.0*(h0 - h))  # J/kg
    u = np.sqrt(u2)
    M = u / a
    return M, g, u

def _find_throat(gas_chamber, p_c, h0, s0):
    # Find p* such that M=1
    p_hi = p_c
    p_lo = max(1.0, 1e-3 * p_c)  # 3–6 decades below chamber should be plenty
    M_hi, *_ = _mach_at_p(gas_chamber, p_hi, h0, s0)  # should be ~0
    M_lo, *_ = _mach_at_p(gas_chamber, p_lo, h0, s0)  # supersonic
    if M_lo < 1.0:
        # If even very low p isn’t supersonic, widen further
        p_lo = 1e-6 * p_c
        M_lo, *_ = _mach_at_p(gas_chamber, p_lo, h0, s0)

    for _ in range(60):
        p_mid = np.sqrt(p_lo * p_hi)  # log-bisection
        M_mid, g_mid, u_mid = _mach_at_p(gas_chamber, p_mid, h0, s0)
        if M_mid > 1.0:
            p_hi = p_mid
        else:
            p_lo = p_mid
        if abs(M_mid - 1.0) < 1e-6 or (p_hi/p_lo - 1) < 1e-6:
            return p_mid, g_mid  # g_mid is the throat state (M≈1)
    return p_mid, g_mid

def _area_ratio(g_star, u_star, g_state, u_state):
    # A/A* = (rho* a*) / (rho u)  ; with M*=1 we can use u*=a*
    return (g_star.density * u_star) / (g_state.density * u_state)

def _find_exit_for_eps(gas_chamber, p_c, h0, s0, eps, p_star, g_star, u_star):
    # Find p_e where A/A* = eps
    # Pressure must be < p*, so bracket between (very low) and p*.
    p_hi = p_star * 0.999
    p_lo = max(1.0, 1e-8 * p_c)

    def area_ratio_at_p(p):
        M, g, u = _mach_at_p(gas_chamber, p, h0, s0)
        return _area_ratio(g_star, u_star, g, u), g, u

    Ahi, _, _ = area_ratio_at_p(p_hi)
    Alo, _, _ = area_ratio_at_p(p_lo)
    # Alo should be huge (deeply expanded)
    for _ in range(60):
        p_mid = np.sqrt(p_lo * p_hi)
        Amid, g_mid, u_mid = area_ratio_at_p(p_mid)
        if Amid > eps:
            # not expanded enough → lower p
            p_lo = p_mid
        else:
            p_hi = p_mid
        if abs(Amid - eps)/eps < 1e-6 or (p_hi/p_lo - 1) < 1e-6:
            return p_mid, g_mid, u_mid
    return p_mid, g_mid, u_mid

def _stoich_of(fuel: Fluid, oxidizer: Fluid, mech: str) -> float:
    gas = ct.Solution(mech)

    def elem_totals(fluid: Fluid):
        mole_fracs = fluid.as_mole_fractions()
        elems = defaultdict(float)
        for sp, x in zip(fluid.propellants, mole_fracs):
            k = gas.species_index(sp)
            for el in gas.element_names:
                n = gas.n_atoms(k, el)
                if n:
                    elems[el] += x * n
        return elems

    def mix_mw(fluid: Fluid):
        mole_fracs = fluid.as_mole_fractions()
        return sum(x * gas.molecular_weights[gas.species_index(sp)]
                   for sp, x in zip(fluid.propellants, mole_fracs))

    # Element totals (per mole mixture)
    E_f = elem_totals(fuel)
    E_o = elem_totals(oxidizer)

    # Oxygen needed for complete combustion
    C, H, S, O_f = E_f.get('C',0), E_f.get('H',0), E_f.get('S',0), E_f.get('O',0)
    O_needed = 2*C + 0.5*H + 2*S - O_f
    if O_needed < 0: O_needed = 0.0

    O_avail = E_o.get('O', 0.0)
    if O_avail <= 0:
        raise ValueError("Oxidizer provides no oxygen atoms.")

    # Moles oxidizer per mole fuel
    b = O_needed / O_avail

    # Stoichiometric O/F mass ratio
    MW_f = mix_mw(fuel)
    MW_o = mix_mw(oxidizer)
    return (b * MW_o) / MW_f

def _cantera_mix_string(fluid: Fluid, mech: str) -> str:
    gas = ct.Solution(mech)

    mole_fracs = fluid.as_mole_fractions()
    parts = []
    for sp, x in zip(fluid.propellants, mole_fracs):
        if sp not in gas.species_names:
            raise ValueError(f"Species {sp} not in mechanism {mech}")
        parts.append(f"{sp}:{x:.6f}")  # keep 6 decimals for Cantera safety

    return ", ".join(parts)

# ================
# Combustion Class
# ================

def mach_at_pressure(gas, s0, h0, p):
    # Set state at fixed S,P and equilibrate (shifting-equilibrium model):
    gas.SP = s0, p
    gas.equilibrate("SP")       # remove this line for frozen model
    h = gas.enthalpy_mass
    a = gas.sound_speed
    u2 = 2*(h0 - h)
    if u2 <= 0.0:
        return 0.0              # not accelerating (or numerical)
    u = u2**0.5
    return u / a

def copy(gas, mech):
    """Return a fresh Solution object with the same state and composition as `gas`."""
    mech = mech or gas.source  # mechanism file if you know it
    gnew = ct.Solution(mech)
    gnew.TPX = gas.T, gas.P, gas.X
    return gnew

class CombustionReaction:
    """Equilibrium flow solver"""

    def __init__(self, optimum: dict[str, float | str], chemrep_map: Optional[dict[str, str]] = None):
        """Add all values in the optimum dict to self"""

        for key, value in optimum.items():
            setattr(self, key, value)
        self.optimum = optimum
        self.chemrep_map = chemrep_map or {}

    @classmethod
    def from_F_eps_Lstar(cls, fu, ox, MR, p_c, F, eps, L_star, A_t, mech, T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15):
        """Calculate optimal values using thrust, exit pressure and L-star"""

        # Compute equivalence ratio from MR
        #MR_stoich = stoich_of(fu, ox, mech=mech)
        #phi = MR/MR_stoich
        g0 = 9.81

        fu_str = _cantera_mix_string(fluid=fu, mech=mech) # get a fuel string compatible with cantera
        ox_str = _cantera_mix_string(fluid=ox, mech=mech) # ditto for ox

        gas_f = ct.Solution(mech); gas_f.TPX = T_fu_in, p_c, fu_str
        gas_o = ct.Solution(mech); gas_o.TPX = T_ox_in, p_c, ox_str

        mt = 1 + MR
        Yf = 1 / mt # find fuel mass fraction
        Yo = MR / mt # find ox mass fraction

        qf = ct.Quantity(gas_f, mass=Yf) # Yf is fuel mass fraction
        qo = ct.Quantity(gas_o, mass=Yo)
        qm = qf + qo  # adiabatic mixing: correct h, Y/X

        gas_c = ct.Solution(mech) # change transport model?
        gas_c.TPX = qm.T, p_c, qm.X  # set pre-reaction mixed state
        gas_c.equilibrate('HP')    # adiabatic, constant pressure

        Ae = eps * A_t                       # exit area from expansion ratio

        inp = EquilNozzleInputs(
            mech=mech,
            chamber_solution=gas_c,         # <-- use the chamber you just computed
            A_throat=A_t,
            A_exit=Ae,
            p_amb=p_amb
        )
        
        input("We got here!")
        noz = EquilibriumNozzle(inp)
        print("=== Throat ===")
        for k, v in noz.throat.items():
            print(f"{k:>8}: {v:.6g}" if isinstance(v, float) else f"{k:>8}: {v}")

        Ae_At = eps
        exit_state = noz.solve_exit_at_area_ratio(Ae_At)
        print("\n=== Exit (A/A* = %.3f) ===" % Ae_At)
        for k, v in exit_state.items():
            if k in ("Ae", "mdot"):
                print(f"{k:>8}: {v:.6g}")
            elif isinstance(v, float):
                print(f"{k:>8}: {v:.6g}")
            else:
                print(f"{k:>8}: {v}")

        F = noz.thrust(exit_state)
        print(f"\nThrust estimate: {F:.6g} N")




"""        # ===================== GENERATED ===============================
        T_c = gas_c.T
        h0  = gas_c.enthalpy_mass     # stagnation enthalpy
        s0  = gas_c.entropy_mass      # stagnation entropy
        rho_c = gas_c.density

        # 1) First order of action, assuming isentropic flow, (h0 and s0, T0 etc unchanged from chamber)
        # Find the pressure and temperature at the throat. 

        # bracket and solve M(p*) = 1
        p_hi = p_c
        p_lo = 0.01 * p_c                # or lower; decrease until M>1
        while mach_at_pressure(gas=copy(gas_c, mech), s0=s0, h0=h0, p=p_lo) <= 1.0:
            p_lo *= 0.5

        # simple bisection
        for _ in range(60):
            pm = 0.5*(p_lo+p_hi)
            M = mach_at_pressure(gas=copy(gas_c, mech), s0=s0, h0=h0, p=pm)
            if M > 1.0:
                p_hi = pm
            else:
                p_lo = pm
        p_t = 0.5*(p_lo+p_hi)

        print(s0)
        input("Here we are")

        # get throat state
        gas_t = ct.Solution(mech)
        gas_t.HP = h0, p_t
        gas_t.equilibrate("SP")           # remove for frozen
        T_t = gas_t.T
        a_t = gas_t.sound_speed
        h_t = gas_t.enthalpy_mass
        u_t = (2*(h0 - h_t))**0.5 # equals a_t (within numerics)
        rho_t = gas_t.density
        G_t = rho_t * u_t"""


        """# 3) Expand to exit area ratio eps (isentropic + equilibrium)
        p_e, g_e, u_e = _find_exit_for_eps(gas_c, p_c, h0, s0, eps, p_star, g_star, u_star)
        T_e, rho_e = g_e.T, g_e.density

        # 4) Performance, areas, mass flow, Isp
        # Thrust coefficient at ambient
        CF_vac = (u_e / c_star) + (p_e / p_c) * eps
        CF_amb = (u_e / c_star) + ((p_e - p_amb) / p_c) * eps

        mdot    = F / (CF_amb * c_star)
        A_t     = c_star * mdot / p_c
        A_e     = eps * A_t
        Isp_vac = CF_vac * c_star / g0
        Isp_amb = CF_amb * c_star / g0

        # Residence time & chamber volume
        t_stay = L_star * A_t * rho_c / mdot
        V_c    = mdot * t_stay / rho_c

        # Split flows by MR
        mdot_f = mdot / (1.0 + MR)
        mdot_o = mdot - mdot_f

        optimum = dict(
            MR=MR, p_c=p_c, F=F, eps=eps, L_star=L_star,
            # Chamber
            T_c=T_c, rho_c=rho_c, h0=h0, s0=s0,
            # Throat
            p_star=p_star, T_star=g_star.T, rho_star=g_star.density, a_star=a_star, G_star=G_star, c_star=c_star,
            # Exit
            p_e=p_e, T_e=T_e, rho_e=rho_e, u_e=u_e,
            # Performance
            CF_vac=CF_vac, CF_amb=CF_amb, Isp_vac=Isp_vac, Isp_amb=Isp_amb,
            # Flow rates & geometry
            mdot=mdot, mdot_fu=mdot_f, mdot_ox=mdot_o,
            A_t=A_t, A_e=A_e, t_stay=t_stay, V_c=V_c,
            # Strings (optional)
            # fuel_str=fu_str, ox_str=ox_str
        )"""

        optimum = dict(T_c=gas_c.T, T_t=gas_t.T, p_t=p_t)
        
        return cls(optimum)



