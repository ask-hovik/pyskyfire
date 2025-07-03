import math
from typing import Sequence, Optional
from rocketcea.cea_obj_w_units import CEA_Obj as CEA_metric
from rocketcea.cea_obj import CEA_Obj as CEA_imperial

import numpy as np
import cantera as ct
import re
import scipy.optimize as opt


# -----------------------------------------------------------------------------
#  Low‑level helpers – kept separate for re‑use/testing
# -----------------------------------------------------------------------------

def _mach_from_area(A: float, At: float, gamma: float, supersonic: bool) -> float:
    """Invert area-ratio → Mach, picking sub- or supersonic branch."""

    def f(M):
        return A / At - 1.0 / M * ((2.0 + (gamma - 1.0) * M ** 2) / (gamma + 1.0))** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))

    if supersonic:
        bracket = (1.0, 100.0)
        x0 = 2.0
    else:
        bracket = (1e-12, 1.0)
        x0 = 0.1

    sol = opt.root_scalar(f, bracket=bracket, x0=x0, method="bisect")
    if not sol.converged:
        raise RuntimeError(f"root‑solve M(A) failed (A/At={A/At:.3f})")
    return sol.root

# -----------------------------------------------------------------------------
#  Thermodynamic helpers copied verbatim from the validated snippet
# -----------------------------------------------------------------------------

_EOS_LOOKUP = {"CH4": ct.Methane, "O2": ct.Oxygen, "H2": ct.Hydrogen}

import cantera as ct

_EOS_LOOKUP = {
    'O2': ct.Oxygen,
    'CH4': ct.Methane,
    'H2': ct.Hydrogen,
}

def flame_T(P, T_fuel_in, T_ox_in, phi, fuel_species, oxidizer_species, mech="gri30.yaml", T_lo=150.0, T_hi=4500.0, tol=0.5, max_iter=60):
    gas_ref = ct.Solution(mech)
    gas_ref.set_equivalence_ratio(1.0, fuel_species, f"{oxidizer_species}:1")
    Y_ref = gas_ref.Y
    m_f_st = Y_ref[gas_ref.species_index(fuel_species)]
    AFR_st = (1.0 - m_f_st) / m_f_st
    m_fuel = 1.0
    m_ox = AFR_st / phi

    def _h_real(species: str, T: float, P: float) -> float:
        eos = _EOS_LOOKUP.get(species)
        if eos is None:
            raise ValueError(f"No real-fluid EOS for {species}")
        s = eos()
        s.TP = T, P
        return s.enthalpy_mass  # J/kg

    H_react = m_fuel * _h_real(fuel_species, T_fuel_in, P) + m_ox * _h_real(oxidizer_species, T_ox_in, P)
    h_target = H_react / (m_fuel + m_ox)
    gas = ct.Solution(mech)
    gas.set_equivalence_ratio(phi, f"{fuel_species}:1", f"{oxidizer_species}:1")
    X0 = gas.X.copy()

    def h_eq(T):
        gas.TPX = T, P, X0
        gas.equilibrate("TP")
        return gas.enthalpy_mass

    f_lo = h_eq(T_lo) - h_target
    f_hi = h_eq(T_hi) - h_target
    if f_lo * f_hi > 0:
        raise RuntimeError("flame_T: failed to bracket root")

    for _ in range(max_iter):
        T_mid = 0.5 * (T_lo + T_hi)
        f_mid = h_eq(T_mid) - h_target
        if abs(f_mid) < tol:
            return T_mid
        if f_lo * f_mid < 0:
            T_hi, f_hi = T_mid, f_mid
        else:
            T_lo, f_lo = T_mid, f_mid
    raise RuntimeError("flame_T: bisection did not converge")


def flame_T_new(P, T_fuel_in, T_ox_in, phi,
            fuel_species, oxidizer_species,
            mech="gri30_highT.yaml",
            T_lo=150.0, T_hi=4500.0,
            tol=0.05, max_iter=60) -> float:
    """
    Compute adiabatic flame temperature for cold reactants using HP-equilibrium.

    Parameters:
    P                : Chamber pressure [Pa]
    T_fuel_in        : Inlet fuel temperature [K]
    T_ox_in          : Inlet oxidizer temperature [K]
    phi              : Equivalence ratio
    fuel_species     : Fuel species name (e.g., 'CH4')
    oxidizer_species: Oxidizer species name (e.g., 'O2')
    mech             : Cantera mechanism file
    T_lo, T_hi       : Not used (kept for compatibility)
    tol, max_iter    : Not used (kept for compatibility)

    Returns:
    T_ad             : Adiabatic flame temperature [K]
    """
    # 1. Reference mixture to get mass fractions
    gas_ref = ct.Solution(mech, transport_model="Multi")
    gas_ref.set_equivalence_ratio(1.0, fuel_species, f"{oxidizer_species}:1")
    Y_ref = gas_ref.Y
    Y_fuel = Y_ref[gas_ref.species_index(fuel_species)]
    Y_ox   = Y_ref[gas_ref.species_index(oxidizer_species)]

    # 2. Real-fluid enthalpy for each stream
    def _h_real(species: str, T: float, P: float) -> float:
        eos = _EOS_LOOKUP.get(species)
        if eos is None:
            raise ValueError(f"No real-fluid EOS for {species}")
        comp = eos()
        comp.TP = T, P
        return comp.enthalpy_mass

    # 3. Compute target specific enthalpy of reactants (mass-weighted)
    h_f = _h_real(fuel_species, T_fuel_in, P)
    h_o = _h_real(oxidizer_species, T_ox_in, P)
    h_target = Y_fuel * h_f + Y_ox * h_o

    # 4. Fresh mixture composition
    gas = ct.Solution(mech, transport_model="Multi")
    gas.set_equivalence_ratio(phi, f"{fuel_species}:1", f"{oxidizer_species}:1")
    X0 = gas.X.copy()

    # 5. Set to target enthalpy and pressure, then equilibrate
    gas.TPX = 300.0, P, X0  # initial guess
    gas.HP = h_target, P
    try:
        gas.equilibrate('HP')
    except Exception as e:
        raise RuntimeError(f"HP-equilibration failed: {e}")

    return gas.T

def cp_cv_eq(gas: ct.Solution, T: float, P: float, dT: float = 0.01):
    gas.TP = T - dT/2, P
    gas.equilibrate("TP")
    h_lo = gas.enthalpy_mass; u_lo = gas.int_energy_mass

    gas.TP = T + dT/2, P
    gas.equilibrate("TP")
    h_hi = gas.enthalpy_mass; u_hi = gas.int_energy_mass

    cp = (h_hi - h_lo) / dT
    cv = (u_hi - u_lo) / dT
    return cp, cv, cp/cv

# -----------------------------------------------------------------------------
#                              Main class
# -----------------------------------------------------------------------------

class CombustionTransport:
    """Quasi-1-D equilibrium nozzle flow solver."""

    def __init__(self, 
                 ox: str, 
                 fu: str, 
                 p_c: float, 
                 p_e: float,  
                 MR = None,
                 phi = None,
                 mech: str = "gri30_highT.yaml", 
                 T_fuel_in: float = 298.15, 
                 T_ox_in: float = 298.15, 
                 tol: float = 1e-4, 
                 max_iter: int = 25, 
                 mode="hybrid"):
        
        self.ox, self.fu = ox, fu
        self.mode = mode
        self.p_c, self.p_e = float(p_c), float(p_e)
        self.gas = ct.Solution(mech, transport_model="Multi")
        self.mech = mech
        self.n_species = self.gas.n_species

        self.T_fuel_in = T_fuel_in
        self.T_ox_in = T_ox_in

        self.tol = tol
        self.max_iter = max_iter
        
        if MR == None and phi != None:
            self.MR = self._find_MR_from_phi(phi)
            self.phi = phi
        elif phi == None and MR != None: 
            self.phi = self._find_phi_from_MR(MR)
            self.MR = MR
            #print(f"MR: {self.MR}")
            #print(f"phi: {self.phi}")
        else:
            raise ValueError("Either MR or phi must be supplied to CombustionTransport")
        
        self._solve_chamber(self.T_fuel_in, self.T_ox_in)
    
    def set_essentials(self, mdot, mdot_fu, mdot_ox):
        self.mdot = mdot
        self.mdot_fu = mdot_fu
        self.mdot_ox = mdot_ox


    def compute_transport(self, contour): 
        
        self.contour = contour
        self._solve_nozzle(self.tol, self.max_iter)

    def _find_MR_from_phi(self, phi):
        tmp = ct.Solution(self.mech, transport_model="Multi")
        tmp.set_equivalence_ratio(1.0, f"{self.fu}:1", f"{self.ox}:1")
        Y = tmp.Y; m_f = Y[tmp.species_index(self.fu)]
        MR_stoich = (1.0 - m_f) / m_f
        return MR_stoich/phi

    def _find_phi_from_MR(self, MR): 
        tmp = ct.Solution(self.mech, transport_model="Multi")
        tmp.set_equivalence_ratio(1.0, f"{self.fu}:1", f"{self.ox}:1")
        Y = tmp.Y; m_f = Y[tmp.species_index(self.fu)]
        MR_stoich = (1.0 - m_f) / m_f
        return MR_stoich/MR
    
    def calculate_c_star(self, gamma, temperature, molecular_weight):
        return (np.sqrt(ct.gas_constant * temperature / (molecular_weight * gamma)) * np.power(2 / (gamma + 1), -(gamma + 1) / (2*(gamma - 1))))

    def _solve_chamber(self, T_fuel_in, T_ox_in):
        self.T_c = flame_T(self.p_c, T_fuel_in, T_ox_in, self.phi, self.fu, self.ox, self.mech)
        self.gas.set_equivalence_ratio(self.phi, f"{self.fu}:1", f"{self.ox}:1")
        self.gas.TP = self.T_c, self.p_c
        self.gas.equilibrate("TP")
        

        self.cp_c, self.cv_c, self.gamma_c = cp_cv_eq(self.gas, self.T_c, self.p_c)
        self.h_c = self.gas.enthalpy_mass
        self.rho_c = self.gas.density
        self.molwt_c = self.gas.mean_molecular_weight
        self.c_star = self.calculate_c_star(self.gamma_c, self.T_c, self.molwt_c)

    def _solve_nozzle(self, tol, max_iter, n_nodes = 100):
        x_min = self.contour.xs[0]
        x_max = self.contour.xs[-1]
        # You had x_domain going from x_max -> x_min before; you can keep that if you like,
        # but for interpolation we typically store these in ascending order:
        self.x_domain = np.linspace(x_min, x_max, n_nodes)

        x_t = self.contour.x_t; 
        A_t = self.contour.A_t

        n = len(self.x_domain)
        gamma_old = np.full(n, self.gamma_c)

        # allocate
        self.M_map = np.empty(n) 
        self.gamma_map = np.empty(n)
        self.T_map = np.empty(n)
        self.p_map = np.empty(n)
        self.cp_map = np.empty(n)
        self.cv_map = np.empty(n)
        self.h_map = np.empty(n)
        self.k_map = np.empty(n)
        self.mu_map = np.empty(n)
        self.Pr_map = np.empty(n)
        self.rho_map = np.empty(n)
        self.Y_map = np.empty((n, self.n_species))

        if self.fu == "H2": 
            cea_fuel = "LH2"
        elif self.fu == "CH4":
            cea_fuel = "CH4"

        fuy = "CH4"
        oxy = "LOX"
        K = CEA_imperial(oxName=self.ox, fuelName=self.fu)
        print(f"WARNING: Due to short timeframe the cea fuel is hardcoded to {fuy} and {oxy}")

        for _ in range(max_iter):
            for i, x in enumerate(self.x_domain):

                A_x = self.contour.A(x)
                supersonic = x >= x_t

                """if self.mode == "hybrid":
                    if x < 0: 
                        cp_cea, mu_cea_mP, k_cea_WcmK, Pr_cea = C.get_Chamber_Transport(Pc=self.p_c, MR=self.MR)
                        k_cea = k_cea_WcmK*100
                        mu_cea = mu_cea_mP*1e-4
                        self.k_map[i] = k_cea
                        self.mu_map[i] = mu_cea
                        self.cp_map[i] = cp_cea
                        self.Pr_map[i] = cp_cea*mu_cea/k_cea

                    elif x >= 0: 
                        cp_cea, mu_cea_mP, k_cea_WcmK, Pr_cea = C.get_Exit_Transport(Pc=self.p_c, MR=self.MR, eps=A_x/A_t, frozen=0, frozenAtThroat=0)
                        k_cea = k_cea_WcmK*100 
                        mu_cea = mu_cea_mP*1e-4 

                        # CEA values
                        self.k_map[i] = k_cea 
                        self.mu_map[i] = mu_cea 
                        self.cp_map[i] = cp_cea 
                        self.Pr_map[i] = cp_cea*mu_cea/k_cea

                    # Cantera values
                    self.M_map[i] = M_i; 
                    self.gamma_map[i] = gamma_i; 
                    self.T_map[i] = T_i; 
                    self.p_map[i] = p_i
                    self.cv_map[i] = cv_i; 
                    self.h_map[i] = self.gas.enthalpy_mass
                    self.rho_map[i] = self.gas.density
                    self.Y_map[i, :] = self.gas.Y"""
                        

                if self.mode == "cea":
                    if x < 0: 
                        #result = get_full(C=K, Pc=self.p_c, MR=self.MR, M=M_cea, eps=A_x/A_t)
                        #input()

                        full_output = K.get_full_cea_output(Pc=self.p_c/1e5, MR=self.MR, subar=None, PcOvPe=None, frozen=0, frozenAtThroat=0, show_transport=1, pc_units='bar', output="metric", show_mass_frac=True, fac_CR=A_x/A_t)
                        #res = extract_full_cea_properties(full_output)
                        #res2 = print_full_output_lines(full_output)
                        #res3 = parse_and_display_properties(full_output)
                        props = parse_cea_properties(full_output)
                        #for k, v in props.items():
                        #    print(f"{k:12s}: {v}")

                        self.k_map[i] = props["conductivity"][1]
                        self.mu_map[i] = props["viscosity"][1]
                        self.cp_map[i] = props["cp"][1]
                        self.gamma_map[i] = props["gamma"][1]
                        self.M_map[i] = props["mach"][1]
                        self.T_map[i] = props["temperature"][1]
                        self.p_map[i] = props["pressure"][1]
                        self.cv_map[i] = props["cv"][1]
                        self.h_map[i] = props["enthalpy"][1]
                        self.rho_map[i] = props["density"][1]
                        self.Pr_map[i] = props["prandtl"][1]
                        

                        """cp_cea, mu_cea_mP, k_cea_WcmK, Pr_cea = C.get_Chamber_Transport(Pc=self.p_c, MR=self.MR)
                        k_cea = k_cea_WcmK*100
                        mu_cea = mu_cea_mP*1e-4

                        molwt_cea, gamma_cea = C.get_Chamber_MolWt_gamma(Pc=self.p_c, MR=self.MR)
                        
                        T_cea, _, _ = C.get_Temperatures(Pc=self.p_c, MR=self.MR, eps=A_x/A_t)
                        h_cea, _, _ = C.get_Enthalpies(Pc=self.p_c, MR=self.MR, eps=A_x/A_t)
                        rho_cea = C.get_Chamber_Density(Pc=self.p_c, MR=self.MR, eps=A_x/A_t)"""
                        
                        """self.k_map[i] = k_cea
                        self.mu_map[i] = mu_cea
                        self.cp_map[i] = cp_cea
                        self.gamma_map[i] = gamma_cea
                        self.M_map[i] = M_cea 
                        self.T_map[i] = T_cea 
                        self.p_map[i] = self.p_c
                        self.cv_map[i] = cp_cea/gamma_cea
                        self.h_map[i] = h_cea
                        self.rho_map[i] = rho_cea
                        self.Pr_map[i] = Pr_cea"""

                    elif x >= 0: 
                        full_output = K.get_full_cea_output(Pc=self.p_c/1e5, MR=self.MR, eps=A_x/A_t, subar=None, PcOvPe=None, frozen=0, frozenAtThroat=0, show_transport=1, pc_units='bar', output="metric", show_mass_frac=True, fac_CR=A_x/A_t)
                        #res = extract_full_cea_properties(full_output)
                        #res2 = print_full_output_lines(full_output)
                        #res3 = parse_and_display_properties(full_output)
                        props = parse_cea_properties(full_output)
                        #for k, v in props.items():
                        #    print(f"{k:12s}: {v}")

                        

                        self.k_map[i] = props["conductivity"][3]
                        self.mu_map[i] = props["viscosity"][3]
                        self.cp_map[i] = props["cp"][3]
                        self.gamma_map[i] = props["gamma"][3]
                        self.M_map[i] = props["mach"][3]
                        self.T_map[i] = props["temperature"][3]
                        self.p_map[i] = props["pressure"][3]
                        self.cv_map[i] = props["cv"][3]
                        self.h_map[i] = props["enthalpy"][3]
                        self.rho_map[i] = props["density"][3]
                        self.Pr_map[i] = props["prandtl"][3]
                        
                        """cp_cea, mu_cea_mP, k_cea_WcmK, Pr_cea = C.get_Exit_Transport(Pc=self.p_c, MR=self.MR, eps=A_x/A_t, frozen=0, frozenAtThroat=0)
                        molwt_cea, gamma_cea = C.get_exit_MolWt_gamma(Pc=self.p_c, MR=self.MR, eps=A_x/A_t, frozen=0, frozenAtThroat=0)
                        M_cea = K.get_MachNumber(Pc=self.p_c*0.000145038, MR=self.MR, eps=A_x/A_t)
                        _, _, T_cea = C.get_Temperatures(Pc=self.p_c, MR=self.MR, eps=A_x/A_t)
                        pcovpe = C.get_PcOvPe(Pc=self.p_c, MR=self.MR, eps=A_x/A_t)
                        p_cea = self.p_c/pcovpe
                        _, _, h_cea = C.get_Enthalpies(Pc=self.p_c, MR=self.MR, eps=A_x/A_t)
                        _, _, rho_cea = C.get_Densities(Pc=self.p_c, MR=self.MR, eps=A_x/A_t)
                        

                        k_cea = k_cea_WcmK*100
                        mu_cea = mu_cea_mP*1e-4
                        self.k_map[i] = k_cea
                        self.mu_map[i] = mu_cea
                        self.cp_map[i] = cp_cea
                        self.gamma_map[i] = gamma_cea
                        self.M_map[i] = M_cea
                        self.T_map[i] = T_cea
                        self.p_map[i] = p_cea
                        self.cv_map[i] = cp_cea/gamma_cea
                        self.h_map[i] = h_cea
                        self.rho_map[i] = rho_cea
                        self.Pr_map[i] = Pr_cea"""

                    self.Y_map[i, :] = self.gas.Y
                    #self.h_map[i] = self.gas.enthalpy_mass

                elif self.mode == "hybrid":
                    if x < 0: 
                        full_output = K.get_full_cea_output(Pc=self.p_c/1e5, MR=self.MR, subar=None, PcOvPe=None, frozen=0, frozenAtThroat=0, show_transport=1, pc_units='bar', output="metric", show_mass_frac=True, fac_CR=A_x/A_t)
                        props = parse_cea_properties(full_output)

                        # Decide which of the four CEA columns to use
                        mf_dict = props["mass_fractions"]
                        T_dict = props["temperature"]
                        p_dict = props["pressure"]
                        M_dict = props["mach"]
                        col = 1
                        M_i = M_dict[col]
                        T_i = T_dict[col]
                        p_i = p_dict[col]

                        # --- 2b)  Build a Cantera Y-vector -------------------
                        Y = np.zeros(self.n_species)
                        for sp, vals in mf_dict.items():
                            if sp in self.gas.species_names:
                                idx = self.gas.species_index(sp)
                                Y[idx] = vals[col]
                        # catch any rounding loss
                        if Y.sum() > 0:
                            Y /= Y.sum()

                        # --- 2c)  Freeze that composition and ask Cantera ----
                        self.gas.TPY = T_i, p_i, Y        # no equilibrate!
                        cp_i, cv_i, gamma_i = cp_cv_eq(self.gas, T=T_i, P=p_i)
                        #cp_i  = self.gas.cp_mass
                        #cv_i  = self.gas.cv_mass
                        #gamma_i = cp_i / cv_i
                        mu_i  = self.gas.viscosity
                        #k_i   = self.gas.thermal_conductivity
                        k_i = props["conductivity"][col]
                        Pr_i  = cp_i * mu_i / k_i
                        h_i   = self.gas.enthalpy_mass
                        rho_i = self.gas.density

                        # --- 2d)  Store --------------------------------------
                        self.M_map[i]     = M_i
                        self.gamma_map[i] = gamma_i
                        self.T_map[i]     = T_i
                        self.p_map[i]     = p_i
                        self.cp_map[i]    = cp_i
                        self.cv_map[i]    = cv_i
                        self.h_map[i]     = h_i
                        self.k_map[i]     = k_i
                        self.mu_map[i]    = mu_i
                        self.Pr_map[i]    = Pr_i
                        self.rho_map[i]   = rho_i
                        self.Y_map[i, :]  = Y

                    elif x >= 0: 
                        full_output = K.get_full_cea_output(Pc=self.p_c/1e5, MR=self.MR, eps=A_x/A_t, subar=None, PcOvPe=None, frozen=0, frozenAtThroat=0, show_transport=1, pc_units='bar', output="metric", show_mass_frac=True, fac_CR=A_x/A_t)
                        props = parse_cea_properties(full_output)

                        
                        mf_dict = props["mass_fractions"]
                        T_dict = props["temperature"]
                        p_dict = props["pressure"]
                        M_dict = props["mach"]
                        col = 3
                        M_i = M_dict[col]
                        T_i = T_dict[col]
                        p_i = p_dict[col]

                        # --- 2b)  Build a Cantera Y-vector -------------------
                        Y = np.zeros(self.n_species)
                        for sp, vals in mf_dict.items():
                            if sp in self.gas.species_names:
                                idx = self.gas.species_index(sp)
                                Y[idx] = vals[col]
                        # catch any rounding loss
                        if Y.sum() > 0:
                            Y /= Y.sum()

                        # --- 2c)  Freeze that composition and ask Cantera ----
                        self.gas.TPY = T_i, p_i, Y        # no equilibrate!
                        cp_i, cv_i, gamma_i = cp_cv_eq(self.gas, T=T_i, P=p_i)

                        #cp_i  = self.gas.cp_mass
                        #cv_i  = self.gas.cv_mass
                        #gamma_i = cp_i / cv_i
                        mu_i  = self.gas.viscosity
                        
                        #k_i   = self.gas.thermal_conductivity
                        k_i = props["conductivity"][col]
                        Pr_i  = cp_i * mu_i / k_i
                        h_i   = self.gas.enthalpy_mass
                        rho_i = self.gas.density

                        # --- 2d)  Store --------------------------------------
                        self.M_map[i]     = M_i
                        self.gamma_map[i] = gamma_i
                        self.T_map[i]     = T_i
                        self.p_map[i]     = p_i
                        self.cp_map[i]    = cp_i
                        self.cv_map[i]    = cv_i
                        self.h_map[i]     = h_i
                        self.k_map[i]     = k_i
                        self.mu_map[i]    = mu_i
                        self.Pr_map[i]    = Pr_i
                        self.rho_map[i]   = rho_i
                        self.Y_map[i, :]  = Y

                elif self.mode == "dual":
                    if x < 0: 
                        full_output = K.get_full_cea_output(Pc=self.p_c/1e5, MR=self.MR, subar=None, PcOvPe=None, frozen=0, frozenAtThroat=0, show_transport=1, pc_units='bar', output="metric", show_mass_frac=True, fac_CR=A_x/A_t)
                        col = 1
                    elif x >= 0: 
                        full_output = K.get_full_cea_output(Pc=self.p_c/1e5, MR=self.MR, eps=A_x/A_t, subar=None, PcOvPe=None, frozen=0, frozenAtThroat=0, show_transport=1, pc_units='bar', output="metric", show_mass_frac=True, fac_CR=A_x/A_t)
                        col = 3
                    
                    props = parse_cea_properties(full_output)

                    # Decide which of the four CEA columns to use
                    mf_dict = props["mass_fractions"]
                    M_dict = props["mach"]
                    T_dict = props["temperature"]
                    p_dict = props["pressure"]

                    
                    M_i = M_dict[col]
                    T_i = T_dict[col]
                    p_i = p_dict[col]

                    # --- 2b)  Build a Cantera Y-vector -------------------
                    Y = np.zeros(self.n_species)
                    for sp, vals in mf_dict.items():
                        if sp in self.gas.species_names:
                            idx = self.gas.species_index(sp)
                            Y[idx] = vals[col]
                    # catch any rounding loss
                    if Y.sum() > 0:
                        Y /= Y.sum()

                    # --- 2c)  Freeze that composition and ask Cantera ----
                    self.gas.TPY = T_i, p_i, Y        # no equilibrate!
                    
                    cp_i, cv_i, gamma_i = cp_cv_eq(self.gas, T=T_i, P=p_i)
                    #cp_i  = self.gas.cp_mass
                    #cv_i  = self.gas.cv_mass
                    #gamma_i = cp_i / cv_i
                    mu_i  = self.gas.viscosity
                    k_i   = self.gas.thermal_conductivity
                    #k_i = props["conductivity"][col]
                    Pr_i  = cp_i * mu_i / k_i
                    h_i   = self.gas.enthalpy_mass
                    rho_i = self.gas.density

                    # --- 2d)  Store --------------------------------------
                    self.M_map[i]     = M_i
                    self.gamma_map[i] = gamma_i
                    self.T_map[i]     = T_i
                    self.p_map[i]     = p_i
                    self.cp_map[i]    = cp_i
                    self.cv_map[i]    = cv_i
                    self.h_map[i]     = h_i
                    self.k_map[i]     = k_i
                    self.mu_map[i]    = mu_i
                    self.Pr_map[i]    = Pr_i
                    self.rho_map[i]   = rho_i
                    self.Y_map[i, :]  = Y

                elif self.mode == "cantera":

                    M_i = _mach_from_area(A_x, A_t, gamma_old[i], supersonic)
                    
                    T_i = self.T_c / (1.0 + 0.5 * (gamma_old[i] - 1.0) * M_i ** 2)
                    p_i = self.p_c * (T_i / self.T_c) ** (gamma_old[i] / (gamma_old[i] - 1.0))
                    #T_ad = flame_T(p_i, T_fuel_in=298.15, T_ox_in=298.15, phi=self.phi, fuel_species=self.fu, oxidizer_species=self.ox, mech=self.mech)
                    #p_i = self.p_c * (T_i / self.T_c) ** (gamma_old[i] / (gamma_old[i] - 1.0))
                    #T_i = T_ad / (1.0 + 0.5 * (gamma_old[i] - 1.0) * M_i ** 2)

                    self.gas.TP = T_i, p_i
                    self.gas.equilibrate("TP")

                    cp_i, cv_i, gamma_i = cp_cv_eq(self.gas, T_i, p_i)
                    # eventually want to have just this one
                    #self.gas.transport_model = "multicomponent"
                    # store
                    self.M_map[i] = M_i; 
                    self.gamma_map[i] = gamma_i; 
                    self.T_map[i] = T_i; 
                    self.p_map[i] = p_i
                    self.cp_map[i] = cp_i; 
                    self.cv_map[i] = cv_i; 
                    self.h_map[i] = self.gas.enthalpy_mass
                    self.k_map[i] = self.gas.thermal_conductivity
                    self.mu_map[i] = self.gas.viscosity
                    self.Pr_map[i] = cp_i*self.gas.viscosity/self.gas.thermal_conductivity
                    self.rho_map[i] = self.gas.density
                    self.Y_map[i, :] = self.gas.Y
                
                else: 
                    raise ValueError("Unrecognized mode input in CombustionTransport")
                

                

                # TODO: A few notes on the above. Unfortunately, some of the cantera transport properties differs ignificantly
                # from the cea ones (mainly k). I believe it is a case of cantera using theoretical correlations rather than experimental ones. 
                # I believe it is possible to drastically increase the accuracy of the cantera estimates by supplying a yaml file
                # with the same data that cea uses, which is based on experiments. I think we could also make a much faster 
                # chamber temperature estimate by extending or even just extrapolating the polynomials down to the liquid vapor line
                # for the different species. Although there may be good reasons not to do this. For now this cantera-cea merge is 
                # _much_ faster than the previous implementation. The only property that is currently fully using the cantera perks 
                # is h. In the future, if the cantera methods match cea more closely, a more precise hot gas side heat transfer 
                # coefficient can be estimated. There is currently a lot of other things in the program that needs attention, 
                # therefore I will keep it like this for now. If we can fully switch to cantera in the future, the simulations may
                # become blazingly fast.  
                

            if np.max(np.abs(self.gamma_map - gamma_old)) < tol:
                return
            
            gamma_old = self.gamma_map
        raise RuntimeError("CombustionTransport: gamma-M iteration did not converge")
    

    # ================
    # internal helpers
    # ================

    def _interp_baseline(self, arr: np.ndarray, x: float) -> float:
        """1-D linear interpolation for scalar baseline arrays."""
        return float(np.interp(x, self.x_domain, arr))

    def _interp_Y(self, x: float) -> np.ndarray:
        """Piece-wise linear interpolation for the Y-vector baseline."""
        if x <= self.x_domain[0]:
            return self.Y_map[0]
        if x >= self.x_domain[-1]:
            return self.Y_map[-1]
        idx = np.searchsorted(self.x_domain, x) - 1
        x0, x1 = self.x_domain[idx], self.x_domain[idx + 1]
        w = (x - x0) / (x1 - x0)
        return (1.0 - w) * self.Y_map[idx] + w * self.Y_map[idx + 1]

    # ------------------------------------------------------------------
    # NEW: state assembly & evaluation helpers
    # ------------------------------------------------------------------
    def _coalesce(
        self,
        x: float,
        T: Optional[float] = None,
        H: Optional[float] = None,
        p: Optional[float] = None,
        Y: Optional[np.ndarray] = None,
    ):
        """
        Gather a complete state description.

        Exactly **one** of T or H may be supplied (raise if both).
        Any quantity left as None is filled from the baseline maps.
        Returns (T, H, p, Y) where either T xor H is not-None.
        """
        if (T is not None) and (H is not None):
            raise ValueError("Pass either T or H, not both.")

        # choose baseline T when neither T nor H given
        if T is None and H is None:
            T = self._interp_baseline(self.T_map, x)

        if p is None:
            p = self._interp_baseline(self.p_map, x)

        if Y is None:
            Y = self._interp_Y(x)

        return T, H, p, Y

    def _eval_state(
        self,
        T: Optional[float],
        H: Optional[float],
        p: float,
        Y: np.ndarray,
    ) -> ct.Solution:
        """
        Push the requested state into `self.gas` (TPY or HPY) and return it.
        """
        g = self.gas
        if H is not None:
            g.HPY = H, p, Y
        else:
            g.TPY = T, p, Y
        return g

    # ------------------------------------------------------------------
    # baseline getters (unchanged behaviour)
    # ------------------------------------------------------------------
    def get_M(self, x: float) -> float:
        return self._interp_baseline(self.M_map, x)

    def get_T(self, x: float,
              *, H: Optional[float] = None,
              p: Optional[float] = None,
              Y: Optional[np.ndarray] = None) -> float:
        """
        Temperature baseline, or value consistent with a supplied enthalpy.

        • original call      → get_T(x)
        • new functionality  → get_T(x, H=… [,p=… ,Y=…])
        """
        if H is None and p is None and Y is None:
            return self._interp_baseline(self.T_map, x)
        T, H, p, Y = self._coalesce(x, None, H, p, Y)
        return self._eval_state(T, H, p, Y).T

    def get_p(self, x: float) -> float:
        """Baseline pressure – unchanged (no overrides allowed)."""
        return self._interp_baseline(self.p_map, x)

    def get_Y(self, x: float) -> np.ndarray:
        """Baseline composition – unchanged (no overrides allowed)."""
        return self._interp_Y(x)

    # ------------------------------------------------------------------
    # scalar property helpers (DRY)
    # ------------------------------------------------------------------
    def _scalar_getter(
        self,
        x: float,
        *,
        baseline_array: np.ndarray,
        field: str,
        T: Optional[float] = None,
        H: Optional[float] = None,
        p: Optional[float] = None,
        Y: Optional[np.ndarray] = None,
    ) -> float:
        if T is None and H is None and p is None and Y is None:
            return self._interp_baseline(baseline_array, x)
        T, H, p, Y = self._coalesce(x, T, H, p, Y)
        return getattr(self._eval_state(T, H, p, Y), field)

    # ------------------------------------------------------------------
    # transport / thermodynamic getters
    # ------------------------------------------------------------------
    def get_cp(self, x: float, *, T=None, H=None, p=None, Y=None) -> float:
        return self._scalar_getter(x, baseline_array=self.cp_map,
                                   field="cp_mass",
                                   T=T, H=H, p=p, Y=Y)

    def get_cv(self, x: float, *, T=None, H=None, p=None, Y=None) -> float:
        return self._scalar_getter(x, baseline_array=self.cv_map,
                                   field="cv_mass",
                                   T=T, H=H, p=p, Y=Y)

    def get_gamma(self, x: float, *, T=None, H=None, p=None, Y=None) -> float:
        if T is None and H is None and p is None and Y is None:
            return self._interp_baseline(self.gamma_map, x)
        T, H, p, Y = self._coalesce(x, T, H, p, Y)
        g = self._eval_state(T, H, p, Y)
        return g.cp_mass / g.cv_mass

    def get_H(self, x: float, *,
              T: Optional[float] = None,
              p: Optional[float] = None,
              Y: Optional[np.ndarray] = None) -> float:
        if T is None and p is None and Y is None:
            return self._interp_baseline(self.h_map, x)
        T, H, p, Y = self._coalesce(x, T, None, p, Y)
        return self._eval_state(T, H, p, Y).enthalpy_mass

    def get_rho(self, x: float, *, T=None, H=None, p=None, Y=None) -> float:
        return self._scalar_getter(x, baseline_array=self.rho_map,
                                   field="density",
                                   T=T, H=H, p=p, Y=Y)

    def get_mu(self, x: float, *, T=None, H=None, p=None, Y=None) -> float:
        return self._scalar_getter(x, baseline_array=self.mu_map,
                                   field="viscosity",
                                   T=T, H=H, p=p, Y=Y)

    def get_k(self, x: float, *, T=None, H=None, p=None, Y=None) -> float:
        return self._scalar_getter(x, baseline_array=self.k_map,
                                   field="thermal_conductivity",
                                   T=T, H=H, p=p, Y=Y)
    def get_Pr(self, x: float, *, T=None, H=None, p=None, Y=None) -> float:
        return self._scalar_getter(x, baseline_array=self.k_map,
                                   field="prandtl",
                                   T=T, H=H, p=p, Y=Y)

    def get_T_aw(self, x, gamma, M_inf, T_inf): 
        """ Adiabatic wall temperature, sometimes called recovery temperature T_r in some sources
        Unsure at what location this gamma is supposed to be evaluated at. Wall or freestream?"""
        Pr = self.get_Pr(x)
        r = Pr**(1/3)
        return T_inf*(1 + r*((gamma -1)/2)*M_inf**2)

class OptimalValues: 
    def __init__(self, ox, fu, F, MR, p_c, p_e, L_star, Isp=None):

        C = CEA_metric(oxName=ox, fuelName=fu,
            isp_units='sec',
            cstar_units='m/s',
            pressure_units='Pa',
            temperature_units='K',
            sonic_velocity_units='m/s',
            enthalpy_units='J/kg',
            density_units='kg/m^3',
            specific_heat_units='J/kg-K',
            viscosity_units='millipoise',
            thermal_cond_units='W/cm-degC',
            fac_CR=None,
            make_debug_prints=False)

        self.C = C

        # Derive variables
        eps = C.get_eps_at_PcOvPe(Pc=p_c, MR=MR, PcOvPe=p_c/p_e, frozen=0, frozenAtThroat=0)
        c_star = C.get_Cstar(Pc=p_c, MR=MR)
        if Isp == None: 
            Isp, _ = C.estimate_Ambient_Isp(Pc=p_c, MR=MR, eps=eps, Pamb=p_e, frozen=0, frozenAtThroat=0)
        
        #Isp = Isp*0.95

        # chamber properties: 
        rho_c, rho_t, rho_e = C.get_Densities(Pc=p_c, MR=MR, eps=eps)

        # calculated
        g = 9.81
        mdot = F/(Isp*g)
        mdot_fu = mdot/(1+MR)
        mdot_ox = mdot - mdot_fu

        A_t = c_star*mdot/p_c
        r_t = np.sqrt(A_t/np.pi)
        t_stay = L_star*A_t*rho_c/mdot
        V_c = mdot*t_stay/rho_c

        A_e = A_t*eps

        # input values
        
        self.MR = MR
        self.p_c = p_c
        self.p_e = p_e
        self.L_star = L_star

        # derived values
        self.Isp = Isp
        self.mdot = mdot
        self.mdot_fu = mdot_fu
        self.mdot_ox = mdot_ox
        self.t_stay = t_stay
        

        # optimal values
        self.eps_opt = eps
        self.A_t_opt = A_t
        self.r_t_opt = r_t
        self.V_c_opt = V_c
        self.A_e_opt = A_e


class OptimalValues2:
    # 1) keep a registry of (required_input_keys, solver_function)
    _solvers = []

    @classmethod
    def register_solver(cls, required_inputs):
        """
        Decorator: register a solver function that
        accepts exactly the given set of input names.
        """
        def decorator(fn):
            cls._solvers.append((set(required_inputs), fn))
            return fn
        return decorator

    def __init__(self, **kwargs):
        """
        kwargs must include at least one solver's required keys;
        and exactly one solver must match.
        """
        provided = set(kwargs)
        # find all solvers whose required_inputs is subset of provided
        matches = [
            fn for reqs, fn in self._solvers
            if reqs <= provided
        ]
        if not matches:
            raise ValueError(f"No solver for inputs {provided}")
        if len(matches) > 1:
            # ambiguous pathway
            sigs = [tuple(sorted(reqs)) for reqs, fn in self._solvers if reqs <= provided]
            raise ValueError(f"Ambiguous input sets {sigs}")
        # exactly one match
        solver = matches[0]
        # solver should return a dict of ALL 16 properties
        results = solver(**kwargs)
        self.__dict__.update(results)

@OptimalValues2.register_solver(
    required_inputs=['ox','fu','F','MR','p_c','p_e','L_star']
)

def _solver_exit_pressure(ox, fu, F, MR, p_c, p_e, L_star):

    C = CEA_metric(oxName=ox, fuelName=fu,
        isp_units='sec',
        cstar_units='m/s',
        pressure_units='Pa',
        temperature_units='K',
        sonic_velocity_units='m/s',
        enthalpy_units='J/kg',
        density_units='kg/m^3',
        specific_heat_units='J/kg-K',
        viscosity_units='millipoise',
        thermal_cond_units='W/cm-degC',
        fac_CR=None,
        make_debug_prints=False)

    # Derive variables
    eps = C.get_eps_at_PcOvPe(Pc=p_c, MR=MR, PcOvPe=p_c/p_e, frozen=0, frozenAtThroat=0)
    c_star = C.get_Cstar(Pc=p_c, MR=MR)
    print(f"c_star: {c_star}")
    Isp, _ = C.estimate_Ambient_Isp(Pc=p_c, MR=MR, eps=eps, Pamb=p_e, frozen=0, frozenAtThroat=0)
    Isp_vac = C.get_Isp(Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)

    eta_Cf = 0.9517
    eta_c_star = 0.9892
    Isp_vac_new = Isp_vac*eta_Cf#*eta_c_star
    Isp_vac_new = 440.3
    print(f"ambient Isp: {Isp}")
    print(f"vacuum Isp:  {Isp_vac_new}")

    # chamber properties: 
    rho_c, rho_t, rho_e = C.get_Densities(Pc=p_c, MR=MR, eps=eps)

    # calculated
    g = 9.80665
    
    mdot = F/(Isp_vac_new*g)
    mdot_fu = mdot/(1+MR)
    mdot_ox = mdot - mdot_fu
    print(f"mdot_target: 16.946211")
    print(f"mdot: {mdot}")

    A_t = c_star*mdot/p_c
    r_t = np.sqrt(A_t/np.pi)
    t_stay = L_star*A_t*rho_c/mdot
    V_c = mdot*t_stay/rho_c

    A_e = A_t*eps

    return dict(ox=ox, fu=fu, F=F, MR=MR, p_c=p_c, p_e=p_e, L_star=L_star, 
                eps=eps, c_star=c_star, Isp=Isp_vac_new, mdot=mdot, mdot_fu=mdot_fu, 
                mdot_ox=mdot_ox, A_t=A_t, r_t=r_t, t_stay=t_stay, V_c=V_c, A_e=A_e)

@OptimalValues2.register_solver(
    required_inputs=['ox','fu','F','MR','p_c','eps','L_star']
)

def _solver_area_ratio(ox, fu, F, MR, p_c, eps, L_star):

    C = CEA_metric(oxName=ox, fuelName=fu,
        isp_units='sec',
        cstar_units='m/s',
        pressure_units='Pa',
        temperature_units='K',
        sonic_velocity_units='m/s',
        enthalpy_units='J/kg',
        density_units='kg/m^3',
        specific_heat_units='J/kg-K',
        viscosity_units='millipoise',
        thermal_cond_units='W/cm-degC',
        fac_CR=None,
        make_debug_prints=False)

    # Derive variables
    c_star = C.get_Cstar(Pc=p_c, MR=MR)
    Isp, _ = C.get_Isp(Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)
    pcovpe = C.get_PcOvPe(Pc=p_c, MR=MR, eps=eps)
    p_e = p_c/pcovpe
    eta_Cf = 0.9517
    eta_c_star = 0.9892

    # chamber properties: 
    rho_c, rho_t, rho_e = C.get_Densities(Pc=p_c, MR=MR, eps=eps)

    # calculated
    g = 9.81
    mdot = F/(Isp*g)
    mdot_fu = mdot/(1+MR)
    mdot_ox = mdot - mdot_fu

    A_t = c_star*mdot/p_c
    r_t = np.sqrt(A_t/np.pi)
    t_stay = L_star*A_t*rho_c/mdot
    V_c = mdot*t_stay/rho_c

    A_e = A_t*eps

    return dict(ox=ox, fu=fu, F=F, MR=MR, p_c=p_c, p_e=p_e, L_star=L_star, 
                eps=eps, c_star=c_star, Isp=Isp, mdot=mdot, mdot_fu=mdot_fu, 
                mdot_ox=mdot_ox, A_t=A_t, r_t=r_t, t_stay=t_stay, V_c=V_c, A_e=A_e)

@OptimalValues2.register_solver(
    required_inputs=['ox','fu','F','MR','p_c','A_e','L_star']
)

def _solver_exit_area(ox, fu, F, MR, p_c, A_e, L_star):
    C = CEA_metric(
        oxName=ox, fuelName=fu,
        isp_units='sec',
        cstar_units='m/s',
        pressure_units='Pa',
        temperature_units='K',
        sonic_velocity_units='m/s',
        enthalpy_units='J/kg',
        density_units='kg/m^3',
        specific_heat_units='J/kg-K',
        viscosity_units='millipoise',
        thermal_cond_units='W/cm-degC',
        fac_CR=None,
        make_debug_prints=False
    )

    # Pre‐compute everything that doesn’t change per eps
    c_star = C.get_Cstar(Pc=p_c, MR=MR)

    def residual(eps):
        # 2) compute performance at that eps/p_e
        Isp = C.get_Isp(Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)
        # 3) mass‐flow & throat area
        mdot = F/(Isp * 9.81)
        A_t  = c_star * mdot / p_c

        # 4) how far off is our guessed A_e?
        return (A_t * eps) - A_e

    # find eps in a reasonable bracket
    eps_opt = opt.brentq(residual, 1e-3, 200)

    # now that we have eps_opt, fully compute everything once
    pcovpe = C.get_PcOvPe(Pc=p_c, MR=MR, eps=eps_opt)
    p_e    = p_c / pcovpe
    Isp = C.get_Isp(Pc=p_c, MR=MR, eps=eps_opt, frozen=0, frozenAtThroat=0)

    rho_c, *_ = C.get_Densities(Pc=p_c, MR=MR, eps=eps_opt)
    mdot      = F/(Isp * 9.81)
    mdot_fu   = mdot/(1+MR)
    mdot_ox   = mdot - mdot_fu

    A_t       = c_star * mdot / p_c
    r_t       = np.sqrt(A_t/np.pi)
    t_stay    = L_star * A_t * rho_c / mdot
    V_c       = mdot * t_stay / rho_c

    return dict(
        ox=ox, fu=fu, F=F, MR=MR, p_c=p_c,
        p_e=p_e, L_star=L_star,
        eps=eps_opt, c_star=c_star, Isp=Isp,
        rho_c=rho_c, mdot=mdot, mdot_fu=mdot_fu,
        mdot_ox=mdot_ox, A_t=A_t, r_t=r_t,
        t_stay=t_stay, V_c=V_c, A_e=A_e,
    )


# -----------------------------------------------------------------------------
#                       Example usage / plotting helper
# -----------------------------------------------------------------------------

"""if __name__ == "__main__":
    import pyskyfire as psf  # assumes your library layout

    # ---- RL10‑like contour inputs ------------------------------------ #
    p_c = 32.7501e5
    p_e = 0.0377e5
    MR = 4.0
    L_star = 0.95
    r_c = 0.123
    ox = "O2"
    fu = "CH4"
    theta_conv = 25
    R_1f, R_2f, R_3f = 1.5, 3, 0.5
    length_fraction = 0.713
    F = 73.4e3
    T_fuel_in=111.6
    T_ox_in=90.1

    V_c = 0.011527815541804631
    r_t = 0.062149375684373835
    eps = 58.16915234505697

    RL10_xs, RL10_rs = psf.regen.contour.get_contour(
        V_c=V_c,
        r_t=r_t,
        area_ratio=eps,
        r_c=r_c,
        theta_conv=theta_conv,
        nozzle="rao",
        R_1f=R_1f,
        R_2f=R_2f,
        R_3f=R_3f,
        length_fraction=length_fraction,
    )

    RL10_contour = psf.regen.Contour(RL10_xs, RL10_rs, name="Replicate Contour")

    cr = CombustionTransport(ox=ox, fu=fu, F=F, phi=None, MR=MR, p_c=p_c, p_e=p_e, eps=eps, L_star=L_star, T_fuel_in=T_fuel_in, T_ox_in=T_ox_in, mech="gri30.yaml" )
    cr.compute_transport(contour=RL10_contour)

    xs = np.array(RL10_contour.xs)

    def compare_with_old(cr, num_nodes=100):
        Compute RocketCEA ('old') values on cr.x_domain and compare with new Cantera results.
        # ---------- RocketCEA objects (user‑unit aware) ------------------
        if cr.fu == "H2": 
            cea_fuel = "LH2"
        elif cr.fu == "CH4":
            cea_fuel = "CH4"
        C = CEA_metric(
            oxName="LOX", fuelName=cea_fuel,
            isp_units='sec',           cstar_units='m/s',
            pressure_units='Pa',       temperature_units='K',
            sonic_velocity_units='m/s', enthalpy_units='J/kg',
            density_units='kg/m^3',    specific_heat_units='J/kg-K',
            viscosity_units='millipoise',          # 1 mP = 1e‑4 Pa·s
            thermal_cond_units='W/cm-degC',        # 1 W/cm‑K = 100 W/m‑K
            fac_CR=None, make_debug_prints=False
        )
        K = CEA_imperial(oxName="LOX", fuelName=cea_fuel)   # for Mach routines

        p_c, MR = cr.p_c, cr.MR
        x_t, A_t, A_c = cr.contour.x_t, cr.contour.A_t, cr.contour.A_c
        x_vals = cr.x_domain

        # ---- chamber transport once, so we can reuse for x < x_t ----

        # storage
        (M_old, gamma_old, T_old, p_old,
        h_old, rho_old, cp_old,
        mu_old, k_old, Pr_old) = ([] for _ in range(10))

        Pr_old2 = []

        for x in x_vals:
            A      = cr.contour.A(x)
            eps    = A / A_t
            sup    = x >  x_t
            throat = abs(x - x_t) < 1e-9

            # Mach number -------------------------------------------------
            if x < x_t:
                Mloc = K.get_Chamber_MachNumber(Pc=p_c*0.000145038, MR=MR, fac_CR=eps)
            elif throat:
                Mloc = 1.0
            else:
                Mloc = K.get_MachNumber(Pc=p_c*0.000145038, MR=MR, eps=eps)
            M_old.append(Mloc)

            # gamma -------------------------------------------------------
            if sup:
                _, gam = C.get_exit_MolWt_gamma(Pc=p_c, MR=MR, eps=eps)
            else:
                _, gam = C.get_Chamber_MolWt_gamma(Pc=p_c, MR=MR, eps=eps)   # chamber γ
            gamma_old.append(gam)

            # temperatures -----------------------------------------------
            if sup:
                _, _, Tloc = C.get_Temperatures(Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)
            else:
                Tloc = Tloc, _, _ = C.get_Temperatures(Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)
            T_old.append(Tloc)

            # pressures --------------------------------------------------
            if sup:
                PcOvPe  = C.get_PcOvPe(Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)
                ploc    = p_c / PcOvPe
            else:
                ploc = p_c
            p_old.append(ploc)

            # enthalpies -------------------------------------------------
            if sup:
                _, _, Hloc = C.get_Enthalpies(Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)
            else:
                Hloc = Hloc, _, _ = C.get_Enthalpies(Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)
            h_old.append(Hloc)

            # densities --------------------------------------------------
            if sup:
                _, _, rho = C.get_Densities(Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)
            else:
                rho = rho, _, _ = C.get_Densities(Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)  # same as ideal gas interpolation used before
            rho_old.append(rho)

            # transport (cp, μ, k, Pr) -----------------------------------
            if sup:
                cp_ex, mu_mP, k_WcmK, Pr = C.get_Exit_Transport(
                    Pc=p_c, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)  # :contentReference[oaicite:1]{index=1}
                mu_ex = mu_mP * 1e-4
                k_ex  = k_WcmK * 100.0
                Pr_ex = Pr
            else:
                cp_ex, mu_ch_mP, k_ch_WcmK, Pr_ch = C.get_Chamber_Transport(Pc=p_c, MR=MR)  # :contentReference[oaicite:0]{index=0}
                mu_ex = mu_ch_mP * 1e-4                # mP  → Pa·s
                k_ex  = k_ch_WcmK * 100.0              # W/cm‑K → W/m‑K
                Pr_ex = Pr_ch

            cp_old.append(cp_ex)
            mu_old.append(mu_ex)
            k_old.append(k_ex)
            Pr_old.append(Pr_ex)
            Pr_old2.append(cp_ex*mu_ex/k_ex)

        # ---- numpy arrays for plotting ---------------------------------
        (M_old, gamma_old, T_old, p_old,
        h_old, rho_old, cp_old,
        mu_old, k_old, Pr_old) = map(np.asarray,
            (M_old, gamma_old, T_old, p_old,
            h_old, rho_old, cp_old,
            mu_old, k_old, Pr_old))

        # ---- plot ------------------------------------------------------
        new_vals = [cr.M_map, cr.gamma_map, cr.T_map, cr.p_map, cr.h_map,  cr.cp_map,
                    cr.k_map, cr.mu_map, cr.Pr_map, cr.rho_map]
        old_vals = [M_old, gamma_old, T_old, p_old, h_old, cp_old,
                    k_old, mu_old, Pr_old, rho_old]
        labels   = ["Mach", "γ", "T", "p", "h", "cp",
                    "k", "μ", "Pr", "ρ"]
        units    = ["",    "",  "[K]", "[Pa]", "[J/kg]", "[J/kg K]",
                    "[W/m K]", "[Pa·s]", "‑", "[kg/m³]"]

        nplots   = len(labels)
        ncols    = 2
        nrows    = math.ceil(nplots / ncols)

        fig, axs = plt.subplots(nrows, ncols, figsize=(10, 2.4*nrows),
                                sharex=True)
        axs      = axs.flatten()

        for i in range(nplots):
            axs[i].plot(x_vals, new_vals[i], label="new (Cantera)")
            axs[i].plot(x_vals, old_vals[i], "--", label="old (CEA)")
            axs[i].set_ylabel(f"{labels[i]} {units[i]}")
            axs[i].legend(loc="best")

        for ax in axs:
            ax.set_xlabel("x [m]")

        plt.tight_layout()
        plt.show()
        
        plt.figure()
        plt.plot(x_vals, Pr_old)
        plt.plot(x_vals, Pr_old2)
        plt.plot(x_vals, cr.Pr)
        plt.show()

    compare_with_old(cr)"""

"""def get_full(C, Pc, MR, M, eps):
    print(f"MR: {MR}")
    full_output = C.get_full_cea_output(Pc=Pc/1e5, MR=MR, subar=None, PcOvPe=None, frozen=0, frozenAtThroat=0, short_output=0, show_transport=1, pc_units='bar', output="metric", show_mass_frac=True, fac_CR=eps)
    #res = extract_full_cea_properties(full_output)
    #res2 = print_full_output_lines(full_output)
    #res3 = parse_and_display_properties(full_output)
    props = parse_cea_properties(full_output)
    for k, v in props.items():
        print(f"{k:12s}: {v}")"""



def _find_anchor_line(lines, keyword):
    """Return index of *first* line containing *keyword* (case‑sensitive)."""
    for i, ln in enumerate(lines):
        if keyword in ln:
            return i
    return -1


def _parse_numbers(line):
    """Extract **exactly four** numeric values from a CEA table row.

    CEA’s format quirks
    -------------------
    A value looks like one of these:
      * ``1.5444 5``   → 1.5444 × 10⁵   (base and exponent as *separate* tokens)
      * ``8.7415-1``   → 8.7415 × 10⁻¹  (negative exponent *glued* to base)
      * ``9.4062 4``   → 9.4062 × 10⁴
    where the *base* always contains a decimal point.
    """
    tokens = line.split()
    values = []
    i = 0
    while i < len(tokens) and len(values) < 4:
        tok = tokens[i]
        # Candidate base token must contain a decimal point
        if '.' in tok:
            # Case 1: exponent glued to the base (look for trailing +/‑digits)
            m = re.match(r'^(-?\d*\.\d+)([+-]\d+)?$', tok)
            if not m:
                i += 1
                continue  # not a proper number token, skip
            base = float(m.group(1))
            exp = 0
            if m.group(2):
                exp = int(m.group(2))
            else:
                # Case 2: exponent is the next standalone integer token
                if i + 1 < len(tokens) and re.match(r'^-?\d+$', tokens[i + 1]):
                    exp = int(tokens[i + 1])
                    i += 1  # consume exponent token
            values.append(base * (10 ** exp))
        i += 1

    return values if len(values) == 4 else None

# -----------------------------------------------------------------------------
# Main parser
# -----------------------------------------------------------------------------

def parse_cea_properties(full_output):
    """Parse a RocketCEA *get_full_cea* text dump and return a dict:

        {
          'pressure'    : [P_injector, P_comb_end, P_throat, P_exit],
          'temperature' : [...],
          ... etc. ...
        }

    Algorithm
    ---------
    1.  Anchor on the first ``INJECTOR`` line, then capture the FIRST matching
        property row after that anchor for each requested variable.
    2.  Anchor on the first ``WITH EQUILIBRIUM REACTIONS`` line to grab
        conductivity & Prandtl rows.
    3.  Use :func:`_parse_numbers` to extract four numeric columns from each
        captured row, handling CEA’s implicit exponents.
    """
    lines = full_output.splitlines()

    #for l in lines:
    #    print(l)
    #input()

    # ------------------------------------------------------------------
    # 1) Locate anchors
    # ------------------------------------------------------------------
    inj_idx = _find_anchor_line(lines, "INJECTOR")
    if inj_idx == -1:
        raise ValueError("'INJECTOR' anchor not found in CEA output")

    eq_idx = _find_anchor_line(lines, "WITH EQUILIBRIUM REACTIONS")
    if eq_idx == -1:
        raise ValueError("'WITH EQUILIBRIUM REACTIONS' block not found")

    # ------------------------------------------------------------------
    # 2) Regex patterns to identify property rows
    # ------------------------------------------------------------------
    chamber_patterns = {
        'pressure':     r"^\s*P,\s*BAR",
        'temperature':  r"^\s*T,\s*K",
        'density':      r"^\s*RHO",
        'enthalpy':     r"^\s*H,\s*(?:KJ|CAL)",
        'cp':           r"^\s*Cp,",              # static Cp row
        'gamma':        r"^\s*GAMMAs",
        'mach':         r"^\s*MACH NUMBER",
        'viscosity':    r"^\s*VISC",
    }

    transport_patterns = {
        'conductivity': r"^\s*CONDUCTIVITY",
        'prandtl':      r"^\s*PRANDTL NUMBER",
    }

    prop_lines = {}



    # ------------------------------------------------------------------
    # 3) Scan chamber‑end section (after INJECTOR) for each pattern
    # ------------------------------------------------------------------
    for key, pat in chamber_patterns.items():
        regex = re.compile(pat, re.I)
        for ln in lines[inj_idx + 1 :]:
            if regex.search(ln):
                prop_lines[key] = ln
                break
        else:
            raise ValueError(f"Could not locate line for '{key}'")


    # ------------------------------------------------------------------
    # 4) Scan equilibrium‑transport section (after WITH EQUILIBRIUM REACTIONS)
    # ------------------------------------------------------------------
    for key, pat in transport_patterns.items():
        regex = re.compile(pat, re.I)
        for ln in lines[eq_idx + 1 :]:
            if regex.search(ln):
                prop_lines[key] = ln
                break
        else:
            raise ValueError(f"Could not locate line for '{key}' in transport block")

    #print(prop_lines["conductivity"])
    #input()
    # ------------------------------------------------------------------
    # 5) Parse numeric columns for each captured row
    # ------------------------------------------------------------------
    parsed = {}
    for key, line in prop_lines.items():
        numbers = _parse_numbers(line)
        if numbers is None:
            raise ValueError(
                f"Failed to parse four numeric values from line for '{key}':\n{line}"
            )
        parsed[key] = numbers


    parsed["conductivity"] = [c * 0.1 for c in parsed["conductivity"]]
    parsed["viscosity"] = [c * 0.0001 for c in parsed["viscosity"]]
    parsed["pressure"] = [c * 1e5 for c in parsed["pressure"]]
    parsed["cp"] = [c * 1000 for c in parsed["cp"]]
    parsed["cv"] = [cp_i / g_i for cp_i, g_i in zip(parsed["cp"], parsed["gamma"])]
    parsed["enthalpy"] = [c * 1000 for c in parsed["enthalpy"]]

    # ───────────────────── 3) Parse MASS FRACTIONS block ──────────────────
    mass_idx = _find_anchor_line(lines, "MASS FRACTIONS")
    mass_fractions: dict[str, list[float]] = {}

    if mass_idx != -1:
        i = mass_idx + 1
        # skip the blank line immediately following the header
        while i < len(lines) and not lines[i].strip():
            i += 1

        # read species rows until the next blank line
        while i < len(lines) and lines[i].strip():
            row = lines[i]
            species = row.split(None, 1)[0].replace('*', '')  # remove radical marker
            values = _parse_numbers(row)
            if values:
                mass_fractions[species] = values
            i += 1
    else:
        # Some very old RocketCEA versions do not emit mass fractions
        mass_fractions = {}

    parsed['mass_fractions'] = mass_fractions

    return parsed

def scale_dict_lists(data: dict[str, list[float]], factor: float) -> dict[str, list[float]]:
    """
    Multiply each 4-element list in *data* by *factor* and
    return a new dict with the same keys.

    Example
    -------
    >>> d = {"p": [1, 2, 3, 4], "T": [10, 20, 30, 40]}
    >>> scale_dict_lists(d, 2)
    {'p': [2, 4, 6, 8], 'T': [20, 40, 60, 80]}
    """
    return {k: [v_i * factor for v_i in v] for k, v in data.items()}
