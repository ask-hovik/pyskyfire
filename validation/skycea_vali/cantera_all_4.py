import math
from typing import Sequence, Optional
from rocketcea.cea_obj_w_units import CEA_Obj as CEA_metric
from rocketcea.cea_obj import CEA_Obj as CEA_imperial

import numpy as np
import cantera as ct
import scipy.optimize as opt
import matplotlib.pyplot as plt

__all__ = ["CombustionReaction"]

# -----------------------------------------------------------------------------
#  Low‑level helpers – kept separate for re‑use/testing
# -----------------------------------------------------------------------------

def _mach_from_area(A: float, At: float, gamma: float, supersonic: bool) -> float:
    """Invert area‑ratio → Mach, picking sub‑ or supersonic branch."""

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

def flame_T(P, T_fuel, T_ox, phi, fuel_species, oxidizer_species, mech="gri30.yaml", T_lo=150.0, T_hi=4500.0, tol=0.5, max_iter=60):
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
            raise ValueError(f"No real‑fluid EOS for {species}")
        s = eos()
        s.TP = T, P
        return s.enthalpy_mass  # J/kg

    H_react = m_fuel * _h_real(fuel_species, T_fuel, P) + m_ox * _h_real(oxidizer_species, T_ox, P)
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

class CombustionReaction:
    """Quasi-1-D equilibrium nozzle flow solver."""

    def __init__(self, 
                 ox: str, 
                 fu: str, 
                 F: float, 
                 p_c: float, 
                 p_e: float, 
                 L_star: float, 
                 contour, 
                 MR = None,
                 phi = None,
                 mech: str = "gri30.yaml", 
                 T_fuel: float = 298.15, 
                 T_ox: float = 298.15, 
                 tol: float = 1e-4, 
                 max_iter: int = 25):
        
        self.ox, self.fu = ox, fu
        self.p_c, self.p_e = float(p_c), float(p_e)
        self.L_star = L_star
        self.contour = contour
        self.gas = ct.Solution(mech)
        self.mech = mech
        
        if MR == None and phi != None:
            self.MR = self._find_MR_from_phi(phi)
            self.phi = phi
        elif phi == None and MR != None: 
            self.phi = self._find_phi_from_MR(MR)
            self.MR = MR
        else:
            raise ValueError("Either MR or phi must be supplied to CombustionReaction")

        self._solve_chamber(T_fuel, T_ox)
        self._solve_nozzle(tol, max_iter)

    def _find_MR_from_phi(self, phi):
        tmp = ct.Solution(self.mech)
        tmp.set_equivalence_ratio(1.0, f"{self.fu}:1", f"{self.ox}:1")
        Y = tmp.Y; m_f = Y[tmp.species_index(self.fu)]
        MR_stoich = (1.0 - m_f) / m_f
        return MR_stoich/phi

    def _find_phi_from_MR(self, MR): 
        tmp = ct.Solution(self.mech)
        tmp.set_equivalence_ratio(1.0, f"{self.fu}:1", f"{self.ox}:1")
        Y = tmp.Y; m_f = Y[tmp.species_index(self.fu)]
        MR_stoich = (1.0 - m_f) / m_f
        return MR_stoich/MR

    def _solve_chamber(self, T_fuel, T_ox):
        self.T_c = flame_T(self.p_c, T_fuel, T_ox, self.phi, self.fu, self.ox, self.mech)
        self.gas.set_equivalence_ratio(self.phi, f"{self.fu}:1", f"{self.ox}:1")
        self.gas.TP = self.T_c, self.p_c
        self.gas.equilibrate("TP")

        self.cp_c, self.cv_c, self.gamma_c = cp_cv_eq(self.gas, self.T_c, self.p_c)
        self.h_c = self.gas.enthalpy_mass

    def _solve_nozzle(self, tol, max_iter, n_nodes = 100):
        x_min = self.contour.xs[0]
        x_max = self.contour.xs[-1]
        # You had x_domain going from x_max -> x_min before; you can keep that if you like,
        # but for interpolation we typically store these in ascending order:
        self.x_domain = np.linspace(x_min, x_max, n_nodes)

        x_t = self.contour.x_t; 
        At = self.contour.A_t

        n = len(self.x_domain)
        gamma_old = np.full(n, self.gamma_c)

        # allocate
        self.M = np.empty(n) 
        self.gamma = np.empty(n)
        self.T = np.empty(n)
        self.p = np.empty(n)
        self.cp = np.empty(n)
        self.cv = np.empty(n)
        self.h = np.empty(n)
        self.k = np.empty(n)
        self.mu = np.empty(n)
        self.Pr = np.empty(n)
        self.rho = np.empty(n)

        C = CEA_metric(
            oxName=self.ox,
            fuelName=self.fu,
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
            make_debug_prints=False,
        )

        for _ in range(max_iter):
            for i, x in enumerate(self.x_domain):
                A_x = self.contour.A(x)
                supersonic = x >= x_t

                M_i = _mach_from_area(A_x, At, gamma_old[i], supersonic)
                T_i = self.T_c / (1.0 + 0.5 * (gamma_old[i] - 1.0) * M_i ** 2)
                p_i = self.p_c * (T_i / self.T_c) ** (gamma_old[i] / (gamma_old[i] - 1.0))

                self.gas.TP = T_i, p_i
                self.gas.equilibrate("TP")

                cp_i, cv_i, gamma_i = cp_cv_eq(self.gas, T_i, p_i)
 
                # store
                self.M[i] = M_i; 
                self.gamma[i] = gamma_i; 
                self.T[i] = T_i; 
                self.p[i] = p_i
                self.cp[i] = cp_i; 
                self.cv[i] = cv_i; 
                self.h[i] = self.gas.enthalpy_mass
                self.k[i] = self.gas.thermal_conductivity
                self.mu[i] = self.gas.viscosity
                self.Pr[i] = cp_i*self.gas.viscosity/self.gas.thermal_conductivity
                self.rho[i] = self.gas.density
                

            if np.max(np.abs(self.gamma - gamma_old)) < tol:
                return
            
            gamma_old = self.gamma
        raise RuntimeError("CombustionReaction: γ-M iteration did not converge")

# -----------------------------------------------------------------------------
#                       Example usage / plotting helper
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    import pyskyfire as psf  # assumes your library layout

    # ---- RL10‑like contour inputs ------------------------------------ #
    p_c = 32.7501e5
    p_e = 0.0377e5
    MR = 5.0
    L_star = 0.95
    r_c = 0.123
    ox = "O2"
    fu = "H2"
    theta_conv = 25
    R_1f, R_2f, R_3f = 1.5, 3, 0.5
    length_fraction = 0.713
    F = 73.4e3

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

    cr = CombustionReaction(ox=ox, fu=fu, F=F, phi=None, MR=MR, p_c=p_c, p_e=p_e, L_star=L_star, contour=RL10_contour)

    xs = np.array(RL10_contour.xs)

    def compare_with_old(cr, num_nodes=100):
        """Compute RocketCEA (‘old’) values on cr.x_domain and compare with new Cantera results."""
        # ---------- RocketCEA objects (user‑unit aware) ------------------
        C = CEA_metric(
            oxName="LOX", fuelName="LH2",
            isp_units='sec',           cstar_units='m/s',
            pressure_units='Pa',       temperature_units='K',
            sonic_velocity_units='m/s', enthalpy_units='J/kg',
            density_units='kg/m^3',    specific_heat_units='J/kg-K',
            viscosity_units='millipoise',          # 1 mP = 1e‑4 Pa·s
            thermal_cond_units='W/cm-degC',        # 1 W/cm‑K = 100 W/m‑K
            fac_CR=None, make_debug_prints=False
        )
        K = CEA_imperial(oxName=cr.ox, fuelName=cr.fu)   # for Mach routines

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
        new_vals = [cr.M, cr.gamma, cr.T, cr.p, cr.h,  cr.cp,
                    cr.k, cr.mu, cr.Pr, cr.rho]
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
        
        """plt.figure()
        plt.plot(x_vals, Pr_old)
        plt.plot(x_vals, Pr_old2)
        plt.plot(x_vals, cr.Pr)
        plt.show()"""

    compare_with_old(cr)
