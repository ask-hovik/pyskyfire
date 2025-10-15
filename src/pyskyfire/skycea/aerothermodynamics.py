from __future__ import annotations

import os
import math
import warnings
from typing import Optional

# Set proper trans.lib for CEA_Wrap
script_dir = os.path.dirname(os.path.abspath(__file__))
trans_path = os.path.join(script_dir, "data", r"trans.lib")
os.environ["CEA_TRANS_LIB"] = trans_path
import CEA_Wrap as cea

from dataclasses import dataclass, asdict
import numpy as np
from bisect import bisect_left


class Aerothermodynamics:  
    def __init__(self, optimum: dict[str, float | str], chemrep_map: Optional[dict[str, str]] = None):
        """Add all values in the optimum dict to self"""
        for key, value in optimum.items():
            setattr(self, key, value)
        self.optimum = optimum
        self.chemrep_map = chemrep_map or {}

    @classmethod
    def from_F_eps_Lstar(cls, fu, ox, MR, p_c, F, eps, L_star, T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15):
        """Calculate optimal values using thrust, exit pressure and L-star"""

        fus = []
        oxs = []
        for prop, frac in zip(fu.propellants, fu.fractions):
            fus.append(cea.Fuel(prop, wt=frac, temp=T_fu_in))

        for prop, frac in zip(ox.propellants, ox.fractions):
            oxs.append(cea.Oxidizer(prop, wt=frac, temp=T_ox_in))
        
        rp = cea.RocketProblem(o_f = MR, 
                               pressure=p_c*0.000145038, # Convert Pa to psi
                               materials=[*oxs, *fus], 
                               sup=eps) 
        R = rp.run()
        
        #names = cea.ThermoInterface.keys()
        #print(names)

        c_star = float(getattr(R, "cstar"))           
        Isp_ideal_amb = float(getattr(R, "isp"))      # perfectly expanded Isp
        Isp_vac = float(getattr(R, "ivac"))           
        rho_c = float(getattr(R, "c_rho"))
        T_c = float(getattr(R, "c_t"))
        T_t = float(getattr(R, "t_t"))
        p_t = float(getattr(R, "t_p"))

        g = 9.81
        p_SL = 1.01325e5 # sea level pressure
        mdot = F/(Isp_vac*g)
        mdot_fu = mdot/(1+MR)
        mdot_ox = mdot - mdot_fu

        A_t = c_star*mdot/p_c
        r_t = np.sqrt(A_t/np.pi)
        t_stay = L_star*A_t*rho_c/mdot
        V_c = mdot*t_stay/rho_c

        A_e = A_t*eps
        r_e = np.sqrt(A_e/np.pi)

        # calculate ambient Isp
        CF_vac = Isp_vac * g / c_star
        CF_amb = CF_vac - (p_amb / p_c) * (A_e/A_t)
        Isp_amb = CF_amb*c_star/g

        CF_SL = CF_vac - (p_SL / p_c) * (A_e/A_t)
        Isp_SL = CF_SL*c_star/g

        optimum = dict(
            fu=fu, ox=ox, MR=MR, p_c=p_c, p_t=p_t, T_c=T_c, T_t=T_t, F=F, eps=eps, L_star=L_star, c_star=c_star,
            p_amb=p_amb, Isp_ideal_amb=Isp_ideal_amb, Isp_vac=Isp_vac, Isp_amb=Isp_amb,
            Isp_SL=Isp_SL, CF_vac=CF_vac, CF_amb=CF_amb, CF_SL=CF_SL, mdot=mdot,
            mdot_fu=mdot_fu, mdot_ox=mdot_ox, t_stay=t_stay, A_t=A_t, A_e=A_e,
            r_t=r_t, r_e=r_e, V_c=V_c, T_fu_in=T_fu_in, T_ox_in=T_ox_in, npts=npts,
        )

        return cls(optimum)

    def compute_aerothermodynamics(self, contour, Nt: int = 64):
        """Create 2-D property maps on (x, T). Column 0 = equilibrium at that x."""
        # ---- x grid ----
        self.x_nodes = np.linspace(contour.xs[0], contour.xs[-1], self.npts)
        A_t = float(contour.A_t)
        eps_nodes = np.array([contour.A(x) / A_t for x in self.x_nodes])
        Nx = len(self.x_nodes)
        Nt = len(self.x_nodes)

        # ---- allocate maps: (Nx, Nt) ----
        self.Nt = Nt
        shape = (Nx, Nt)
        self.M_map     = np.empty(shape)
        self.T_map     = np.empty(shape)
        self.p_map     = np.empty(shape)
        self.rho_map   = np.empty(shape)
        self.cp_map    = np.empty(shape)
        self.gamma_map = np.empty(shape)
        self.h_map     = np.empty(shape)
        self.a_map     = np.empty(shape)
        self.mu_map    = np.empty(shape)
        self.k_map     = np.empty(shape)
        self.Pr_map    = np.empty(shape)
        self.X_map     = []              # still per-x dict (equilibrium)
        self.T_grid    = np.empty(shape) # per-row temperature grid

        # convenience aliases
        fu = self.fu
        ox = self.ox
        fus = []
        oxs = []
        for prop, frac in zip(fu.propellants, fu.fractions):
            fus.append(cea.Fuel(prop, wt=frac, temp=self.T_fu_in))

        for prop, frac in zip(ox.propellants, ox.fractions):
            oxs.append(cea.Oxidizer(prop, wt=frac, temp=self.T_ox_in))
        

        for i, (x, eps) in enumerate(zip(self.x_nodes, eps_nodes)):
            output = True
            if output is True:
                print(f'\rPrecomputing map: {math.ceil(i/self.npts*100)}%', end='', flush=True)

            # ---- equilibrium at this x (column 0) ----
            if x < 0:
                rp = cea.RocketProblem(
                    o_f=self.MR, pressure=self.p_c*0.000145038, materials=[*fus, *oxs], sub=eps)
                R = rp.run()
                self.X_map.append(getattr(R, "prod_c", {}))
            else:
                rp = cea.RocketProblem(
                    o_f=self.MR, pressure=self.p_c*0.000145038, materials=[*fus, *oxs], sup=eps)
                R = rp.run()
                self.X_map.append(getattr(R, "prod_e", {}))

            # store equilibrium (col 0)
            self.M_map[i, 0]     = getattr(R, "mach", None)
            self.T_map[i, 0]     = getattr(R, "t", None)
            self.p_map[i, 0]     = getattr(R, "p", None)
            self.rho_map[i, 0]   = getattr(R, "rho", None)
            self.cp_map[i, 0]    = getattr(R, "cp", None)
            self.gamma_map[i, 0] = getattr(R, "gamma", None)
            self.h_map[i, 0]     = getattr(R, "h", None)
            self.a_map[i, 0]     = getattr(R, "son", None)
            self.mu_map[i, 0]    = getattr(R, "visc", None)
            self.k_map[i, 0]     = getattr(R, "cond", None)
            self.Pr_map[i, 0]    = getattr(R, "pran", None)

            # ---- build this row’s temperature grid: T_eq(x) → 250 K ----
            T_eq_i = float(self.T_map[i, 0])   # equilibrium T at this x
            T_low  = 200
            # make it descending or ascending; we’ll treat as ascending for interp
            row_Ts = np.linspace(T_eq_i, T_low, Nt)
            self.T_grid[i, :] = row_Ts

            # ---- fill columns 1..Nt-1 with TP evaluations at fixed pressure ----
            p_Pa_i = float(self.p_map[i, 0]) * 1e5  # bar → Pa (CEA typically returns bar)
            for j in range(1, Nt):
                Tj = row_Ts[j]
                tp = cea.TPProblem(
                    pressure=p_Pa_i * 0.000145038,  # Pa → psi
                    temperature=Tj,
                    materials=[*fus, *oxs],
                    o_f=self.MR
                )
                Rt = tp.run()
                self.M_map[i, j]     = getattr(Rt, "mach", None)
                self.T_map[i, j]     = getattr(Rt, "t", None)
                self.p_map[i, j]     = getattr(Rt, "p", None)
                self.rho_map[i, j]   = getattr(Rt, "rho", None)
                self.cp_map[i, j]    = getattr(Rt, "cp", None)
                self.gamma_map[i, j] = getattr(Rt, "gamma", None)
                self.h_map[i, j]     = getattr(Rt, "h", None)
                self.a_map[i, j]     = getattr(Rt, "son", None)
                self.mu_map[i, j]    = getattr(Rt, "visc", None)
                self.k_map[i, j]     = getattr(Rt, "cond", None)
                self.Pr_map[i, j]    = getattr(Rt, "pran", None)

    def _bilinear_map(self, prop_map: np.ndarray, x: float, T: float) -> float | None:
        xs = self.x_nodes
        if not (xs[0] <= x <= xs[-1]):
            return None

        i1 = bisect_left(xs, x)
        if i1 == 0:
            i0 = 0; wx = 0.0
        elif i1 >= len(xs):
            i0 = len(xs) - 1; wx = 0.0
        else:
            i0 = i1 - 1
            x0, x1 = xs[i0], xs[i1]
            wx = (x - x0) / (x1 - x0)

        def interp_row(i):
            Tr = self.T_grid[i, :]
            Zr = prop_map[i, :]
            # ensure ascending T for np.interp
            if Tr[0] > Tr[-1]:
                Tr = Tr[::-1]
                Zr = Zr[::-1]
            if not (Tr[0] <= T <= Tr[-1]):
                return None
            return float(np.interp(T, Tr, Zr))

        v0 = interp_row(i0)
        v1 = interp_row(i1)
        if (v0 is None) or (v1 is None):
            return None
        return (1.0 - wx) * v0 + wx * v1

        printout = False # TODO: remove printout when stable
        if printout:
            print(f"T:     \n{self.T_map}")
            print(f"p:     \n{self.p_map}")
            print(f"rho:   \n{self.rho_map}")
            print(f"cp:    \n{self.cp_map}")
            print(f"gamma: \n{self.gamma_map}")
            print(f"h:     \n{self.h_map}")
            print(f"a:     \n{self.a_map}")
            print(f"mu:    \n{self.mu_map}")
            print(f"k:     \n{self.k_map}")
            print(f"X_map:     \n{self.X_map}")   
            s = input() 

    # ---------- Getter functionality ----------
    def _interp_scalar(self, x: float, xs: np.ndarray, ys: np.ndarray) -> float:
        """Linear interpolation with endpoint clamping."""
        if x <= xs[0]:  return float(ys[0])
        if x >= xs[-1]: return float(ys[-1])
        i = bisect_left(xs, x)
        x0, x1 = xs[i-1], xs[i]
        y0, y1 = ys[i-1], ys[i]
        w = (x - x0) / (x1 - x0)
        return float(y0 * (1.0 - w) + y1 * w)

    def _interp_X_dict(self, x: float) -> dict[str, float]:
        """
        Interpolate mole-fraction dict at position x from self.X_map (list[dict]).
        Uses the union of species at the two bracketing nodes, fills missing with 0,
        clips negatives, and renormalizes to sum = 1.
        """
        xs = self.x_nodes
        Xlist = self.X_map
        if not len(xs) or not len(Xlist):
            return {}

        # Clamp OOB to nearest endpoint
        if x <= xs[0]:
            d = Xlist[0] or {}
            S = sum(d.values()) or 1.0
            return {k: max(0.0, v)/S for k, v in d.items()}
        if x >= xs[-1]:
            d = Xlist[-1] or {}
            S = sum(d.values()) or 1.0
            return {k: max(0.0, v)/S for k, v in d.items()}

        # Bracket + blend
        i = bisect_left(xs, x)
        x0, x1 = xs[i-1], xs[i]
        w = (x - x0) / (x1 - x0)
        d0 = Xlist[i-1] or {}
        d1 = Xlist[i]   or {}
        species = set(d0) | set(d1)

        out = {}
        for s in species:
            v = (1.0 - w) * d0.get(s, 0.0) + w * d1.get(s, 0.0)
            out[s] = 0.0 if v < 0.0 else v

        # Renormalize
        S = sum(out.values())
        if S > 0.0:
            for s in out:
                out[s] /= S
        return out

    def _evaluate_tp(self, *, T: float, p: float):
        """
        Equilibrium evaluation at (T, p) using your reactants and MR (no frozen option).
        """
        #print(f"\ntemp: {T}")
        tp = cea.TPProblem(
            pressure=p*0.000145038, # Pa to psi
            temperature=T,
            materials=[*self.fu, *self.ox],
            o_f=self.MR
        )
        warnings.warn("Computation outside precomputed range")
        return tp.run()

    def _interp_eq_column(self, map2d: np.ndarray, x: float) -> float:
        """Interpolate along x using the equilibrium column (col 0)."""
        return self._interp_scalar(x, self.x_nodes, map2d[:, 0])

    def _evaluate_hp(self, *, h_target_Jkg: float, p: float):
        # STILL UNDER CONSTRUCTION
        """
        Evaluate an equilibrium HP state at local pressure p [Pa] and mixture specific enthalpy h [J/kg].
        Uses the 'chemical_representation + hf' override path.
        Strategy: assign zero enthalpy to all reactants EXCEPT one 'adjuster' fuel species,
        whose molar enthalpy is chosen so that the mixture enthalpy equals h_target_Jkg.

        Requirements:
          - self.chemrep_map must provide exploded formulas for any species we override.
          - We'll pick the FIRST fuel in self.fu as the 'adjuster'.
        """
        if not self.fu:
            raise RuntimeError("No fuel materials found.")
        adjuster = self.fu[0]
        adj_name = getattr(adjuster, "name", None) or getattr(adjuster, "NAME", None)
        if adj_name is None:
            # CEA_Wrap objects typically have .name
            raise RuntimeError("Cannot determine name of first fuel for enthalpy adjustment.")

        # ---- build overall reactant mass fractions (normalized) ----
        # Group-internal weights:
        fu_w = np.array([getattr(m, "wt", 0.0) for m in self.fu], dtype=float)
        ox_w = np.array([getattr(m, "wt", 0.0) for m in self.ox], dtype=float)
        if fu_w.sum() <= 0 or ox_w.sum() <= 0:
            raise RuntimeError("Fuel/oxidizer component weights must be positive.")

        fu_w /= fu_w.sum()
        ox_w /= ox_w.sum()

        # Overall stream mass fractions given MR = m_ox / m_fu:
        w_fu_stream = 1.0 / (1.0 + self.MR)
        w_ox_stream = self.MR / (1.0 + self.MR)
        w_overall = []
        names = []
        is_fuel_flags = []
        for (m, wf) in zip(self.fu, fu_w):
            nm = getattr(m, "name", None) or getattr(m, "NAME", None)
            names.append(nm); is_fuel_flags.append(True)
            w_overall.append(w_fu_stream * wf)
        for (m, wo) in zip(self.ox, ox_w):
            nm = getattr(m, "name", None) or getattr(m, "NAME", None)
            names.append(nm); is_fuel_flags.append(False)
            w_overall.append(w_ox_stream * wo)
        w_overall = np.asarray(w_overall, dtype=float)

        # ---- compute adjuster molar mass from exploded formula ----
        if adj_name not in self.chemrep_map:
            raise KeyError(f"Missing chemical_representation for adjuster '{adj_name}'. "
                           f"Add to self.chemrep_map.")
        mw_adjuster = _mw_from_exploded(self.chemrep_map[adj_name])  # kg/mol

        # Mass fraction of the adjuster in the overall mixture:
        try:
            idx_adj = names.index(adj_name)  # first occurrence among fuels
        except ValueError:
            raise RuntimeError("Adjuster fuel name not found in assembled reactant list.")
        w_adj = w_overall[idx_adj]
        if w_adj <= 0.0:
            raise RuntimeError("Adjuster species has zero mass fraction; cannot carry enthalpy.")

        # ---- choose assigned molar enthalpies H_i (J/mol) ----
        # All others = 0; adjuster gets H_adj so that h_mix = w_adj * H_adj / mw_adj:
        H_adj = h_target_Jkg * mw_adjuster / w_adj  # J/mol
        assigned_H = np.zeros_like(w_overall)
        assigned_H[idx_adj] = H_adj

        # ---- construct CEA_Wrap materials with chemical_representation + hf ----
        # We preserve inlet temperatures from self.fu/self.ox objects.
        materials = []
        # fuels first, keeping original order
        for (m, wf) in zip(self.fu, fu_w):
            nm = getattr(m, "name", None) or getattr(m, "NAME", None)
            Tm = getattr(m, "temp", 298.15)
            if nm not in self.chemrep_map:
                raise KeyError(f"Missing chemical_representation for fuel '{nm}'.")
            hf = H_adj if (nm == adj_name) else 0.0
            materials.append(cea.Fuel(
                name=nm, wt=float(wf*100.0), temp=float(Tm),
                chemical_representation=self.chemrep_map[nm], hf=float(hf)
            ))
        # oxidizers next
        for (m, wo) in zip(self.ox, ox_w):
            nm = getattr(m, "name", None) or getattr(m, "NAME", None)
            Tm = getattr(m, "temp", 298.15)
            if nm not in self.chemrep_map:
                raise KeyError(f"Missing chemical_representation for oxidizer '{nm}'.")
            materials.append(cea.Oxidizer(
                name=nm, wt=float(wo*100.0), temp=float(Tm),
                chemical_representation=self.chemrep_map[nm], hf=0.0
            ))

        # ---- run HP at local pressure ----
        hp = cea.HPProblem(
            pressure=p * 0.000145038,   # Pa -> psi
            materials=materials,
            o_f=self.MR
        )
        return hp.run()

    def _get_prop(self, *, x: float, T: float | None, h: float | None, map_attr: str, res_attr: str) -> float:
        if (T is not None) and (h is not None):
            raise ValueError("Provide only one of T or h.")
        prop_map = getattr(self, map_attr)

        # 1) No T/h: equilibrium column along x
        if (T is None) and (h is None):
            return self._interp_eq_column(prop_map, x)

        # Local static pressure from the equilibrium column (authoritative)
        p_pa = self.get_p(x)

        # 2) T provided: try map, else do a live TP at (T, p(x))
        if T is not None:
            v = self._bilinear_map(prop_map, x, T)
            if v is not None:
                return v
            R = self._evaluate_tp(T=T, p=p_pa)
            return float(getattr(R, res_attr, float("nan")))

        # 3) h provided (J/kg): run a single HP evaluation at (h, p(x))
        R = self._evaluate_hp(h_target_Jkg=h, p=p_pa)
        return float(getattr(R, res_attr, float("nan")))


    # ----------- public getters (maps + composition) -----------
    def get_T(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="T_map", res_attr="t")

    def get_p(self, x: float, T: float | None = None, h: float | None = None) -> float:
        # Pressure is map-only; ignore T/h by design.
        return self._interp_eq_column(self.p_map, x) * 1e5  # bar → Pa

    def get_M(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="M_map", res_attr="mach")

    def get_rho(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="rho_map", res_attr="rho")

    def get_cp(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="cp_map", res_attr="cp") * 1e3  # kJ/kg-K -> J/kg-K

    def get_gamma(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="gamma_map", res_attr="gamma")

    def get_h(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="h_map", res_attr="h") * 1e3    # kJ/kg -> J/kg

    def get_H(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="h_map", res_attr="h") * 1e3    # alias to J/kg

    def get_a(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="a_map", res_attr="son")

    def get_mu(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="mu_map", res_attr="visc")

    def get_k(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="k_map", res_attr="cond")

    def get_Pr(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(x=x, T=T, h=h, map_attr="Pr_map", res_attr="pran")

    def get_X(self, x: float) -> dict[str, float]:
        """Interpolated mole-fraction dict at position x."""
        return self._interp_X_dict(x)


# Map "public" property keys -> (CEA result attribute, unit_scale)
# unit_scale applies only when we are extracting from a TP run (not from maps).
_PROP_SPEC = {
    "T":     ("t",    1.0),      # K
    "rho":   ("rho",  1.0),      # kg/m^3
    "cp":    ("cp",   1e3),      # kJ/kg-K -> J/kg-K
    "gamma": ("gamma",1.0),
    "h":     ("h",    1e3),      # kJ/kg -> J/kg
    "H":     ("h",    1e3),      # alias of h
    "a":     ("son",  1.0),      # m/s
    "mu":    ("visc", 1.0),      # Pa·s
    "k":     ("cond", 1.0),      # W/m-K
    "Pr":    ("pran", 1.0),
    # "p" is handled specially (we take it from the precomputed map, in Pa)
}

# --- minimal periodic table for MW (kg/mol). Extend as needed. ---
_PERIODIC = {
    "H": 1.00794e-3, "C": 12.0107e-3, "N": 14.0067e-3, "O": 15.9994e-3,
    "F": 18.9984e-3, "Cl": 35.453e-3, "Ar": 39.948e-3, "He": 4.0026e-3,
    "Ne": 20.1797e-3, "S": 32.065e-3, "Si": 28.0855e-3, "B": 10.811e-3,
}
def _mw_from_exploded(chemrep: str) -> float:
    """
    Compute molar mass [kg/mol] from an exploded formula like 'C 2 H 6 O 1'.
    Supports one- or two-letter symbols; coefficients must be integers.
    """
    toks = chemrep.split()
    if len(toks) % 2 != 0:
        raise ValueError(f"Bad chemical_representation: {chemrep!r}")
    i = 0; mw = 0.0
    while i < len(toks):
        sym = toks[i]; n = int(float(toks[i+1])); i += 2
        if sym not in _PERIODIC:
            raise KeyError(f"Element {sym!r} not in periodic table; extend _PERIODIC.")
        mw += _PERIODIC[sym] * n
    return mw
