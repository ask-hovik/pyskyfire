from __future__ import annotations

import os

#ident = "RocketCEA"
ident = "CEAexec"
#ident = "CEA_Wrap"

if ident == "RocketCEA": 
    # RocketCEA
    os.environ["CEA_TRANS_LIB"]  = r"C:\Users\askho\MyFiles\Area\flamme_aerospace\testing\trans_thermo\RocketCEA\trans.lib"
    os.environ["CEA_THERMO_LIB"] = r"C:\Users\askho\MyFiles\Area\flamme_aerospace\testing\trans_thermo\RocketCEA\thermo.lib"

elif ident == "CEAexec":
    # CEAexec
    os.environ["CEA_TRANS_LIB"] = r"C:\Users\askho\MyFiles\Area\flamme_aerospace\testing\trans_thermo\CEAexec\trans.lib"
    os.environ["CEA_THERMO_LIB"] = r"C:\Users\askho\MyFiles\Area\flamme_aerospace\testing\trans_thermo\CEAexec\thermo.lib"

elif ident == "CEA_Wrap":
    # CEA_Wrap
    os.environ["CEA_TRANS_LIB"] = r"C:\Users\askho\MyFiles\Area\flamme_aerospace\testing\trans_thermo\CEA_Wrap\trans.lib"
    os.environ["CEA_THERMO_LIB"] = r"C:\Users\askho\MyFiles\Area\flamme_aerospace\testing\trans_thermo\CEA_Wrap\thermo.lib"

import CEA_Wrap as cea

from dataclasses import dataclass, asdict
import numpy as np
from bisect import bisect_left

@dataclass
class EngineOptimum:
    fu: str
    ox: str
    MR: float
    p_c: float
    F: float
    eps: float
    L_star: float
    p_amb: float
    Isp_ideal_amb: float
    Isp_vac: float
    Isp_amb: float
    Isp_SL: float
    CF_vac: float
    CF_amb: float
    CF_SL: float
    mdot_fu: float
    mdot_ox: float
    t_stay: float
    A_t: float
    A_e: float
    r_t: float
    r_e: float
    V_c: float
    npts: float # number of evaluated points along the nozzle for the aerothermodynamic properties
    
class Aerothermodynamics: 

    def __init__(self, optimum: EngineOptimum):
        """Add all values in EngineOptimum to self"""
        for key, value in asdict(optimum).items():
            setattr(self, key, value)

    @classmethod
    def from_F_eps_Lstar(cls, fu, ox, T_fu_in, T_ox_in, MR, p_c, F, eps, L_star, p_amb, npts):
        """Calculate optimal values using thrust, exit pressure and L-star"""

        cea_fu = cea.Fuel(fu, temp=T_fu_in)
        cea_ox = cea.Oxidizer(ox, temp=T_ox_in)

        rp = cea.RocketProblem(o_f = MR, 
                               pressure=p_c*0.000145038, # Convert Pa to psi
                               materials=[cea_fu, cea_ox], 
                               sup=eps) 
        R = rp.run()

        c_star = float(getattr(R, "cstar"))           
        Isp_ideal_amb = float(getattr(R, "isp"))      # perfectly expanded Isp
        Isp_vac = float(getattr(R, "ivac"))           
        rho_c = float(getattr(R, "c_rho"))

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


        optimum = EngineOptimum(cea_fu, cea_ox, MR, p_c, F, eps, L_star, p_amb, 
                                Isp_ideal_amb, Isp_vac, Isp_amb, Isp_SL, 
                                CF_vac, CF_amb, CF_SL, mdot_fu, mdot_ox, 
                                t_stay, A_t, A_e, r_t, r_e, V_c, npts)

        return cls(optimum)

    @classmethod
    def from_F_pamb_Lstar():
        """Calculate optimal values using thrust, area ratio and L-star"""
        print("not implemented yet")

    def compute_aerothermodynamics(self, contour):
        """Create maps of properties along the engine contour"""
        self.x_nodes = np.linspace(contour.xs[0], contour.xs[-1], self.npts)
        A_t = float(contour.A_t) #Keep in mind that this is not _necessarily_ the same as self.A_t
        eps_nodes = np.array([contour.A(x) / A_t for x in self.x_nodes])

        # Pre-allocate memory for the arrays 
        self.T_map = np.empty_like(self.x_nodes)
        self.p_map = np.empty_like(self.x_nodes)
        self.rho_map = np.empty_like(self.x_nodes)
        self.cp_map = np.empty_like(self.x_nodes)
        self.gamma_map = np.empty_like(self.x_nodes)
        self.h_map = np.empty_like(self.x_nodes)
        self.a_map = np.empty_like(self.x_nodes)
        self.mu_map = np.empty_like(self.x_nodes)
        self.k_map = np.empty_like(self.x_nodes)
        self.X_map = []

        for i, (x, eps) in enumerate(zip(self.x_nodes, eps_nodes)):
            if x < 0:
                rp = cea.RocketProblem(
                    o_f=self.MR,
                    pressure=self.p_c*0.000145038,    # Convert Pa to psi
                    materials=[self.fu, self.ox],
                    sub=eps)
                R = rp.run()

                #print(R)

                self.T_map[i]     = getattr(R, "c_t", None)
                self.p_map[i]     = getattr(R, "c_p", None)
                self.rho_map[i]   = getattr(R, "c_rho", None)
                self.cp_map[i]    = getattr(R, "c_cp", None)
                self.gamma_map[i] = getattr(R, "c_gamma", None)
                self.h_map[i]     = getattr(R, "c_h", None)
                self.a_map[i]     = getattr(R, "c_son", None)
                self.mu_map[i]    = getattr(R, "c_visc", None)
                self.k_map[i]     = getattr(R, "c_cond", None)
                self.X_map.append(getattr(R, "prod_c", {}))   # mole fractions → X

            elif x >= 0:
                rp = cea.RocketProblem(
                    o_f=self.MR,
                    pressure=self.p_c*0.000145038,    # Convert Pa to psi
                    materials=[self.fu, self.ox],
                    sup=eps)
                
                R = rp.run()

                self.T_map[i]     = getattr(R, "t", None)
                self.p_map[i]     = getattr(R, "p", None)
                self.rho_map[i]   = getattr(R, "rho", None)
                self.cp_map[i]    = getattr(R, "cp", None)
                self.gamma_map[i] = getattr(R, "gamma", None)
                self.h_map[i]     = getattr(R, "h", None)
                self.a_map[i]     = getattr(R, "son", None)
                self.mu_map[i]    = getattr(R, "visc", None)
                self.k_map[i]     = getattr(R, "cond", None)
                self.X_map.append(getattr(R, "prod_e", {}))   # mole fractions → X

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
        tp = cea.TPProblem(
            pressure=p,
            temperature=T,
            materials=[self.fu, self.ox],
            o_f=self.MR
        )
        return tp.run()

    def _get_prop(self, *, x: float, T: float | None, map_attr: str, res_attr: str) -> float:
        """
        Generic getter:
        - if T is None → interpolate precomputed map
        - else → equilibrium TPProblem at (T, p(x)) and return `res_attr`
        """
        xs = self.x_nodes
        ymap = getattr(self, map_attr)
        if T is None:
            return self._interp_scalar(x, xs, ymap)
        
        p = self._interp_scalar(x, xs, self.p_map)  # in whatever unit your map stores
        R = self._evaluate_tp(T=T, p=p)
        return float(getattr(R, res_attr, float("nan")))

    # ----------- public getters (maps + composition) -----------
    def get_T(self, x: float, T: float | None = None) -> float:
        """
        If T is None: interpolate mapped temperature.
        If T is provided: equilibrium TP at (T, p(x)) and return result 't'.
        """
        return self._get_prop(x=x, T=T, map_attr="T_map", res_attr="t")

    def get_p(self, x: float, T: float | None = None) -> float:
        """
        Pressure getter (interpolated only)
        """
        return self._interp_scalar(x, self.x_nodes, self.p_map)*1e5 # correct bar -> Pa

    def get_rho(self, x: float, T: float | None = None) -> float:
        return self._get_prop(x=x, T=T, map_attr="rho_map", res_attr="rho")

    def get_cp(self, x: float, T: float | None = None) -> float:
        return self._get_prop(x=x, T=T, map_attr="cp_map", res_attr="cp")*1e3

    def get_gamma(self, x: float, T: float | None = None) -> float:
        return self._get_prop(x=x, T=T, map_attr="gamma_map", res_attr="gamma")

    def get_h(self, x: float, T: float | None = None) -> float:
        return self._get_prop(x=x, T=T, map_attr="h_map", res_attr="h")

    def get_a(self, x: float, T: float | None = None) -> float:
        return self._get_prop(x=x, T=T, map_attr="a_map", res_attr="son")

    def get_mu(self, x: float, T: float | None = None) -> float:
        return self._get_prop(x=x, T=T, map_attr="mu_map", res_attr="visc")

    def get_k(self, x: float, T: float | None = None) -> float:
        return self._get_prop(x=x, T=T, map_attr="k_map", res_attr="cond")

    def get_X(self, x: float) -> dict[str, float]:
        """Interpolated mole-fraction dict at position x."""
        return self._interp_X_dict(x)
