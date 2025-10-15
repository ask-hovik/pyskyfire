# nozzle_equil_solver.py
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

# You need: pip install cantera scipy
import cantera as ct
from scipy.optimize import root
from scipy.optimize import least_squares
from cantera import CanteraError
from typing import Optional, Union, Dict

@dataclass
class EquilNozzleInputs:
    mech: str

    # EITHER supply a ready-made chamber Solution...
    chamber_solution: Optional[ct.Solution] = None

    # ...OR supply a chamber state explicitly (any one of X as str/dict/array is fine)
    chamber_T: Optional[float] = None
    chamber_p: Optional[float] = None
    chamber_X: Optional[Union[str, Dict[str, float], np.ndarray]] = None

    # Geometry / environment
    A_throat: float = 0.0              # [m^2] (required)
    A_exit: Optional[float] = None     # [m^2]
    p_amb: Optional[float] = None      # [Pa]

    # Numerics
    fd_eps: float = 1e-4               # for a_eq finite-diff

class EquilibriumNozzle:
    def __init__(self, inp: EquilNozzleInputs):
        self.inp = inp
        self.gas = ct.Solution(inp.mech)

        # --- Initialize chamber state ----------------------------------------
        if inp.chamber_solution is not None:
            # Copy the exact chamber state you computed (e.g., HP-equilibrated)
            gC = inp.chamber_solution
            self.gas.TPX = gC.T, gC.P, gC.X
        else:
            if inp.chamber_T is None or inp.chamber_p is None or inp.chamber_X is None:
                raise ValueError("Provide either chamber_solution, or (chamber_T, chamber_p, chamber_X).")
            self.gas.TPX = inp.chamber_T, inp.chamber_p, inp.chamber_X

        # IMPORTANT: Chamber should be an equilibrium state. If you passed HP-equilibrated
        # gas_c already, you can skip this; otherwise ensure chamber is TP-equilibrated:
        # self.gas.equilibrate("TP")

        # Stagnation invariants (V≈0 in chamber)
        self.h0 = self.gas.enthalpy_mass
        self.s0 = self.gas.entropy_mass

        # Solve throat first (sets mdot)
        self.throat = self.solve_throat()
        self.mdot = self.throat["rho"] * self.throat["a_eq"] * inp.A_throat

    # --- Utilities ------------------------------------------------------------
    def _set_elements_by_mixing_rule(self):
        """
        Construct an arbitrary species mixture that matches the requested element totals,
        then let equilibrium find a consistent composition. If you already have propellant
        streams, replace this with your mixing logic and equilibrium at (T_c, p_c).
        """
        # Start with all species present a tiny bit to avoid zeroing, then equilibrate.
        X = np.ones(self.gas.n_species)
        self.gas.TPX = self.inp.T_c, self.inp.p_c, X
        self.gas.equilibrate("TP")  # rough start; true element totals enforced below

        # If you *must* enforce exact element totals, you'd solve μ-based constrained problem;
        # practically, for propulsion we set the propellant mix upstream and equilibrate TP here.

    def _equilibrate_TP(self, T, p):
        self.gas.TP = T, p
        self.gas.equilibrate("TP")
        rho = self.gas.density
        h   = self.gas.enthalpy_mass
        s   = self.gas.entropy_mass
        cp  = self.gas.cp_mass
        cv  = self.gas.cv_mass
        # Cantera universal gas constant (J/kmol-K). Fallback for older versions:
        try:
            R_univ = ct.gas_constant
        except AttributeError:
            R_univ = ct.UniversalGasConstant
        Rmix = R_univ / self.gas.mean_molecular_weight  # J/(kg·K)
        gamma = cp / cv
        a_frozen = math.sqrt(gamma * Rmix * T)
        return dict(T=T, p=p, rho=rho, h=h, s=s, cp=cp, cv=cv,
                    gamma=gamma, Rmix=Rmix, a_frozen=a_frozen)

    def _a_eq(self, T, p):
        """
        Equilibrium speed of sound a_eq = sqrt((∂p/∂ρ)_s,eq).
        Approximate by small symmetric isentropic TP steps:
        - Move slightly in p, adjust T to keep s = s0 (via 1D root), re-equilibrating each time;
        - Compute dp/dρ at constant s from two neighboring points.
        """
        base = self._equilibrate_TP(T, p)
        s_target = base["s"]

        def T_for_s_at_p(p_target, T_seed):
            def f(Ttrial):
                st = self._equilibrate_TP(Ttrial, p_target)["s"]
                return st - s_target
            sol = root(lambda x: f(x[0]), x0=np.array([T_seed]))
            return float(sol.x[0])

        # Choose dp around p
        dp = max(self.inp.fd_eps * p, 1.0)  # at least 1 Pa
        # Adjust T+ and T- s.t. s = s_target
        T_plus  = T_for_s_at_p(p + dp, T)   # re-equilibrates inside
        plus = self._equilibrate_TP(T_plus, p + dp)
        T_minus = T_for_s_at_p(p - dp, T)
        minus = self._equilibrate_TP(T_minus, p - dp)

        drho = plus["rho"] - minus["rho"]
        dp2 = (p + dp) - (p - dp)
        dpdro_s = dp2 / drho
        aeq = math.sqrt(max(dpdro_s, 0.0))
        return aeq

    # --- Throat solve (M=1) ---------------------------------------------------
    def solve_throat(self):
        """
        Unknowns: (T*, p*) with constraints:
          1) s(T*,p*,eq) = s0
          2) h0 - h(T*,p*,eq) - a_eq(T*,p*)^2 / 2 = 0   [since M=1 ⇒ V = a_eq]
        """
        T0, p0 = self.inp.T_c, self.inp.p_c
        x0 = np.array([0.9*T0, 0.6*p0])  # decent seeds: cooler & lower p than chamber

        def F(x):
            T, p = float(x[0]), float(x[1])
            st = self._equilibrate_TP(T, p)
            aeq = self._a_eq(T, p)
            f1 = st["s"] - self.s0
            f2 = self.h0 - st["h"] - 0.5 * aeq * aeq
            return np.array([f1, f2], dtype=float)

        sol = root(F, x0, method="hybr")
        Tt, pt = sol.x
        st = self._equilibrate_TP(Tt, pt)
        aeq = self._a_eq(Tt, pt)
        return dict(
            T=Tt, p=pt, rho=st["rho"], h=st["h"], s=st["s"],
            a_eq=aeq, V=aeq, M=1.0
        )

    def solve_exit_at_area_ratio(self, Ae_over_At: float):
        Ae = Ae_over_At * self.inp.A_throat
        mdot = self.mdot
        h0, s0 = self.h0, self.s0

        # Seeds: start near throat, lower p, lower T
        Tseed = max(500.0, 0.75 * self.throat["T"])
        pseed = max(200.0, 0.25 * self.throat["p"])
        x0 = np.array([math.log(Tseed), math.log(pseed)], dtype=float)

        # Characteristic scales to balance residuals
        Hs = abs(h0) if abs(h0) > 1 else 1.0
        Ss = abs(s0) if abs(s0) > 1 else 1.0
        Ms = abs(mdot) if abs(mdot) > 1e-12 else 1.0  # not used directly but here for completeness

        def residuals(x):
            T = math.exp(float(x[0]))
            p = math.exp(float(x[1]))
            try:
                st = self._equilibrate_TP(T, p)
            except CanteraError:
                # push solver away from infeasible guesses
                return np.array([1e6, 1e6], dtype=float)

            # Energy → V from h0
            dh = h0 - st["h"]
            if dh <= 0:
                # nonphysical (imaginary V): penalize; keeps solver in the right region
                return np.array([1e3 + dh, 1e3], dtype=float)

            V = math.sqrt(2.0 * dh)
            # Residuals (scaled)
            r_isent = (st["s"] - s0) / Ss
            r_cont  = (st["rho"] * V * Ae - mdot) / (mdot if mdot != 0 else 1.0)
            return np.array([r_isent, r_cont], dtype=float)

        sol = least_squares(
            residuals, x0, method="trf",
            ftol=1e-12, xtol=1e-12, gtol=1e-12,
            max_nfev=4000
        )

        lnT, lnp = sol.x
        Te, pe = math.exp(lnT), math.exp(lnp)
        st = self._equilibrate_TP(Te, pe)
        Ve = math.sqrt(max(2.0 * (h0 - st["h"]), 0.0))
        aeq = self._a_eq(Te, pe)
        Me = Ve / max(aeq, 1e-12)

        return dict(
            T=Te, p=pe, V=Ve, M=Me, rho=st["rho"], h=st["h"], s=st["s"], a_eq=aeq,
            Ae=Ae, mdot=mdot,
            success=bool(sol.success), nfev=sol.nfev, cost=sol.cost
        )



    def thrust(self, exit_state):
        """1-D thrust (no loss model): F = mdot*V_e + (p_e - p_amb) A_e"""
        if self.inp.p_amb is None:
            raise ValueError("p_amb is required for thrust.")
        return self.mdot * exit_state["V"] + (exit_state["p"] - self.inp.p_amb) * exit_state["Ae"]
