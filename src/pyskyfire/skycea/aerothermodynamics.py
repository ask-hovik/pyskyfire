"""CEA-backed gas properties along a rocket-engine contour.

All public values use base SI units. The old constructor and getter API is
retained, but the implementation has only one design path and one state path.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field

import cea
import numpy as np

G0 = 9.81
SEA_LEVEL_PRESSURE = 1.01325e5
PA_TO_BAR = 1e-5
BAR_TO_PA = 1e5
KJ_TO_J = 1e3
MILLIPOISE_TO_PA_S = 1e-4
MW_PER_CM_K_TO_W_PER_M_K = 0.1


@dataclass(slots=True)
class GasState:
    T: float = math.nan
    p: float = math.nan
    M: float = math.nan
    rho: float = math.nan
    cp: float = math.nan
    gamma: float = math.nan
    h: float = math.nan
    a: float = math.nan
    mu: float = math.nan
    k: float = math.nan
    Pr: float = math.nan
    mw: float = math.nan
    X: dict = field(default_factory=dict)

    def __repr__(self):
        return f"GasState(T={self.T:.1f} K, p={self.p:.3e} Pa, M={self.M:.3f}, h={self.h:.4e} J/kg)"


class Aerothermodynamics:
    def __init__(self, optimum, transport="equilibrium", cache=True):
        self.optimum = optimum
        self.__dict__.update(optimum)

        if transport not in ("equilibrium", "frozen"):
            raise ValueError("transport must be 'equilibrium' or 'frozen'")

        self.transport = transport
        cea.init()

        fuel = np.asarray(self.fu.fractions, dtype=float)
        oxidizer = np.asarray(self.ox.fractions, dtype=float)
        
        if fuel.sum() <= 0 or oxidizer.sum() <= 0:
            raise ValueError("Fuel and oxidizer fractions must have positive sums")
        fuel /= fuel.sum()
        oxidizer /= oxidizer.sum()

        self.reactant_names = [*self.fu.propellants, *self.ox.propellants]
        self.reactant_weights = np.concatenate((
            fuel / (1 + self.MR),
            oxidizer * self.MR / (1 + self.MR),
        ))
        self.reactant_temperatures = np.asarray(
            [self.T_fu_in] * len(self.fu.propellants)
            + [self.T_ox_in] * len(self.ox.propellants),
            dtype=float,
        )

        reactants = [
            cea.Reactant(name, temperature=float(T))
            for name, T in zip(self.reactant_names, self.reactant_temperatures)
        ]
        reactant_mix = cea.Mixture(reactants)
        product_mix = cea.Mixture(reactants, products_from_reactants=True)
        self._reactant_enthalpy = reactant_mix.calc_property(
            cea.ENTHALPY, self.reactant_weights, self.reactant_temperatures
        ) / cea.R

        self._eq_solver = cea.EqSolver(product_mix, reactants=reactant_mix, transport=True)
        self._eq_solution = cea.EqSolution(self._eq_solver)
        self._rocket_solver = cea.RocketSolver(product_mix, reactants=reactant_mix,
                                               transport=True)
        self._rocket_solution = cea.RocketSolution(self._rocket_solver)
        self.contour = None
        self._cache_enabled = cache
        self._cache = {}

    @classmethod
    def _from_design(cls, fu, ox, MR, p_c, *, eps=None, p_e=None,
                     F=None, mdot=None, A_t=None,
                     L_star=None, V_c=None, t_stay=None,
                     T_fu_in=298.15, T_ox_in=298.15,
                     p_amb=1.013e5, npts=15, **kw):
        """The single shared path used by every sizing constructor."""
        if (eps is None) == (p_e is None):
            raise ValueError("Provide exactly one of eps or p_e")
        if sum(v is not None for v in (F, mdot, A_t)) != 1:
            raise ValueError("Provide exactly one of F, mdot or A_t")
        if sum(v is not None for v in (L_star, V_c, t_stay)) != 1:
            raise ValueError("Provide exactly one of L_star, V_c or t_stay")

        cea.init()
        fuel = np.asarray(fu.fractions, dtype=float)
        oxidizer = np.asarray(ox.fractions, dtype=float)
        if fuel.sum() <= 0 or oxidizer.sum() <= 0:
            raise ValueError("Fuel and oxidizer fractions must have positive sums")
        fuel /= fuel.sum()
        oxidizer /= oxidizer.sum()

        names = list(fu.propellants) + list(ox.propellants)
        weights = np.concatenate((fuel / (1 + MR), oxidizer * MR / (1 + MR)))
        temperatures = np.asarray(
            [T_fu_in] * len(fu.propellants) + [T_ox_in] * len(ox.propellants)
        )
        reactants = [
            cea.Reactant(name, temperature=float(T))
            for name, T in zip(names, temperatures)
        ]
        reactant_mix = cea.Mixture(reactants)
        product_mix = cea.Mixture(reactants, products_from_reactants=True)
        solver = cea.RocketSolver(product_mix, reactants=reactant_mix, transport=True)
        solution = cea.RocketSolution(solver)
        reactant_h = reactant_mix.calc_property(
            cea.ENTHALPY, weights, temperatures
        ) / cea.R

        station = {"supar": float(eps)} if eps is not None else {"pi_p": p_c / p_e}
        solver.solve(solution, weights, p_c * PA_TO_BAR, hc=reactant_h, **station)
        if not solution.converged:
            raise RuntimeError("Design-point CEA rocket solve did not converge")

        exit_i = int(solution.num_pts) - 1
        c_star = float(np.atleast_1d(solution.c_star)[0])
        Isp_vac = float(np.atleast_1d(solution.Isp_vacuum)[exit_i]) / G0
        Isp_ideal_amb = float(np.atleast_1d(solution.Isp)[exit_i]) / G0
        rho_c = float(np.atleast_1d(solution.density)[0])
        T_c = float(np.atleast_1d(solution.T)[0])
        T_t = float(np.atleast_1d(solution.T)[1])
        p_t = float(np.atleast_1d(solution.P)[1])       # retained in bar
        eps = float(eps) if eps is not None else float(np.atleast_1d(solution.ae_at)[exit_i])

        if F is not None:
            mdot = F / (Isp_vac * G0)
        elif mdot is not None:
            F = mdot * Isp_vac * G0
        else:
            mdot = p_c * A_t / c_star
            F = mdot * Isp_vac * G0

        A_t = c_star * mdot / p_c
        if V_c is not None:
            L_star = V_c / A_t
        elif t_stay is not None:
            L_star = t_stay * mdot / (A_t * rho_c)

        mdot_fu = mdot / (1 + MR)
        mdot_ox = mdot - mdot_fu
        A_e = A_t * eps
        t_stay = L_star * A_t * rho_c / mdot
        V_c = mdot * t_stay / rho_c
        CF_vac = Isp_vac * G0 / c_star
        CF_amb = CF_vac - p_amb / p_c * eps
        CF_SL = CF_vac - SEA_LEVEL_PRESSURE / p_c * eps

        optimum = dict(
            fu=fu, ox=ox, MR=MR, p_c=p_c, p_t=p_t, T_c=T_c, T_t=T_t,
            F=F, eps=eps, L_star=L_star, c_star=c_star, p_amb=p_amb,
            Isp_ideal_amb=Isp_ideal_amb, Isp_vac=Isp_vac,
            Isp_amb=CF_amb * c_star / G0, Isp_SL=CF_SL * c_star / G0,
            CF_vac=CF_vac, CF_amb=CF_amb, CF_SL=CF_SL,
            mdot=mdot, mdot_fu=mdot_fu, mdot_ox=mdot_ox,
            t_stay=t_stay, A_t=A_t, A_e=A_e,
            r_t=math.sqrt(A_t / math.pi), r_e=math.sqrt(A_e / math.pi), V_c=V_c,
            T_fu_in=T_fu_in, T_ox_in=T_ox_in, npts=npts,
        )
        return cls(optimum, **kw)

    @classmethod
    def from_F_eps_Lstar(cls, fu, ox, MR, p_c, F, eps, L_star,
                         T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15, **kw):
        return cls._from_design(fu, ox, MR, p_c, F=F, eps=eps, L_star=L_star,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, npts=npts, **kw)

    @classmethod
    def from_F_pe_Lstar(cls, fu, ox, MR, p_c, F, p_e, L_star,
                        T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15, **kw):
        return cls._from_design(fu, ox, MR, p_c, F=F, p_e=p_e, L_star=L_star,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, npts=npts, **kw)

    @classmethod
    def from_mdot_eps_Lstar(cls, fu, ox, MR, p_c, mdot, eps, L_star,
                            T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15, **kw):
        return cls._from_design(fu, ox, MR, p_c, mdot=mdot, eps=eps, L_star=L_star,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, npts=npts, **kw)

    @classmethod
    def from_mdot_pe_Lstar(cls, fu, ox, MR, p_c, mdot, p_e, L_star,
                           T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15, **kw):
        return cls._from_design(fu, ox, MR, p_c, mdot=mdot, p_e=p_e, L_star=L_star,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, npts=npts, **kw)

    @classmethod
    def from_At_eps_Lstar(cls, fu, ox, MR, p_c, A_t, eps, L_star,
                          T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15, **kw):
        return cls._from_design(fu, ox, MR, p_c, A_t=A_t, eps=eps, L_star=L_star,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, npts=npts, **kw)

    @classmethod
    def from_rt_eps_Lstar(cls, fu, ox, MR, p_c, r_t, eps, L_star, **kw):
        return cls.from_At_eps_Lstar(fu, ox, MR, p_c, math.pi * r_t**2,
                                     eps, L_star, **kw)

    @classmethod
    def from_F_eps_Vc(cls, fu, ox, MR, p_c, F, eps, V_c,
                      T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15, **kw):
        return cls._from_design(fu, ox, MR, p_c, F=F, eps=eps, V_c=V_c,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, npts=npts, **kw)

    @classmethod
    def from_F_eps_tstay(cls, fu, ox, MR, p_c, F, eps, t_stay,
                         T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15, **kw):
        return cls._from_design(fu, ox, MR, p_c, F=F, eps=eps, t_stay=t_stay,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, npts=npts, **kw)

    @classmethod
    def from_mass_flows_eps_Lstar(cls, fu, ox, p_c, mdot_fu, mdot_ox,
                                  eps, L_star, **kw):
        if mdot_fu <= 0:
            raise ValueError("mdot_fu must be positive")
        return cls.from_mdot_eps_Lstar(fu, ox, mdot_ox / mdot_fu, p_c,
                                       mdot_fu + mdot_ox, eps, L_star, **kw)

    @classmethod
    def from_optimum(cls, optimum, **kw):
        return cls(optimum, **kw)

    def compute_aerothermodynamics(self, contour, Nt=64):
        if Nt != 64:
            warnings.warn("Nt is ignored; station properties are solved live.",
                          DeprecationWarning, stacklevel=2)
        self.contour = contour
        self.x_nodes = np.linspace(contour.xs[0], contour.xs[-1], self.npts)
        self.clear_cache()
        return self

    attach_contour = compute_aerothermodynamics

    def equilibrium_state(self, problem, value, pressure):
        """Solve TP, HP or SP using base-SI inputs."""
        problem = problem.lower()
        if problem == "tp":
            cea_problem, cea_value = cea.TP, value
        elif problem == "hp":
            cea_problem, cea_value = cea.HP, value / cea.R
        elif problem == "sp":
            cea_problem, cea_value = cea.SP, value / cea.R
        else:
            raise ValueError(f"Unsupported problem type {problem!r}")

        self._eq_solver.solve(self._eq_solution, cea_problem, cea_value,
                              pressure * PA_TO_BAR, self.reactant_weights)
        sol = self._eq_solution
        if not sol.converged:
            raise RuntimeError(f"CEA {problem.upper()} solve did not converge")

        equilibrium = self.transport == "equilibrium"
        p = float(sol.P) * BAR_TO_PA
        return GasState(
            T=float(sol.T), p=p, rho=float(sol.density),
            cp=float(sol.cp_eq if equilibrium else sol.cp_fr) * KJ_TO_J,
            gamma=float(sol.gamma_s), h=float(sol.enthalpy) * KJ_TO_J,
            a=math.sqrt(float(sol.gamma_s) * p / sol.density),
            mu=float(sol.viscosity) * MILLIPOISE_TO_PA_S,
            k=float(sol.conductivity_eq if equilibrium else sol.conductivity_fr)
              * MW_PER_CM_K_TO_W_PER_M_K,
            Pr=float(sol.Pr_eq if equilibrium else sol.Pr_fr),
            mw=float(sol.MW), X=dict(sol.mole_fractions),
        )

    _eq_state = equilibrium_state

    def get_state(self, x, T=None, h=None):
        """Return all properties at x, optionally at an imposed T or h."""
        if T is not None and h is not None:
            raise ValueError("Provide only one of T or h")

        key = (float(x), None if T is None else float(T),
               None if h is None else float(h), self.transport)
        if self._cache_enabled and key in self._cache:
            return self._cache[key]

        if T is not None or h is not None:
            state = self.equilibrium_state(
                "tp" if T is not None else "hp",
                T if T is not None else h,
                self.get_state(x).p,
            )
        else:
            if self.contour is None:
                raise RuntimeError("Attach a contour before querying station properties")

            area_ratio = self.contour.A(x) / float(self.contour.A_t)
            if area_ratio <= 1 + 1e-10:
                station = {}
            elif x < 0:
                station = {"subar": float(area_ratio)}
            else:
                station = {"supar": float(area_ratio)}

            self._rocket_solver.solve(
                self._rocket_solution, self.reactant_weights,
                self.p_c * PA_TO_BAR, hc=self._reactant_enthalpy, **station
            )
            sol = self._rocket_solution
            if not sol.converged:
                raise RuntimeError(f"CEA rocket solve failed at x={x:.6g} m")

            equilibrium = self.transport == "equilibrium"
            cp = "cp_eq" if equilibrium else "cp_fr"
            conductivity = "conductivity_eq" if equilibrium else "conductivity_fr"
            Pr = "Pr_eq" if equilibrium else "Pr_fr"
            state = GasState(
                T=float(np.atleast_1d(sol.T)[-1]),
                p=float(np.atleast_1d(sol.P)[-1]) * BAR_TO_PA,
                M=float(np.atleast_1d(sol.Mach)[-1]),
                rho=float(np.atleast_1d(sol.density)[-1]),
                cp=float(np.atleast_1d(getattr(sol, cp))[-1]) * KJ_TO_J,
                gamma=float(np.atleast_1d(sol.gamma_s)[-1]),
                h=float(np.atleast_1d(sol.enthalpy)[-1]) * KJ_TO_J,
                a=float(np.atleast_1d(sol.sonic_velocity)[-1]),
                mu=float(np.atleast_1d(sol.viscosity)[-1]) * MILLIPOISE_TO_PA_S,
                k=float(np.atleast_1d(getattr(sol, conductivity))[-1])
                  * MW_PER_CM_K_TO_W_PER_M_K,
                Pr=float(np.atleast_1d(getattr(sol, Pr))[-1]),
                mw=float(np.atleast_1d(sol.MW)[-1]),
                X={name: float(np.atleast_1d(values)[-1])
                   for name, values in sol.mole_fractions.items()},
            )

        if self._cache_enabled:
            self._cache[key] = state
        return state

    def clear_cache(self): self._cache.clear()

    @property
    def cache_size(self): return len(self._cache)

    def get_T(self, x, T=None, h=None): return self.get_state(x, T, h).T
    def get_p(self, x, T=None, h=None): return self.get_state(x).p
    def get_M(self, x, T=None, h=None): return self.get_state(x, T, h).M
    def get_rho(self, x, T=None, h=None): return self.get_state(x, T, h).rho
    def get_cp(self, x, T=None, h=None): return self.get_state(x, T, h).cp
    def get_gamma(self, x, T=None, h=None): return self.get_state(x, T, h).gamma
    def get_h(self, x, T=None, h=None): return self.get_state(x, T, h).h
    def get_H(self, x, T=None, h=None): return self.get_state(x, T, h).h
    def get_a(self, x, T=None, h=None): return self.get_state(x, T, h).a
    def get_mu(self, x, T=None, h=None): return self.get_state(x, T, h).mu
    def get_k(self, x, T=None, h=None): return self.get_state(x, T, h).k
    def get_Pr(self, x, T=None, h=None): return self.get_state(x, T, h).Pr
    def get_molecular_weight(self, x, T=None, h=None): return self.get_state(x, T, h).mw
    def get_X(self, x, T=None, h=None): return self.get_state(x, T, h).X


def verify_state_scaling(aero, T=3000.0, P=2e6, rtol=1e-4):
    """Check TP -> HP and, if available, TP -> SP round trips."""
    reference = aero.equilibrium_state("tp", T, P)
    checks = {"hp": reference.h}
    try:
        checks["sp"] = float(aero._eq_solution.entropy) * KJ_TO_J
    except Exception:
        pass

    results = {}
    for problem, value in checks.items():
        try:
            state = aero.equilibrium_state(problem, value, reference.p)
            results[problem] = {
                "T_ref": reference.T,
                "T_got": state.T,
                "ok": math.isclose(state.T, reference.T, rel_tol=rtol),
            }
        except Exception as exc:
            results[problem] = {"T_ref": reference.T, "error": repr(exc), "ok": False}

    failed = [problem for problem, result in results.items() if not result["ok"]]
    if failed:
        raise AssertionError(f"State normalization failed for {failed}: {results}")
    return results
