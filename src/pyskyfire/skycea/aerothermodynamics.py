"""CEA-backed gas properties along a rocket-engine contour.

All public values use base SI units. The old constructor and getter API is
retained, but the implementation has only one design path and one state path.
"""

from __future__ import annotations

import math
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
# Default lower boundary for direct temperature-conditioned CEA solves [K].
# Below this configurable boundary, properties are continued along their
# right-hand C1 tangent instead of being clamped to a constant value.
MIN_CEA_TEMPERATURE = 200.0

_LOG_TANGENT_FIELDS = (
    "rho",
    "cp",
    "gamma",
    "a",
    "mu",
    "k",
    "Pr",
    "mw",
)


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
    def __init__(
        self,
        optimum,
        transport="equilibrium",
        cache=True,
        minimum_cea_temperature=MIN_CEA_TEMPERATURE,
    ):
        self.optimum = optimum
        self.__dict__.update(optimum)

        if transport not in ("equilibrium", "frozen"):
            raise ValueError("transport must be 'equilibrium' or 'frozen'")

        self.transport = transport
        self._temperature_tangent_cache = {}
        self.minimum_cea_temperature = minimum_cea_temperature
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

        reactant_mix, product_mix = self._make_mixtures()
        self._reactant_enthalpy = reactant_mix.calc_property(
            cea.ENTHALPY, self.reactant_weights, self.reactant_temperatures
        ) / cea.R

        self._make_solvers(reactant_mix, product_mix)
        self.contour = None
        self._cache_enabled = cache
        self._cache = {}

    @property
    def minimum_cea_temperature(self):
        """Lowest temperature sent directly to CEA [K]."""
        return self._minimum_cea_temperature

    @minimum_cea_temperature.setter
    def minimum_cea_temperature(self, value):
        value = float(value)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(
                "minimum_cea_temperature must be finite and positive"
            )
        self._minimum_cea_temperature = value
        tangent_cache = getattr(self, "_temperature_tangent_cache", None)
        if tangent_cache is not None:
            tangent_cache.clear()

    def _make_mixtures(self):
        """Build the CEA mixtures represented by this object's Python state."""
        reactants = [
            cea.Reactant(name, temperature=float(T))
            for name, T in zip(self.reactant_names, self.reactant_temperatures)
        ]
        reactant_mix = cea.Mixture(reactants)
        product_mix = cea.Mixture(reactants, products_from_reactants=True)
        return reactant_mix, product_mix

    def _make_solvers(self, reactant_mix=None, product_mix=None):
        """Create the native CEA handles, which cannot themselves be pickled."""
        cea.init()
        if reactant_mix is None or product_mix is None:
            reactant_mix, product_mix = self._make_mixtures()

        self._eq_solver = cea.EqSolver(product_mix, reactants=reactant_mix, transport=True)
        self._eq_solution = cea.EqSolution(self._eq_solver)
        self._rocket_solver = cea.RocketSolver(product_mix, reactants=reactant_mix, transport=True)
        self._rocket_solution = cea.RocketSolution(self._rocket_solver)

    def __getstate__(self):
        """Return persistent state without native handles or runtime cache."""
        state = self.__dict__.copy()
        for name in (
            "_eq_solver",
            "_eq_solution",
            "_rocket_solver",
            "_rocket_solution",
            "_cache",
            "_temperature_tangent_cache",
        ):
            state.pop(name, None)
        return state

    def __setstate__(self, state):
        """Restore Python state and recreate the process-local CEA handles."""
        self.__dict__.update(state)
        if not hasattr(self, "_minimum_cea_temperature"):
            self._minimum_cea_temperature = MIN_CEA_TEMPERATURE
        self._cache = {}
        self._temperature_tangent_cache = {}
        self._make_solvers()

    @classmethod
    def _from_design(cls, fu, ox, MR, p_c, *, eps=None, p_e=None,
                     F=None, mdot=None, A_t=None,
                     L_star=None, V_c=None, t_stay=None,
                     T_fu_in=298.15, T_ox_in=298.15,
                     p_amb=1.013e5, **kw):
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
        # CEA's Isp is the exit velocity, i.e. the impulse of a nozzle expanded
        # exactly to its own exit pressure. It carries no pressure-thrust term
        # and does not depend on p_amb.
        Isp_optimum = float(np.atleast_1d(solution.Isp)[exit_i]) / G0
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
            Isp_optimum=Isp_optimum, Isp_vac=Isp_vac,
            Isp_amb=CF_amb * c_star / G0, Isp_SL=CF_SL * c_star / G0,
            CF_vac=CF_vac, CF_amb=CF_amb, CF_SL=CF_SL,
            mdot=mdot, mdot_fu=mdot_fu, mdot_ox=mdot_ox,
            t_stay=t_stay, A_t=A_t, A_e=A_e,
            r_t=math.sqrt(A_t / math.pi), r_e=math.sqrt(A_e / math.pi), V_c=V_c,
            T_fu_in=T_fu_in, T_ox_in=T_ox_in,
        )
        return cls(optimum, **kw)

    @classmethod
    def from_F_eps_Lstar(cls, fu, ox, MR, p_c, F, eps, L_star,
                         T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, **kw):
        return cls._from_design(fu, ox, MR, p_c, F=F, eps=eps, L_star=L_star,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, **kw)

    @classmethod
    def from_F_pe_Lstar(cls, fu, ox, MR, p_c, F, p_e, L_star,
                        T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, **kw):
        return cls._from_design(fu, ox, MR, p_c, F=F, p_e=p_e, L_star=L_star,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, **kw)

    @classmethod
    def from_mdot_eps_Lstar(cls, fu, ox, MR, p_c, mdot, eps, L_star,
                            T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, **kw):
        return cls._from_design(fu, ox, MR, p_c, mdot=mdot, eps=eps, L_star=L_star,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, **kw)

    @classmethod
    def from_mdot_pe_Lstar(cls, fu, ox, MR, p_c, mdot, p_e, L_star,
                           T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, **kw):
        return cls._from_design(fu, ox, MR, p_c, mdot=mdot, p_e=p_e, L_star=L_star,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, **kw)

    @classmethod
    def from_At_eps_Lstar(cls, fu, ox, MR, p_c, A_t, eps, L_star,
                          T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, **kw):
        return cls._from_design(fu, ox, MR, p_c, A_t=A_t, eps=eps, L_star=L_star,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, **kw)

    @classmethod
    def from_rt_eps_Lstar(cls, fu, ox, MR, p_c, r_t, eps, L_star, **kw):
        return cls.from_At_eps_Lstar(fu, ox, MR, p_c, math.pi * r_t**2,
                                     eps, L_star, **kw)

    @classmethod
    def from_F_eps_Vc(cls, fu, ox, MR, p_c, F, eps, V_c,
                      T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, **kw):
        return cls._from_design(fu, ox, MR, p_c, F=F, eps=eps, V_c=V_c,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, **kw)

    @classmethod
    def from_F_eps_tstay(cls, fu, ox, MR, p_c, F, eps, t_stay,
                         T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, **kw):
        return cls._from_design(fu, ox, MR, p_c, F=F, eps=eps, t_stay=t_stay,
                                T_fu_in=T_fu_in, T_ox_in=T_ox_in,
                                p_amb=p_amb, **kw)

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

    def compute_aerothermodynamics(self, contour):
        """Attach chamber geometry for live, coordinate-driven CEA queries."""
        self.contour = contour
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

    @staticmethod
    def _forward_derivative(values, step):
        """Second-order right-hand derivative at the first sample."""
        value_0, value_1, value_2 = values
        return (-3.0 * value_0 + 4.0 * value_1 - value_2) / (2.0 * step)

    def _temperature_tangent(self, x, pressure):
        """Build or retrieve the C1 continuation at the configured boundary."""
        boundary = self.minimum_cea_temperature
        key = (float(x), float(pressure), boundary, self.transport)
        cached = self._temperature_tangent_cache.get(key)
        if cached is not None:
            return cached

        step = max(0.1, 1.0e-3 * boundary)
        states = tuple(
            self.equilibrium_state("tp", boundary + offset * step, pressure)
            for offset in range(3)
        )
        base = states[0]
        linear_slopes = {
            "h": self._forward_derivative(
                tuple(state.h for state in states), step
            )
        }
        log_slopes = {}
        for name in _LOG_TANGENT_FIELDS:
            samples = np.asarray(
                [getattr(state, name) for state in states], dtype=float
            )
            if np.all(np.isfinite(samples)) and np.all(samples > 0.0):
                log_slopes[name] = self._forward_derivative(
                    tuple(np.log(samples)), step
                )
            else:
                linear_slopes[name] = self._forward_derivative(
                    tuple(samples), step
                )

        tangent = (base, linear_slopes, log_slopes)
        self._temperature_tangent_cache[key] = tangent
        return tangent

    def _extrapolate_temperature_state(self, x, temperature, pressure):
        """Continue a temperature-conditioned CEA state below its boundary."""
        base, linear_slopes, log_slopes = self._temperature_tangent(x, pressure)
        delta_temperature = float(temperature) - self.minimum_cea_temperature
        values = {
            name: getattr(base, name)
            for name in GasState.__dataclass_fields__
            if name != "X"
        }
        values["T"] = float(temperature)
        values["p"] = float(pressure)
        for name, slope in linear_slopes.items():
            values[name] = getattr(base, name) + slope * delta_temperature
        for name, slope in log_slopes.items():
            exponent = np.clip(slope * delta_temperature, -700.0, 700.0)
            values[name] = getattr(base, name) * math.exp(float(exponent))
        values["X"] = dict(base.X)
        return GasState(**values)

    def _extrapolate_enthalpy_state(self, x, enthalpy, pressure):
        """Invert the same tangent when an HP query falls below its boundary."""
        base, linear_slopes, _ = self._temperature_tangent(x, pressure)
        enthalpy_slope = linear_slopes["h"]
        if not math.isfinite(enthalpy_slope) or enthalpy_slope <= 0.0:
            raise RuntimeError(
                "CEA temperature-boundary enthalpy tangent is not invertible"
            )
        temperature = self.minimum_cea_temperature + (
            float(enthalpy) - base.h
        ) / enthalpy_slope
        if temperature >= self.minimum_cea_temperature:
            raise RuntimeError(
                "CEA HP solve failed for a state above the configured "
                "temperature boundary"
            )
        state = self._extrapolate_temperature_state(x, temperature, pressure)
        state.h = float(enthalpy)
        return state

    def get_state(self, x, T=None, h=None):
        """Return all properties at x, optionally at an imposed T or h."""
        if T is not None and h is not None:
            raise ValueError("Provide only one of T or h")

        # Only the contour state is stable enough to reuse. Temperature- and
        # enthalpy-conditioned states are Newton trial points that almost never
        # recur and made this cache grow without bound during coupled solves.
        cacheable = self._cache_enabled and T is None and h is None
        key = (float(x), None, None, self.transport)
        if cacheable and key in self._cache:
            return self._cache[key]

        if T is not None or h is not None:
            pressure = self.get_state(x).p
            if T is not None and float(T) < self.minimum_cea_temperature:
                state = self._extrapolate_temperature_state(
                    x, float(T), pressure
                )
            elif T is not None:
                state = self.equilibrium_state(
                    "tp",
                    float(T),
                    pressure,
                )
            else:
                try:
                    state = self.equilibrium_state("hp", h, pressure)
                except RuntimeError:
                    state = self._extrapolate_enthalpy_state(
                        x, h, pressure
                    )
                else:
                    if state.T < self.minimum_cea_temperature:
                        state = self._extrapolate_enthalpy_state(
                            x, h, pressure
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

        if cacheable:
            self._cache[key] = state
        return state

    def clear_cache(self):
        self._cache.clear()
        tangent_cache = getattr(self, "_temperature_tangent_cache", None)
        if tangent_cache is not None:
            tangent_cache.clear()

    def reset_solver_state(self):
        """Recreate native CEA handles and clear all cached station states.

        NASA CEA's nonlinear solution handles retain their previous solution
        as an initial state. That is useful during a single axial march, but
        separate simulations of the same inputs must not depend on which
        simulation ran immediately before them. Resetting the native handles
        provides a deterministic boundary between such simulations without
        rebuilding this object's design configuration.
        """
        self._make_solvers()
        self.clear_cache()

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
    def get_v(self, x, T=None, h=None):
        state = self.get_state(x, T, h)
        return state.M * state.a
    def get_T0(self, x, T=None, h=None):
        """Total temperature [K], from conserved total enthalpy.

        Inverted through the equilibrium equation of state rather than taken
        from the frozen-gamma isentropic relation: gamma is a local property
        here, so that relation drifts above the chamber value down the nozzle.

        The inversion is done at the *local static* pressure, the same
        convention the boundary-layer reference and recovery states use, which
        keeps this comparable with ``T_aw``. It is therefore not the classical
        isentropic stagnation temperature, and it is not constant along an
        equilibrium nozzle: total enthalpy is conserved, but as the pressure
        falls the equilibrium shifts towards dissociation, so the same enthalpy
        corresponds to a lower temperature. It never exceeds its chamber value.
        """
        state = self.get_state(x, T, h)
        H0 = state.h + 0.5 * (state.M * state.a) ** 2
        try:
            return self.get_state(x, h=H0).T
        except Exception:  # noqa: BLE001
            return state.T * (1.0 + 0.5 * (state.gamma - 1.0) * state.M ** 2)
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
