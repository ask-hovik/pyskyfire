
"""
film_solver.py

Contour-based Grisson-style film cooling model for pyskyfire.

Key design choices in this version
----------------------------------
- Geometry is read directly from thrust_chamber.contour.
- Film-cooling inputs live on thrust_chamber.film_cooling.
- ThrustChamber is expected to resolve film_cooling.x_fraction -> film_cooling.x
  before this solver is used.
- Coolant liquid/vapor properties are pulled from CoolProp through the existing
  coolant transport object on the selected cooling circuit.
- H2O and CO2 mole fractions used by the current Grisson gas-radiation model
  are taken directly from FilmCooling, not from combustion_transport.

Important assumptions
---------------------
- The injected coolant is a pure fluid.
- Gas-side mean molecular weight is still expected to be available from
  thrust_chamber.combustion_transport.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence, Any

from scipy.optimize import brentq

try:
    from CoolProp.CoolProp import PropsSI
except Exception:  # pragma: no cover
    PropsSI = None


@dataclass
class GasProperties:
    T_g: float
    G_ch: float
    mu_g: float
    Cp_g: float
    Pr: float
    M_g: float
    gamma: float
    mole_fraction_H2O: float
    mole_fraction_CO2: float
    chamber_pressure: float
    turbulence_intensity: float = 0.0


@dataclass
class LiquidCoolantProperties:
    T_injection: float
    T_saturation: float
    latent_heat: float
    Cp_liquid: float
    mu_liquid: float
    mu_vapor: float
    rho_liquid: float
    rho_vapor: float
    surface_tension: float
    absorptivity: float
    M_c: float
    Cp_vapor: float


@dataclass
class LiquidFilmResults:
    x: list[float] = field(default_factory=list)
    Gamma: list[float] = field(default_factory=list)
    m_vap: list[float] = field(default_factory=list)
    h_conv: list[float] = field(default_factory=list)
    Q_conv: list[float] = field(default_factory=list)
    Q_rad: list[float] = field(default_factory=list)
    T_film: list[float] = field(default_factory=list)
    x_dryout: Optional[float] = None
    Gamma_dryout: float = 0.0


@dataclass
class GaseousFilmResults:
    x: list[float] = field(default_factory=list)
    Mbl: list[float] = field(default_factory=list)
    T_aw: list[float] = field(default_factory=list)
    T_w: list[float] = field(default_factory=list)
    h_conv: list[float] = field(default_factory=list)
    Q_rad: list[float] = field(default_factory=list)

class ContourGeometryEvaluator:
    """Thin adapter exposing the geometry that the film model actually needs."""

    def __init__(self, contour):
        self.contour = contour

    def diameter(self, x: float) -> float:
        return 2.0 * float(self.contour.r(x))

    def area(self, x: float) -> float:
        return float(self.contour.A(x))

    def dD_dx(self, x: float) -> float:
        return 2.0 * float(self.contour.dr_dx(x))

    @property
    def D_chamber(self) -> float:
        return 2.0 * float(self.contour.r_c)

    @property
    def D_throat(self) -> float:
        return 2.0 * float(self.contour.r_t)

    @property
    def A_chamber(self) -> float:
        return float(self.contour.A_c)

    @property
    def A_throat(self) -> float:
        return float(self.contour.A_t)


def _require(value: Any, name: str) -> Any:
    if value is None:
        raise AttributeError(f"Required quantity '{name}' was not available.")
    return value


def _try_call(obj: Any, method_name: str, *args, default=None):
    fn = getattr(obj, method_name, None)
    if fn is None:
        return default
    return fn(*args)


def _extract_gas_molecular_weight(combustion_transport, x: float) -> float:
    candidates = [
        lambda: _try_call(combustion_transport, "get_MW", x),
        lambda: _try_call(combustion_transport, "get_molecular_weight", x),
        lambda: _try_call(combustion_transport, "get_molar_mass", x),
        lambda: getattr(combustion_transport, "MW", None),
        lambda: getattr(combustion_transport, "molecular_weight", None),
    ]
    for getter in candidates:
        try:
            value = getter()
        except TypeError:
            value = None
        if value is not None:
            return float(value)
    raise AttributeError(
        "Could not extract gas molecular weight from combustion_transport. "
        "Add a clean accessor such as get_molecular_weight(x)."
    )


def _coolprop_fluid_name(coolant_transport) -> str:
    for attr_name in ("fluid_name", "coolprop_name", "fluid", "_fluid_name", "_coolprop_name"):
        value = getattr(coolant_transport, attr_name, None)
        if isinstance(value, str) and value:
            return value

    fluid_obj = getattr(coolant_transport, "fluid", None)
    if isinstance(fluid_obj, str) and fluid_obj:
        return fluid_obj

    propellants = getattr(fluid_obj, "propellants", None)
    fractions = getattr(fluid_obj, "fractions", None)
    if propellants and len(propellants) == 1 and fractions and len(fractions) == 1 and abs(float(fractions[0]) - 1.0) < 1e-12:
        return str(propellants[0])

    raise AttributeError(
        "Could not infer a CoolProp fluid name from coolant_transport. "
        "Expose a string like coolant_transport.fluid_name = 'ethanol'."
    )


def _coolprop_scalar(output: str, name1: str, value1: float, name2: str, value2: float, fluid: str) -> float:
    if PropsSI is None:
        raise ImportError("CoolProp is required by film_solver.py but could not be imported.")
    return float(PropsSI(output, name1, float(value1), name2, float(value2), fluid))


def _build_gas_properties(thrust_chamber, film_cooling) -> GasProperties:
    contour = thrust_chamber.contour
    gas = thrust_chamber.combustion_transport

    x_ref = float(film_cooling.x)
    T_g = float(gas.get_T(x_ref))
    mu_g = float(gas.get_mu(x_ref))
    Cp_g = float(gas.get_cp(x_ref))
    Pr = float(gas.get_Pr(x_ref))
    gamma = float(gas.get_gamma(x_ref))
    chamber_pressure = float(gas.get_p(x_ref))
    mdot_g = float(gas.mdot)
    G_ch = mdot_g / float(contour.A(x_ref))

    M_g = _extract_gas_molecular_weight(gas, x_ref)

    return GasProperties(
        T_g=T_g,
        G_ch=G_ch,
        mu_g=mu_g,
        Cp_g=Cp_g,
        Pr=Pr,
        M_g=M_g,
        gamma=gamma,
        mole_fraction_H2O=float(film_cooling.mole_fraction_H2O),
        mole_fraction_CO2=float(film_cooling.mole_fraction_CO2),
        chamber_pressure=chamber_pressure,
        turbulence_intensity=float(film_cooling.turbulence_intensity),
    )


def _build_coolant_properties(thrust_chamber, boundary_conditions, circuit_index: int) -> LiquidCoolantProperties:
    coolant_transport = thrust_chamber.cooling_circuits[circuit_index].coolant_transport
    fluid = _coolprop_fluid_name(coolant_transport)

    T_inj = float(boundary_conditions.T_coolant_in)
    p_inj = float(boundary_conditions.p_coolant_in)

    T_sat = _coolprop_scalar("T", "P", p_inj, "Q", 0.0, fluid)
    h_l_sat = _coolprop_scalar("Hmass", "P", p_inj, "Q", 0.0, fluid)
    h_v_sat = _coolprop_scalar("Hmass", "P", p_inj, "Q", 1.0, fluid)
    latent_heat = h_v_sat - h_l_sat

    Cp_liquid = _coolprop_scalar("Cpmass", "T", T_inj, "P", p_inj, fluid)
    mu_liquid = _coolprop_scalar("VISCOSITY", "T", T_inj, "P", p_inj, fluid)
    rho_liquid = _coolprop_scalar("Dmass", "T", T_inj, "P", p_inj, fluid)

    mu_vapor = _coolprop_scalar("VISCOSITY", "P", p_inj, "Q", 1.0, fluid)
    rho_vapor = _coolprop_scalar("Dmass", "P", p_inj, "Q", 1.0, fluid)
    surface_tension = _coolprop_scalar("SURFACE_TENSION", "P", p_inj, "Q", 0.0, fluid)
    Cp_vapor = _coolprop_scalar("Cpmass", "P", p_inj, "Q", 1.0, fluid)
    M_c = 1000.0 * _coolprop_scalar("MOLAR_MASS", "P", p_inj, "Q", 1.0, fluid)

    return LiquidCoolantProperties(
        T_injection=T_inj,
        T_saturation=T_sat,
        latent_heat=latent_heat,
        Cp_liquid=Cp_liquid,
        mu_liquid=mu_liquid,
        mu_vapor=mu_vapor,
        rho_liquid=rho_liquid,
        rho_vapor=rho_vapor,
        surface_tension=surface_tension,
        absorptivity=0.0,
        M_c=M_c,
        Cp_vapor=Cp_vapor,
    )


class GasEmittanceCalculator:
    _H2O_TABLE = [
        [500,  0.600, 0.010, 1.60],
        [1000, 0.500, 0.025, 1.45],
        [1500, 0.430, 0.048, 1.35],
        [2000, 0.370, 0.080, 1.25],
        [2500, 0.325, 0.120, 1.18],
        [3000, 0.290, 0.175, 1.12],
    ]

    _CO2_TABLE = [
        [500,  0.290, 0.0020, 1.40],
        [1000, 0.220, 0.0060, 1.30],
        [1500, 0.175, 0.0140, 1.22],
        [2000, 0.145, 0.0280, 1.16],
        [2500, 0.125, 0.0500, 1.12],
        [3000, 0.110, 0.0850, 1.08],
    ]

    def __init__(self, gas: GasProperties):
        self.gas = gas

    def _parabolic_interpolate(self, table: list[list[float]], T: float) -> tuple[float, float, float]:
        temps = [row[0] for row in table]
        n_pts = len(temps)

        if T <= temps[0]:
            return table[0][1], table[0][2], table[0][3]
        if T >= temps[-1]:
            return table[-1][1], table[-1][2], table[-1][3]

        idx = 1
        for i in range(1, n_pts - 1):
            if T <= temps[i]:
                idx = i
                break
        else:
            idx = n_pts - 2

        i0, i1, i2 = idx - 1, idx, idx + 1
        T0, T1, T2 = temps[i0], temps[i1], temps[i2]

        def quad_interp(v0, v1, v2):
            L0 = (T - T1) * (T - T2) / ((T0 - T1) * (T0 - T2))
            L1 = (T - T0) * (T - T2) / ((T1 - T0) * (T1 - T2))
            L2 = (T - T0) * (T - T1) / ((T2 - T0) * (T2 - T1))
            return v0 * L0 + v1 * L1 + v2 * L2

        eps_f = quad_interp(table[i0][1], table[i1][1], table[i2][1])
        c = quad_interp(table[i0][2], table[i1][2], table[i2][2])
        n = quad_interp(table[i0][3], table[i1][3], table[i2][3])
        return eps_f, c, n

    def _species_emittance(self, table: list[list[float]], T: float, rho_opt: float) -> float:
        eps_f, c, n = self._parabolic_interpolate(table, T)
        if rho_opt <= 0.0:
            return 0.0
        return eps_f / (1.0 + (rho_opt / c) ** (-n)) ** (1.0 / n)

    def gas_emittance(self, T: float, D: float) -> float:
        L_ft = D * 3.2808
        p_atm = self.gas.chamber_pressure * 9.8692e-6

        rho_H2O = self.gas.mole_fraction_H2O * p_atm * L_ft
        rho_CO2 = self.gas.mole_fraction_CO2 * p_atm * L_ft

        eps_H2O = self._species_emittance(self._H2O_TABLE, T, rho_H2O)
        eps_CO2 = self._species_emittance(self._CO2_TABLE, T, rho_CO2)

        rho_sum = rho_H2O + rho_CO2
        eps_overlap = 0.5 * (
            self._species_emittance(self._H2O_TABLE, T, rho_sum)
            + self._species_emittance(self._CO2_TABLE, T, rho_sum)
        )

        eps_total = eps_H2O + eps_CO2 - eps_overlap
        return max(0.0, min(1.0, eps_total))


class LiquidFilmSolver:
    STEFAN_BOLTZMANN = 5.6704e-8

    def __init__(self, thrust_chamber, boundary_conditions, circuit_index: int = 0):
        self.thrust_chamber = thrust_chamber
        self.boundary_conditions = boundary_conditions
        self.circuit_index = int(circuit_index)

        self.film = _require(getattr(thrust_chamber, "film_cooling", None), "thrust_chamber.film_cooling")
        if self.film.x is None:
            raise ValueError(
                "thrust_chamber.film_cooling.x is None. Resolve x_fraction -> x "
                "inside ThrustChamber before constructing the film solver."
            )

        self.geom_eval = ContourGeometryEvaluator(thrust_chamber.contour)
        self.gas = _build_gas_properties(thrust_chamber, self.film)
        self.coolant = _build_coolant_properties(thrust_chamber, boundary_conditions, self.circuit_index)
        self.coolant.absorptivity = float(self.film.liquid_absorptivity)
        self.emittance_calc = GasEmittanceCalculator(self.gas)

        self._L_eff = self.coolant.latent_heat + self.coolant.Cp_liquid * (self.coolant.T_saturation - self.coolant.T_injection)
        self._Gamma_0 = float(self.film.coolant_mass_flow_rate) / float(self.film.film_injection_perimeter)

        if self.coolant.M_c < self.gas.M_g:
            self._K_M_exp = 0.60
        else:
            self._K_M_exp = 0.35
        self._K_M = (self.gas.M_g / self.coolant.M_c) ** self._K_M_exp
        self._K_t = 1.0 + 4.0 * self.gas.turbulence_intensity

    def _local_gas_flux(self, x: float) -> float:
        A_ref = self.geom_eval.A_chamber
        A_x = self.geom_eval.area(x)
        return self.gas.G_ch * A_ref / A_x

    def _h0_colburn(self, x: float, x_eff: float) -> float:
        G = self._local_gas_flux(x)
        Re_x = G * x_eff / self.gas.mu_g
        St = 0.023 * Re_x ** (-0.2) * self.gas.Pr ** (-0.6)
        return St * G * self.gas.Cp_g

    def _transpiration_correction(self, h0: float, m_vap: float) -> float:
        if m_vap <= 0.0:
            return h0

        def residual(h):
            H = self.gas.Cp_g * self._K_M * m_vap / h
            return h * H / math.log(1.0 + H) - h0

        try:
            return brentq(residual, h0 * 1e-6, h0 * (1.0 + 1e-3), xtol=1e-6, maxiter=100)
        except ValueError:
            H_est = self.gas.Cp_g * self._K_M * m_vap / h0
            return h0 * math.log(1.0 + H_est) / H_est if H_est > 0 else h0

    def _radiative_heat_flux(self, x: float) -> float:
        D = self.geom_eval.diameter(x)
        eps_g = self.emittance_calc.gas_emittance(self.gas.T_g, D)
        return self.coolant.absorptivity * eps_g * self.STEFAN_BOLTZMANN * self.gas.T_g ** 4

    def _evaporation_rate(self, x: float, Gamma: float, x_inj: float) -> tuple[float, float, float, float]:
        T_film = self.coolant.T_saturation
        x_eff = x - x_inj + 0.9 * (self.gas.mu_g / self._local_gas_flux(x)) ** 1.25
        h0 = self._K_t * self._h0_colburn(x, max(x_eff, 1e-9))
        Q_rad = self._radiative_heat_flux(x)

        Q_conv_0 = h0 * (self.gas.T_g - T_film)
        m_vap_est = max(0.0, (Q_conv_0 + Q_rad) / self._L_eff)

        h_corrected = self._transpiration_correction(h0, m_vap_est)
        Q_conv = h_corrected * (self.gas.T_g - T_film)
        m_vap = max(0.0, (Q_conv + Q_rad) / self._L_eff)

        h_corrected = self._transpiration_correction(h0, m_vap)
        Q_conv = h_corrected * (self.gas.T_g - T_film)
        m_vap = max(0.0, (Q_conv + Q_rad) / self._L_eff)

        return m_vap, h_corrected, Q_conv, Q_rad

    def solve(self, x_array: Sequence[float]) -> LiquidFilmResults:
        results = LiquidFilmResults()
        x_inj = float(self.film.x)
        Gamma = self._Gamma_0

        start_idx = 0
        for i, x in enumerate(x_array):
            if x >= x_inj:
                start_idx = i
                break

        prev_x = float(x_array[start_idx])
        results.x.append(prev_x)
        results.Gamma.append(Gamma)

        m_vap, h, Q_conv, Q_rad = self._evaporation_rate(prev_x, Gamma, x_inj)
        results.m_vap.append(m_vap)
        results.h_conv.append(h)
        results.Q_conv.append(Q_conv)
        results.Q_rad.append(Q_rad)
        results.T_film.append(self.coolant.T_saturation)

        for i in range(start_idx + 1, len(x_array)):
            x = float(x_array[i])
            dx = x - prev_x
            Gamma_new = Gamma - m_vap * dx

            if Gamma_new <= 0.0:
                x_dryout = prev_x + Gamma / m_vap if m_vap > 0.0 else x
                results.x_dryout = x_dryout
                results.Gamma_dryout = 0.0
                break

            Gamma = Gamma_new
            m_vap, h, Q_conv, Q_rad = self._evaporation_rate(x, Gamma, x_inj)

            results.x.append(x)
            results.Gamma.append(Gamma)
            results.m_vap.append(m_vap)
            results.h_conv.append(h)
            results.Q_conv.append(Q_conv)
            results.Q_rad.append(Q_rad)
            results.T_film.append(self.coolant.T_saturation)

            prev_x = x
        else:
            results.x_dryout = None
            results.Gamma_dryout = Gamma

        return results


class GaseousFilmSolver:
    STEFAN_BOLTZMANN = 5.6704e-8

    def __init__(self, thrust_chamber, boundary_conditions, circuit_index: int = 0):
        self.thrust_chamber = thrust_chamber
        self.boundary_conditions = boundary_conditions
        self.circuit_index = int(circuit_index)

        self.film = _require(getattr(thrust_chamber, "film_cooling", None), "thrust_chamber.film_cooling")
        if self.film.x is None:
            raise ValueError(
                "thrust_chamber.film_cooling.x is None. Resolve x_fraction -> x "
                "inside ThrustChamber before constructing the film solver."
            )

        self.geom_eval = ContourGeometryEvaluator(thrust_chamber.contour)
        self.gas = _build_gas_properties(thrust_chamber, self.film)
        self.coolant = _build_coolant_properties(thrust_chamber, boundary_conditions, self.circuit_index)
        self.coolant.absorptivity = float(self.film.liquid_absorptivity)
        self.emittance_calc = GasEmittanceCalculator(self.gas)

        self._K_M = (self.coolant.M_c / self.gas.M_g) ** 0.14
        self._K_t = 1.0 + 10.2 * self.gas.turbulence_intensity

    def _local_gas_flux(self, x: float) -> float:
        A_ref = self.geom_eval.A_chamber
        A_x = self.geom_eval.area(x)
        return self.gas.G_ch * A_ref / A_x

    def _recovery_temperature(self, x: float) -> float:
        r = self.gas.Pr ** (1.0 / 3.0)
        T_0 = self.gas.T_g
        gamma = self.gas.gamma
        D = self.geom_eval.diameter(x)
        D_t = self.geom_eval.D_throat

        area_ratio = (D / D_t) ** 2
        T_s = T_0 * (2.0 / (gamma + 1.0)) if area_ratio <= 1.05 else T_0
        return T_0 - (1.0 - r) * (T_0 - T_s)

    def _convective_h_gaseous(self, x: float, Mbl: float) -> float:
        G = self._local_gas_flux(x)
        return 0.1963 * self._K_t * G * self.gas.Cp_g * (self.gas.mu_g / Mbl) ** 0.25

    def _radiative_heat_flux(self, x: float) -> float:
        D = self.geom_eval.diameter(x)
        eps_g = self.emittance_calc.gas_emittance(self.gas.T_g, D)
        return eps_g * self.STEFAN_BOLTZMANN * self.gas.T_g ** 4

    def _initial_conditions(self, x_i: float, Mc_total: float) -> tuple[float, float]:
        x_inj = float(self.film.x)
        xi = x_i - x_inj
        G_i = self._local_gas_flux(x_i)

        K = 0.1963 * self._K_t * G_i * (self.gas.mu_g / Mc_total) ** 0.25 / Mc_total
        Xi = K * xi

        Mbl_i = Mc_total * (1.0 + 0.325 * Xi ** 0.8) if Xi > 0.0 else Mc_total
        T_aw_i = self.coolant.T_saturation
        return Mbl_i, T_aw_i

    def _derivatives(self, x: float, Mbl: float, T_aw: float) -> tuple[float, float]:
        G = self._local_gas_flux(x)
        D = self.geom_eval.diameter(x)

        dMbl_dx_e = 0.1963 * self._K_t * G * (self.gas.mu_g / max(Mbl, 1e-10)) ** 0.25
        dD_dx = self.geom_eval.dD_dx(x)
        dMbl_dx_c = -Mbl * (1.0 / D) * dD_dx if D > 0.0 else 0.0
        dMbl_dx = dMbl_dx_e + dMbl_dx_c

        Mc = self.film.coolant_mass_flow_rate / (math.pi * D) if D > 0.0 else 0.0

        T_r = self._recovery_temperature(x)
        denom = Mbl + (self.coolant.Cp_vapor / (self._K_M * self.gas.Cp_g) - 1.0) * Mc
        denom = max(denom, 1e-10)

        dTaw_dx = dMbl_dx_e * (T_r - T_aw) / denom
        Q_rad = self._radiative_heat_flux(x)
        dTaw_dx += Q_rad / (self.gas.Cp_g * max(Mbl, 1e-10))

        return dMbl_dx, dTaw_dx

    def solve(self, x_array: Sequence[float], x_dryout: float, Mc_total: float) -> GaseousFilmResults:
        results = GaseousFilmResults()

        start_idx = 0
        for i, x in enumerate(x_array):
            if x >= x_dryout:
                start_idx = i
                break

        Mbl, T_aw = self._initial_conditions(x_dryout, Mc_total)

        prev_x = float(x_array[start_idx])
        h = self._convective_h_gaseous(prev_x, Mbl)
        Q_rad = self._radiative_heat_flux(prev_x)
        T_w = T_aw + Q_rad / h if h > 0.0 else T_aw

        results.x.append(prev_x)
        results.Mbl.append(Mbl)
        results.T_aw.append(T_aw)
        results.T_w.append(T_w)
        results.h_conv.append(h)
        results.Q_rad.append(Q_rad)

        for i in range(start_idx + 1, len(x_array)):
            x = float(x_array[i])
            dx = x - prev_x

            k1_Mbl, k1_Taw = self._derivatives(prev_x, Mbl, T_aw)
            Mbl_pred = Mbl + dx * k1_Mbl
            Taw_pred = T_aw + dx * k1_Taw
            k2_Mbl, k2_Taw = self._derivatives(x, Mbl_pred, Taw_pred)

            Mbl = Mbl + 0.5 * dx * (k1_Mbl + k2_Mbl)
            T_aw = T_aw + 0.5 * dx * (k1_Taw + k2_Taw)

            h = self._convective_h_gaseous(x, Mbl)
            Q_rad = self._radiative_heat_flux(x)
            T_w = T_aw + Q_rad / h if h > 0.0 else T_aw

            results.x.append(x)
            results.Mbl.append(Mbl)
            results.T_aw.append(T_aw)
            results.T_w.append(T_w)
            results.h_conv.append(h)
            results.Q_rad.append(Q_rad)

            prev_x = x

        return results


class GrissonFilmCoolingModel:
    """
    Full contour-based Grisson film cooling model attached to a pyskyfire ThrustChamber.
    """

    def __init__(self, thrust_chamber, boundary_conditions, circuit_index: int = 0):
        self.thrust_chamber = thrust_chamber
        self.boundary_conditions = boundary_conditions
        self.circuit_index = int(circuit_index)

        self._liquid_solver = LiquidFilmSolver(thrust_chamber, boundary_conditions, circuit_index)
        self._gaseous_solver = GaseousFilmSolver(thrust_chamber, boundary_conditions, circuit_index)

    def solve(self, x_array: Sequence[float]) -> tuple[LiquidFilmResults, GaseousFilmResults]:
        liquid_results = self._liquid_solver.solve(x_array)

        if liquid_results.x_dryout is None:
            gaseous_results = GaseousFilmResults()
        else:
            x_dryout = liquid_results.x_dryout
            perimeter_dryout = math.pi * self._liquid_solver.geom_eval.diameter(x_dryout)
            Mc_flux = self._liquid_solver.film.coolant_mass_flow_rate / perimeter_dryout
            gaseous_results = self._gaseous_solver.solve(x_array, x_dryout, Mc_flux)

        return liquid_results, gaseous_results