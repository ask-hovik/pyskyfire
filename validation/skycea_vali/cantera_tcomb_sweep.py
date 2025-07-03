#!/usr/bin/env python3
# phi_sweep_comparison.py
#
# Compare adiabatic flame temperature vs. equivalence ratio (phi)
# for CH4/O2 using Cantera (robust multiphase) and RocketCEA.

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
from rocketcea.cea_obj_w_units import CEA_Obj

# --- helper functions from robust_flame_T.py --------------------------
def h_real(species, T, P):
    eos = {"H2": ct.Methane, "O2": ct.Oxygen}.get(species)
    if eos is None:
        raise ValueError(f"No EOS for pure {species}")
    s = eos()
    s.TP = T, P
    return s.enthalpy_mass

def flame_T_cantera(P, T_fuel, T_ox, phi,
                    fuel_species="H2", oxidizer_species="O2",
                    mech="gri30.yaml", tol=0.5, max_iter=60):
    ref = ct.Solution(mech)
    ref.set_equivalence_ratio(1.0, fuel_species, f"{oxidizer_species}:1")
    Yref = ref.Y
    m_f_st = Yref[ref.species_index(fuel_species)]
    AF_st = (1 - m_f_st) / m_f_st

    m_fuel = 1.0
    m_ox = AF_st / phi

    H_react = (m_fuel * h_real(fuel_species, T_fuel, P) +
               m_ox   * h_real(oxidizer_species, T_ox, P))
    m_tot = m_fuel + m_ox
    h_target = H_react / m_tot

    gas = ct.Solution(mech)
    gas.set_equivalence_ratio(phi, fuel_species, f"{oxidizer_species}:1")
    X0 = gas.X.copy()

    def h_eq(T):
        gas.TPX = T, P, X0
        gas.equilibrate("TP")
        return gas.enthalpy_mass

    T_lo, T_hi = 150.0, 4500.0
    f_lo = h_eq(T_lo) - h_target
    f_hi = h_eq(T_hi) - h_target
    if f_lo * f_hi > 0:
        raise RuntimeError("Cannot bracket flame temperature")

    for _ in range(max_iter):
        T_mid = 0.5*(T_lo + T_hi)
        f_mid = h_eq(T_mid) - h_target
        if abs(f_mid) < tol:
            return T_mid
        if f_lo * f_mid < 0:
            T_hi, f_hi = T_mid, f_mid
        else:
            T_lo, f_lo = T_mid, f_mid

    raise RuntimeError("Bisection did not converge")

def flame_T_cea(P, MR):
    cea = CEA_Obj(oxName='O2', fuelName='H2',
                  isp_units='sec', cstar_units='m/s',
                  pressure_units='Pa', temperature_units='K',
                  sonic_velocity_units='m/s', enthalpy_units='J/kg',
                  density_units='kg/m^3', specific_heat_units='J/kg-K',
                  viscosity_units='millipoise',
                  thermal_cond_units='W/cm-degC',
                  fac_CR=None, make_debug_prints=False)
    return cea.get_Tcomb(Pc=P, MR=MR)

# --- main sweep ----------------------------------------------------------
if __name__ == "__main__":
    P       = 32e5     # chamber pressure [Pa]
    T_fuel  = 298.15    # inlet CH4 temperature [K]
    T_ox    = 298.15    # inlet O2 temperature [K]
    phi_vals = np.linspace(0.1, 5, 50)

    T_cantera = []
    T_cea     = []
    for phi in phi_vals:
        T_cantera.append(flame_T_cantera(P, T_fuel, T_ox, phi))
        MR = 7.94 / phi    # pure O2/CH4 stoich mass ratio =4
        T_cea.append(flame_T_cea(P, MR))

    # plot
    plt.plot(phi_vals, T_cantera, linewidth=1)
    plt.plot(phi_vals, T_cea, linestyle="--", linewidth=1)
    plt.xlabel("Equivalence ratio, φ")
    plt.ylabel("Adiabatic flame temperature [K]")
    plt.title("CH₄/O₂ flame Tₐ vs. φ: Cantera vs. RocketCEA")
    plt.legend(["Cantera", "RocketCEA"])
    plt.grid(True)
    plt.tight_layout()
    plt.show()
