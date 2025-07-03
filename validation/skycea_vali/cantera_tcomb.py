#!/usr/bin/env python3
# robust_flame_T.py  –  phase‑agnostic adiabatic flame temperature
#
# Works for CH4 + O2 whether they enter as liquid, two‑phase, or gas.

import cantera as ct
import math
from rocketcea.cea_obj_w_units import CEA_Obj as CEA_metric

def h_real(species, T, P):
    """Specific enthalpy [J/kg] from Cantera multiphase EOS."""

    eos = {"CH4": ct.Methane, "O2": ct.Oxygen}.get(species)
    if eos is None:
        raise ValueError(f"No EOS for pure {species}")
    
    s = eos()
    s.TP = T, P
    return s.enthalpy_mass

def flame_T(P, 
            T_fuel, 
            T_ox, 
            phi,
            fuel_species="CH4", 
            oxidizer_species="O2",
            mech="gri30.yaml",
            tol=0.5, 
            max_iter=60):
    """
    Return adiabatic flame temperature [K] for any inlet phase.
    """
    # --- stoichiometric mass ratio (kg oxid / kg fuel) -----------------
    ref = ct.Solution(mech)
    ref.set_equivalence_ratio(1.0, fuel_species, f"{oxidizer_species}:1")
    Yref = ref.Y
    m_f_st = Yref[ref.species_index(fuel_species)]
    AF_st  = (1 - m_f_st) / m_f_st

    m_fuel = 1.0
    m_ox   = AF_st / phi

    # --- reactant enthalpy ---------------------------------------------
    H_react = (m_fuel*h_real(fuel_species, T_fuel, P) + m_ox*h_real(oxidizer_species, T_ox, P)) #J
    # the above line fetches the specific(?) enthalpies of each fluid, multiplies it with the mass ratio for that species, and then adds them together
    # so it is sort of a total enthalpy at the inlet state for the mixture, for 1 kg of fuel to OF kg of oxidizer. 

    m_tot = m_fuel + m_ox
    h_target = H_react / m_tot                                  # J/kg
    # the above divides the sort of total enthalpy by the reference masses total


    # --- combustion mixture object (composition fixed by phi) ----------
    # --- combustion mixture object (composition fixed by phi) ----------
    gas = ct.Solution(mech)
    # Record the *reactant* composition (set_equivalence_ratio sets X)
    gas.set_equivalence_ratio(phi, f"{fuel_species}:1", f"{oxidizer_species}:1")
    X0 = gas.X.copy()     # freeze reactant mole fractions

    # helper to evaluate equilibrium enthalpy at (T, P)
    def h_eq(T):
        # reset state to reactants at this trial T and P
        gas.TPX = T, P, X0
        # find equilibrium composition at constant T, P
        gas.equilibrate("TP")
        # return the resulting mixture enthalpy
        return gas.enthalpy_mass

    # bracket the root
    T_lo, T_hi = 150.0, 4500.0
    
    f_lo = h_eq(T_lo) - h_target
    f_hi = h_eq(T_hi) - h_target
    if f_lo * f_hi > 0:
        raise RuntimeError("Failed to bracket the flame temperature")

    # bisection
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

# CEA: 
def flame_T_cea(P, MR):

    # CEA
    Pc_Pa      = P           # 50 barchamber pressure
    phi        = 1.0            # stoichiometric CH₄/O₂
    MR         = 4 / phi       # O/F mass ratio for RocketCEA (pure O₂ → MR=4)

    cea = CEA_metric(oxName='LOX', fuelName='CH4',
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

    s = cea.get_Tcomb(Pc=Pc_Pa, # number or list of chamber pressures
                    MR=MR)
    
    return s


# ----------------------------------------------------------------------
if __name__ == "__main__":
    P = 50e5
    MR = 4.0
    phi = 1

    Tad = flame_T(P=P, T_fuel=111.6, T_ox=90.1, phi=phi)
    Tcea = flame_T_cea(P=P, MR=MR)
    print(f"Cantera Tcomb: {Tad}")

    print(f"CEA Tcomb: {Tcea}")

