import cantera as ct
from rocketcea.cea_obj import CEA_Obj
import numpy as np
import math
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# 0) SETTINGS
# -----------------------------------------------------------------------------
Pc_Pa    = 5.0e6               # 5 MPa
eps      = 40.0                # nozzle area ratio (only for CEA)
dT       = 0.01                 # ΔT [K] for finite diff

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


#T_guess  = 3500.0              # initial guess for Cantera equil solve

# conversions
PSI_PER_PASCAL      = 1.0/6894.757
BTU_LBM_TO_J_PER_KG = 1055.0559/0.45359237
BTU_LBM_R_TO_J_KG_K = BTU_LBM_TO_J_PER_KG/(5.0/9.0)
Pc_psia = Pc_Pa * PSI_PER_PASCAL

# phi sweep
phis = np.linspace(0.5, 2.0, 25)   # from 0.5 to 2.0 in 16 steps

# prepare storage
cp_cantera = []; cv_cantera = []; gamma_cantera = []
cp_cea     = []; cv_cea     = []; gamma_cea     = []

# instantiate cantera gas and CEA object
gas = ct.Solution('gri30.yaml')
cea = CEA_Obj(oxName='LOX', fuelName='CH4')

for phi in phis:
    # ---- Cantera equilibrium at this phi ----
    gas.set_equivalence_ratio(phi, 'CH4:1', 'O2:1')
    gas.equilibrate('HP')          # constant enthalpy, constant pressure

    T_eq = flame_T(P=Pc_Pa, T_fuel=111.6, T_ox=90.1, phi=phi)


    # frozen cp/cv at eq. comp:
    cp_fz = gas.cp_mass
    cv_fz = gas.cv_mass

    # finite-difference for eq. cp & cv
    gas.TP = T_eq - dT/2, Pc_Pa; gas.equilibrate('TP')
    h_lo = gas.enthalpy_mass; u_lo = gas.int_energy_mass
    gas.TP = T_eq + dT/2, Pc_Pa; gas.equilibrate('TP')
    h_hi = gas.enthalpy_mass; u_hi = gas.int_energy_mass

    cp_eq = (h_hi - h_lo)/dT
    cv_eq = (u_hi - u_lo)/dT
    gamma_eq = cp_eq/cv_eq

    cp_cantera.append(cp_eq)
    cv_cantera.append(cv_eq)
    gamma_cantera.append(gamma_eq)

    # ---- RocketCEA equilibrium at this phi ----
    MR = 4.0/phi   # pure O2/CH4 stoich = 4 mass ratio
    # equilibrium cp [BTU/lbm-R] → [J/kg-K]
    cp_btu = cea.get_Chamber_Cp(Pc=Pc_psia, MR=MR, eps=eps, frozen=0)
    cp_SI  = cp_btu * BTU_LBM_R_TO_J_KG_K

    # gamma from CEA
    molwt, gam = cea.get_Chamber_MolWt_gamma(Pc=Pc_psia, MR=MR)
    cv_SI  = cp_SI/gam

    cp_cea.append(cp_SI)
    cv_cea.append(cv_SI)
    gamma_cea.append(gam)

# PLOTTING
plt.figure()
plt.plot(phis, cp_cantera, label='Cantera cp')
plt.plot(phis, cp_cea,     label='CEA cp')
plt.xlabel('phi')
plt.ylabel('cp [J/kg/K]')
plt.legend()
plt.title('Equilibrium cp vs phi')

plt.figure()
plt.plot(phis, cv_cantera, label='Cantera cv')
plt.plot(phis, cv_cea,     label='CEA cv')
plt.xlabel('phi')
plt.ylabel('cv [J/kg/K]')
plt.legend()
plt.title('Equilibrium cv vs phi')

plt.figure()
plt.plot(phis, gamma_cantera, label='Cantera gamma')
plt.plot(phis, gamma_cea,     label='CEA gamma')
plt.xlabel('phi')
plt.ylabel('gamma (cp/cv)')
plt.legend()
plt.title('Equilibrium gamma vs phi')

plt.show()
