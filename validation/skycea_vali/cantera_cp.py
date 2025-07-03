import cantera as ct
from rocketcea.cea_obj import CEA_Obj
import math

# -----------------------------------------------------------------------------
# 1) DEFINE CONDITIONS & CONVERSIONS
# -----------------------------------------------------------------------------
Pc_Pa      = 5.0e6           # 5 MPa chamber pressure
phi        = 1.0             # stoichiometric CH₄/O₂
MR         = 4.0 / phi       # O/F mass ratio for RocketCEA (pure O₂ → MR=4)
PSI_PER_PASCAL   = 1.0/6894.757
Pc_psia    = Pc_Pa * PSI_PER_PASCAL
BTU_LBM_R_TO_J_KG_K = 1055.0559 / 0.45359237 / (5.0/9.0)

# -----------------------------------------------------------------------------
# 2) CANETRA (equilibrium cp via finite ΔT + frozen cp at eq. comp.)
# -----------------------------------------------------------------------------
gas = ct.Solution('gri30.yaml')       # GRI‑Mech 3.0 mechanism

# 2a) find the equilibrium temperature at Pc, φ
gas.TP = 3513.21, Pc_Pa               # initial guess T=3000 K
gas.set_equivalence_ratio(phi, 'CH4:1', 'O2:1')
gas.equilibrate('TP')
T_eq = gas.T

# 2b) frozen-specific-heat at that equilibrium composition
cp_frozen = gas.cp_mass              # J kg⁻¹ K⁻¹

# 2c) finite-difference along the equilibrium path
dT = 1.0
gas.TP = T_eq - dT/2.0, Pc_Pa; gas.equilibrate('TP'); h_lo = gas.enthalpy_mass
gas.TP = T_eq + dT/2.0, Pc_Pa; gas.equilibrate('TP'); h_hi = gas.enthalpy_mass
cp_eq = (h_hi - h_lo)/dT

# -----------------------------------------------------------------------------
# 3) ROCKET CEA
# -----------------------------------------------------------------------------
cea = CEA_Obj(oxName='LOX', fuelName='CH4')
cp_btu = cea.get_Chamber_Cp(Pc=Pc_psia, MR=MR, eps=40.0, frozen=0)
s = cea.get_full_cea_output(Pc=Pc_psia, MR=MR, eps=40.0, frozen=0)
print(s)
cp_btu_frz = cea.get_Chamber_Cp(Pc=Pc_psia, MR=MR, eps=40.0, frozen=1)
cp_cea = cp_btu * BTU_LBM_R_TO_J_KG_K
cp_cea_frz = cp_btu_frz * BTU_LBM_R_TO_J_KG_K

# -----------------------------------------------------------------------------
# 4) COMPARE
# -----------------------------------------------------------------------------
print(f"Equilibrium T              = {T_eq:8.1f}  K")
print(f"Cantera  frozen Cp         = {cp_frozen:8.1f}  J/(kg·K)")
print(f"Cantera  equilibrium Cp    = {cp_eq:8.1f}  J/(kg·K)")
print(f"RocketCEA frozen Cp        = {cp_cea_frz:8.1f}  J/(kg·K)")
print(f"RocketCEA equilibrium Cp   = {cp_cea:8.1f}  J/(kg·K)")
print(f"Rel. diff (Cantera vs CEA) = {(cp_eq-cp_cea)/cp_cea:+6.3%}")
print(f"Δ Cp (frozen → eq, Cantera)= {(cp_eq-cp_frozen)/cp_frozen:+6.3%}")
