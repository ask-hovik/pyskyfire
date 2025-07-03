import cantera as ct
import numpy as np

# Get all of the Species objects defined in the GRI 3.0 mechanism
species = {S.name: S for S in ct.Species.list_from_file("gri30_highT.yaml")}

# Create an ideal gas phase object with species representing complete combustion
complete_species = [species[S] for S in ("CH4", "O2", "CO2", "H2O")]
gas1 = ct.Solution(thermo="ideal-gas", species=complete_species)

phi = np.linspace(0.2, 2.5, 500)
T_complete = np.zeros(phi.shape)
for i in range(len(phi)):
    gas1.TP = 300, ct.one_atm
    gas1.set_equivalence_ratio(phi[i], "CH4", "O2:1")
    gas1.equilibrate("HP")
    T_complete[i] = gas1.T

# Create an IdealGas object including incomplete combustion species
gas2 = ct.Solution(thermo="ideal-gas", species=species.values())
T_incomplete = np.zeros(phi.shape)
for i in range(len(phi)):
    gas2.TP = 300, ct.one_atm
    gas2.set_equivalence_ratio(phi[i], "CH4", "O2:1")
    gas2.equilibrate("HP")
    T_incomplete[i] = gas2.T

# --- plotting in the same style as plot.py -----------------------------
import matplotlib.pyplot as plt
import matplotlib as mpl



# --- plotting in the same style as plot.py -----------------------------
import matplotlib.pyplot as plt

# 1. Apply the style and font settings taken from plot.py
plt.style.use('ggplot')
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Latin Modern Roman"],
    "mathtext.fontset": "custom",             # tell mathtext to use our set
    "mathtext.rm": "Latin Modern Math",       # upright
    "mathtext.it": "Latin Modern Math:italic",# italic
    "mathtext.bf": "Latin Modern Math:bold",  # bold
    "font.size": 10,
    "savefig.dpi": 600,     # high-resolution exports
})

# 2. Create the figure / axis and make the plot
fig, ax = plt.subplots(figsize=(8, 4))          # (width, height) in inches
ax.plot(phi, T_complete,   lw=2, label="Complete Combustion")
ax.plot(phi, T_incomplete, lw=2, label="Incomplete Combustion")

# 3. Labels, legend, grid — again matching plot.py conventions
ax.set_xlabel(r"Equivalence ratio, $\varphi$")
ax.set_ylabel("Temperature [K]")
ax.legend(frameon=True, framealpha=0.9)
ax.grid(True)

fig.tight_layout()          # neat spacing à la plot.py
plt.show()

