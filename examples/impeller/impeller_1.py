from pyskyfire.pump.impeller_new import Impeller
from pyskyfire.viz.impeller_viz import plot_impeller

imp = Impeller(Q=0.05, H=20.0, n=3500)
print(imp)

plot_impeller(imp)