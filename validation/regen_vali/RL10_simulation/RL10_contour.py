import pyskyfire as psf
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj


""" First import RL10 engine info, and plot """
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_file_path = os.path.join(script_dir, 'RL10_data\\RL10_contour.csv')

xs_metric = []
rs_metric = []

# Open and read the CSV file
with open(csv_file_path, mode='r', newline='') as csv_file:
    reader = csv.DictReader(csv_file)
    for row in reader:
        xs_metric.append(float(row['x_metric']))
        rs_metric.append(float(row['r_metric']))

# Create the Contour object using the read data
RL10_contour = psf.regen.Contour(xs_metric, rs_metric, name="RL10 Contour")

wall = psf.regen.Wall(material = psf.common.material.StainlessSteel304, thickness = 0.3e-3) # Use the built in C106 copper data


# The following are the inputs that generate a contour remarkably similar to the actual RL10 contour. 
p_c = 32.7501e5
p_e = 0.0377e5 
MR = 5.0
L_star = 0.95
r_c = 0.123

fu = {"cea": "LH2", 
      "cantera": "H2",
      "coolprop": "hydrogen"}
#fu = {"cea": "CH4", 
#      "cantera": "CH4",
#      "coolprop": "methane"}

ox = {"cea": "LOX", 
      "cantera": "O2",
      "coolprop": "oxygen"}

theta_conv = 25
R_1f = 1.5
R_2f = 3
R_3f = 0.5
#eps = 58
length_fraction = 0.713
#length_fraction = 0.8
g = 9.81
F = 73.4e3
T_coolant_in = 33 # Kelvin
p_coolant_in = 69e5 # 68 bar

optimals = psf.skycea.OptimalValues(ox["cea"], fu["cea"], F, MR, p_c, p_e, L_star)

V_c = optimals.V_c_opt
r_t = optimals.r_t_opt
eps = optimals.eps_opt

RL10_xs, RL10_rs = psf.regen.contour_2.get_contour_2(V_c = V_c, 
                                            r_t = r_t, 
                                            area_ratio = eps, 
                                            AR_c = 1.45, #r_c = r_c, 
                                            theta_conv = theta_conv,
                                            nozzle = "rao", 
                                            R_1f=R_1f,
                                            R_2f=R_2f,
                                            R_3f=R_3f,
                                            length_fraction = length_fraction
                                            )

# TODO: add step to plot calculated vs real contour

RL10_replicate_contour = psf.regen.Contour(RL10_xs, RL10_rs, name = "Replicate Contour")

psf.regen.plot.plot_contour([RL10_contour, RL10_replicate_contour])
plt.show()
