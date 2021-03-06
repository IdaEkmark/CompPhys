import numpy as np
import matplotlib.pyplot as plt

dt = 1e-3
folder = '3/'
filename_T = "T_dt" + str(dt)
filename_P = "P_dt" + str(dt) + "_init"

array_T = np.genfromtxt(folder + filename_T + '.csv', delimiter=',', skip_header=1)
array_P = np.genfromtxt(folder + filename_P + '.csv', delimiter=',', skip_header=1)

t_T = array_T[:, 0]
T_inst   = array_T[:, 1] - 273.15
t_P = array_P[:, 0]
P_inst   = array_P[:, 1] / (6.24e-7)

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(t_T, T_inst)
ax2.plot(t_P, P_inst)
plt.show()