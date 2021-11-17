#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
dt = 2e-2
filename = "E_dt" + str(dt)
filename_p = "Epot_dt" + str(dt)
filename_k = "Ekin_dt" + str(dt)
folder = "2/"
N = 256

array_p = np.genfromtxt(folder + filename_p + '.csv', delimiter=',', skip_header=1)
array_k = np.genfromtxt(folder + filename_k + '.csv', delimiter=',', skip_header=1)

t_p = array_p[:, 0]
E_p = array_p[:, 1]
t_k = array_k[:, 0]
E_k = array_k[:, 1]
T = 2*np.mean(E_k)/(3*N) / (8.6e-5)

fig, ax = plt.subplots(figsize=(11,7))
ax.plot(t_p, E_p, label = '$E_{potential}$')
ax.plot(t_k, E_k, label = '$E_{kinetic}$')
ax.plot(t_k, E_k+E_p, label = '$E_{total}$')

ax.set_xlabel('$t$ (ps)')
ax.set_ylabel('$E$ (eV/unit cell)')
ax.grid()
ax.legend(loc='best')
ax.set_title('$T=$' + str(int(round(T,0))) + ', $dt=$' + str(dt))
plt.show()


fig.savefig(folder + filename + '.pdf')

