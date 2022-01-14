#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

axisfontsize = 14
fig, ax = plt.subplots(figsize=(9,5))

x = np.zeros((6,401))
pDens = np.zeros((6,401))

for i in range(6):
    filename = "pDensPos" + str(i) + "0fs"
    array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

    x[i,:] = array[:, 0]
    pDens[i,:] = array[:, 1]

    ax.plot(x[i,:], pDens[i,:], label=("$|\Psi(x,t)|^2$ at $t = $" + str(i*10) + " fs"))

ax.set_xlabel('$x$ (Å)', fontsize=axisfontsize)
ax.set_ylabel('Probability density (Å$^{-1}$)', fontsize=axisfontsize)
ax.set_title('Probability density as a function of position x', fontsize=axisfontsize)
ax.legend()

ax.set_xlim([-2, 12])
ax.grid()

plt.show()
fig.savefig('2/' + filename + '.pdf')