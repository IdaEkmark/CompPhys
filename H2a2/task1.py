#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
filename = "TP"
axisfontsize = 14
array = np.genfromtxt('1/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
P = array[:, 1]

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(T, P)

ax.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax.set_ylabel('$P$', fontsize=axisfontsize)
ax.grid()

plt.show()
fig.savefig('1/' + filename + '.pdf')

filename = "TU"
axisfontsize = 14
array = np.genfromtxt('1/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
U = array[:, 1]
C = np.gradient(U)
fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12,5))
ax1.plot(T, U)
ax2.plot(T, C)

ax1.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax1.set_ylabel('$U$ (eV)', fontsize=axisfontsize)
ax1.grid()
ax2.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax2.set_ylabel('$C$ (eV/K)', fontsize=axisfontsize)
ax2.grid()

plt.show()
fig.savefig('1/' + filename + '.pdf')
