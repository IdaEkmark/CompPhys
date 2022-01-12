#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

# Position space
filename = "pDensPos"
axisfontsize = 14
array = np.genfromtxt('1/' + filename + '.csv', delimiter=',', skip_header=1)

x = array[:, 0]
pDens = array[:, 1]

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(x, pDens)

ax.set_xlabel('$x$ (Å)', fontsize=axisfontsize)
ax.set_ylabel('Probability density', fontsize=axisfontsize)
ax.set_title('Probability density as a function of position x', fontsize=axisfontsize)

ax.set_xlim([-5, 5])
ax.grid()

plt.show()
fig.savefig('1/' + filename + '.pdf')

# Momentum space, theoretical
filename = "pDensMomTheory"
axisfontsize = 14
array = np.genfromtxt('1/' + filename + '.csv', delimiter=',', skip_header=1)

p = array[:, 0]
pDens = array[:, 1]

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(p, pDens)

ax.set_xlabel('$p$ (eV / (Å/fs))', fontsize=axisfontsize)
ax.set_ylabel('Probability density', fontsize=axisfontsize)
ax.set_title('Probability density as a function of momentum p', fontsize=axisfontsize)

ax.set_xlim([0, 10])
ax.grid()

plt.show()
fig.savefig('1/' + filename + '.pdf')

# Momentum space, calculated
filename = "pDensMomCalc"
axisfontsize = 14
array = np.genfromtxt('1/' + filename + '.csv', delimiter=',', skip_header=1)

p = array[:, 0]
pDens = array[:, 1]

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(p, pDens)

ax.set_xlabel('$p$ (eV / (Å/fs))', fontsize=axisfontsize)
ax.set_ylabel('Probability density', fontsize=axisfontsize)
ax.set_title('Calculated probability density as a function of momentum p', fontsize=axisfontsize)

ax.set_xlim([0, 10])
ax.grid()

plt.show()
fig.savefig('1/' + filename + '.pdf')