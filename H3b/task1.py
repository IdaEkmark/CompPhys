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

filename1 = "pDensPosCalc"
axisfontsize = 14
array = np.genfromtxt('1/' + filename1 + '.csv', delimiter=',', skip_header=1)

xC = array[:, 0]
pDensC = array[:, 1]

fig, ax = plt.subplots(figsize=(7,5))
#ax.plot(xC, pDensC, label="Calculated $|\psi(x)|^2$")
ax.plot(x, pDens, label="Theoretical $|\psi(x)|^2$")

ax.set_xlabel('$x$ (Å)', fontsize=axisfontsize)
ax.set_ylabel('Probability density (Å$^{-1}$)', fontsize=axisfontsize)
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


filename = "pDensMomCalc"
axisfontsize = 14
array = np.genfromtxt('1/' + filename + '.csv', delimiter=',', skip_header=1)

pC = array[:, 0]
pDensC = array[:, 1]


fig, ax = plt.subplots(figsize=(7,5))
ax.plot(p, pDens, label='Theoretical $|\phi(p)|^2$')
ax.plot(pC, pDensC, '.', label='Calculated $|\phi(p)|^2$')

ax.set_xlabel('$p$ (eV / (Å/fs))', fontsize=axisfontsize)
ax.set_ylabel('Probability density ((Å/fs) / eV)', fontsize=axisfontsize)
ax.set_title('Probability density as a function of momentum p', fontsize=axisfontsize)
ax.legend()

ax.set_xlim([0, 10])
ax.grid()

plt.show()
fig.savefig('1/' + 'pDensMom' + '.pdf')


'''
# Momentum space, calculated

fig, ax = plt.subplots(figsize=(7,5))


ax.set_xlabel('$p$ (eV / (Å/fs))', fontsize=axisfontsize)
ax.set_ylabel('Probability density', fontsize=axisfontsize)
ax.set_title('Calculated probability density as a function of momentum p', fontsize=axisfontsize)

ax.set_xlim([0, 10])
ax.grid()

plt.show()
fig.savefig('1/' + filename + '.pdf')
'''