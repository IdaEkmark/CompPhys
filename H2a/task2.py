#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
axisfontsize = 14

filename = "TPmean"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

filename = "Trmean"
array2 = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
P = array[:, 1]
r = array2[:, 1]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7,5))
ax1.plot(T, P)

ax1.set_xlabel('$T$ ($^o$C)', fontsize=axisfontsize)
ax1.set_ylabel('$P$', fontsize=axisfontsize)
ax1.grid()

ax2.plot(T, r)

ax2.set_xlabel('$T$ ($^o$C)', fontsize=axisfontsize)
ax2.set_ylabel('$r$', fontsize=axisfontsize)
ax2.grid()

plt.show()
filename = 'TPrmean'
fig.savefig('2/' + filename + '.pdf')

filename = "TP"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

filename = "Tr"
array2 = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
P = array[:, 1]
r = array2[:, 1]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7,5))
ax1.plot(T, P)

ax1.set_xlabel('$T$ ($^o$C)', fontsize=axisfontsize)
ax1.set_ylabel('$P$', fontsize=axisfontsize)
ax1.grid()

ax2.plot(T, r)

ax2.set_xlabel('$T$ ($^o$C)', fontsize=axisfontsize)
ax2.set_ylabel('$r$', fontsize=axisfontsize)
ax2.grid()

plt.show()
filename = 'TPr'
fig.savefig('2/' + filename + '.pdf')

filename = "TU"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

filename = "TC"
array2 = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
U = array[:, 1]
C = array2[:, 1]
fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12,5))
ax1.plot(T, U)
ax2.plot(T, C)

ax1.set_xlabel('$T$ ($^o$C)', fontsize=axisfontsize)
ax1.set_ylabel('$U$ (eV)', fontsize=axisfontsize)
ax1.grid()
ax2.set_xlabel('$T$ ($^o$C)', fontsize=axisfontsize)
ax2.set_ylabel('$C$ (eV/K)', fontsize=axisfontsize)
ax2.grid()

plt.show()
fig.savefig('2/' + filename + '.pdf')
