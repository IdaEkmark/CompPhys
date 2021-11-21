#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
filename = "Epots_a0"
axisfontsize = 14
array = np.genfromtxt('1/' + filename + '.csv', delimiter=',', skip_header=1)

a0 = array[:, 0]
a0cube = np.power(a0, 3)
Epot = array[:, 1]/64.0

fig, ax = plt.subplots(figsize=(7.5,5))
ax.plot(a0cube, Epot)

ax.set_xlabel('$a_0^3$ ($\AA^3$)', fontsize=axisfontsize)
ax.set_ylabel('$E_{pot}$ (eV/unit cell)', fontsize=axisfontsize)
ax.grid()

ax.set_xlim([63,69])
ax.set_ylim([-13.45,-13.4])
fig.savefig('1/' + filename + 'cube' + '.pdf')

fig, ax = plt.subplots(figsize=(7.5,5))
ax.plot(a0, Epot)

ax.set_xlabel('$a_0$ ($\AA$)', fontsize=axisfontsize)
ax.set_ylabel('$E_{pot}$ (eV/unit cell)', fontsize=axisfontsize)

ax.set_xlim([3.95,4.1])
ax.set_ylim([-13.45,-13.4])
fig.savefig('1/' + filename + 'Zoom' + '.pdf')
