#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
filename = "Epot_a0"
array = np.genfromtxt('1/' + filename + '.csv', delimiter=',', skip_header=1)

a0 = array[0, :]
a0cube = a0 ** 3
Epot = array[1, :]

fig, ax = plt.subplots(figsize=(11,7))
ax.plot(a0cube, Epot)

ax.set_xlabel('$a_0^3$ (Å$^3$)')
ax.set_ylabel('$E_{\rm pot}$')
ax.grid()

#ax.set_xlim([63,69])
fig.savefig('6/' + filename + 'cube' + '.pdf')

fig, ax = plt.subplots(figsize=(11,7))
ax.plot(a0cube, Epot)

ax.set_xlabel('$a_0$ (Å)')
ax.set_ylabel('$E_{\rm pot}$')

ax.set_xlim([3,5])
fig.savefig('6/' + filename + '.pdf')
