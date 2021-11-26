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

p = np.poly1d(np.polyfit(a0, Epot, 2))
xp = np.linspace(4,4.1, 1000)
minlab = 'Minimum at ' + str(np.round(xp[np.argmin(p(xp))]**3,1)) + ' $\mathrm{\AA}^3$'

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(np.array([xp[np.argmin(p(xp))]**3, xp[np.argmin(p(xp))]**3]), np.array([-14, -12]), '--', color = 'silver', label=minlab)
ax.plot(xp**3, p(xp), label='Polyfit')
ax.plot(a0cube, Epot, '.', color = 'k', label='Data')
ax.legend()

ax.set_xlabel('$a_0^3$ ($\mathrm{\AA}^3$)', fontsize=axisfontsize)
ax.set_ylabel('$E_\mathrm{pot}$ (eV/unit cell)', fontsize=axisfontsize)
ax.grid()

ax.set_xlim([64,68])
ax.set_ylim([-13.445,-13.415])
ax.set_xticks([64, 65, 66, 67, 68])
ax.set_yticks([-13.42, -13.43, -13.44])
plt.tight_layout()
plt.show()
fig.savefig('1/' + filename + 'cube' + '.pdf')

minlab = 'Minimum at ' + str(np.round(xp[np.argmin(p(xp))],2)) + ' $\mathrm{\AA}$'
fig, ax = plt.subplots(figsize=(7,5))
ax.plot(np.array([xp[np.argmin(p(xp))], xp[np.argmin(p(xp))]]), np.array([-14, -12]), '--', color = 'silver', label=minlab)
ax.plot(xp, p(xp), label='Polyfit')
ax.plot(a0, Epot, '.', color = 'k', label='Data')
ax.legend()

ax.set_xlabel('$a_0$ ($\mathrm{\AA}$)', fontsize=axisfontsize)
ax.set_ylabel('$E_\mathrm{pot}$ (eV/unit cell)', fontsize=axisfontsize)

ax.set_xlim([4,4.08])
ax.set_ylim([-13.445,-13.415])
ax.set_xticks([4.0, 4.02, 4.04, 4.06, 4.08])
ax.set_yticks([-13.42, -13.43, -13.44])
ax.grid()
plt.tight_layout()
plt.show()
fig.savefig('1/' + filename + '.pdf')
