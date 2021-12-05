#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
filename = "densityfunction"
axisfontsize = 14
array = np.genfromtxt('1/' + filename + '.csv', delimiter=',', skip_header=1)

r = array[:, 0]
rho = array[:, 1]
r2 = np.linspace(0, max(r))
Z2 = 2
rho2 = Z2**3 * 4 * r2**2 * np.exp(-2*Z2*r2)
Z3 = 27/16
rho3 = Z3**3 * 4 * r2**2 * np.exp(-2*Z3*r2)
#rho4 = (exp(-2*(sqrt(x1^2+x2^2+x3^2) + sqrt(x4^2+x5^2+x6^2))+1/(2*(0.1+1/(sqrt((x1-x4)^2+(x2-x5)^2+(x3-x6)^2))))))**2

fig, ax = plt.subplots(figsize=(7,5))
ax.scatter(r, rho)
ax.plot(r2, rho2, label='Z=2')
ax.plot(r2, rho3, label='Z=27/16')

ax.legend(fontsize=axisfontsize)
ax.set_xlabel('$r$ ($\mathrm{\AA}$)', fontsize=axisfontsize)
ax.set_ylabel('$\\rho$ (??)', fontsize=axisfontsize)
ax.grid()

plt.show()
#fig.savefig('1/' + filename + 'cube' + '.pdf')
