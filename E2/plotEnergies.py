#!/usr/bin/env python
###############################################################################
# Modified version of E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

filename = "KPE"
folder = '4'
array = np.genfromtxt(folder + '/' + filename + '_1e-1.csv', delimiter=',', skip_header=1)
modeNums = [1, 2, 3, 4, 5]
maxtIndex = -2

fig, ax = plt.subplots(figsize=(11,7))
for k in modeNums:
    #ax.plot(array[:maxtIndex + 1, 0], array[:maxtIndex + 1, 1 + 3*(k-1)], label = "Kinetic energy, mode " + str(k))
    #ax.plot(array[:maxtIndex + 1, 0], array[:maxtIndex + 1, 2 + 3*(k-1)], label = "Potential energy, mode " + str(k))
    ax.plot(array[:maxtIndex + 1, 0], array[:maxtIndex + 1, 3 + 3*(k-1)], label = "Total energy, mode " + str(k))

ax.set_xlabel('Time (FPUT units)')
ax.set_ylabel('Energy (FPUT units)')
ax.set_title('$\\alpha = 0.1$')
ax.legend(loc='best')
ax.grid()

fig.savefig(folder + '/' + 'E' + '_1to5_1e-1.pdf')
