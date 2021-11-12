#!/usr/bin/env python
###############################################################################
# Modified version of E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

filename = "KPE_0.01"
folder = '3/'
array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)
modeNums = range(1,33)
maxtIndex = 2500

fig, ax = plt.subplots(figsize=(11, 7))

for k in modeNums:
    #ax.plot(array[:, 0], array[:, 1 + 3*(k-1)], label = "Kinetic energy")
    #ax.plot(array[:, 0], array[:, 2 + 3*(k-1)], label = "Potential energy")
    ax.plot(array[:, 0], array[:, 3 + 3*(k-1)])#, label = "Total energy, mode " + str(k))

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Energy (eV)')
ax.set_title('Total energy, $\\alpha = 0.01$')
#ax.legend(loc='best')
ax.set_yscale('log')
ax.grid()

fig.savefig(folder + 'TotalEnergy_alpha0.01.pdf')
'''
maxtIndex = 2500
array1 = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)
array2 = np.genfromtxt('3/' + filename + '.csv', delimiter=',', skip_header=1)
k = 1
diff = array1[:maxtIndex + 1, 1 + 3 * (k - 1)]-array2[:maxtIndex + 1, 1 + 3 * (k - 1)]
print(str(diff))
'''