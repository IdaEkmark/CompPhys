#!/usr/bin/env python
###############################################################################
# Modified version of E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

filename = "KPE"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)
modeNums = [1, 2, 3, 4, 5]
maxtIndex = 2500

fig, ax = plt.subplots(figsize=(11,7))
for k in modeNums:
    ax.plot(array[:maxtIndex + 1, 0], array[:maxtIndex + 1, 1 + 3*(k-1)], label = "Kinetic energy, mode " + str(k))
    ax.plot(array[:maxtIndex + 1, 0], array[:maxtIndex + 1, 2 + 3*(k-1)], label = "Potential energy, mode " + str(k))
    ax.plot(array[:maxtIndex + 1, 0], array[:maxtIndex + 1, 3 + 3*(k-1)], label = "Total energy, mode " + str(k))

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Energy (eV)')
ax.set_title('$\Delta t = 0.1$')
ax.legend(loc='best')
ax.grid()

fig.savefig('2/' + filename + '_1to5.pdf')
'''
maxtIndex = 2500
array1 = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)
array2 = np.genfromtxt('3/' + filename + '.csv', delimiter=',', skip_header=1)
k = 1
diff = array1[:maxtIndex + 1, 1 + 3 * (k - 1)]-array2[:maxtIndex + 1, 1 + 3 * (k - 1)]
print(str(diff))
'''