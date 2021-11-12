#!/usr/bin/env python
###############################################################################
# Modified version of E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

filename = "KPE_0.01"
array = np.genfromtxt('3/' + filename + '.csv', delimiter=',', skip_header=1)
modeNums = range(1,33)
maxtIndex = 2500

for k in modeNums:
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.plot(array[:maxtIndex + 1, 0], array[:maxtIndex + 1, 1 + 3*(k-1)], label = "Kinetic energy")
    ax.plot(array[:maxtIndex + 1, 0], array[:maxtIndex + 1, 2 + 3*(k-1)], label = "Potential energy")
    ax.plot(array[:maxtIndex + 1, 0], array[:maxtIndex + 1, 3 + 3*(k-1)], label = "Total energy")

    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (eV)')
    ax.set_title('Mode ' + str(k) + ', $\alpha = 0.01$')
    ax.legend(loc='best')
    ax.grid()

    fig.savefig('3/' + filename + '_mode' + str(int(k)) + '.pdf')
'''
maxtIndex = 2500
array1 = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)
array2 = np.genfromtxt('3/' + filename + '.csv', delimiter=',', skip_header=1)
k = 1
diff = array1[:maxtIndex + 1, 1 + 3 * (k - 1)]-array2[:maxtIndex + 1, 1 + 3 * (k - 1)]
print(str(diff))
'''