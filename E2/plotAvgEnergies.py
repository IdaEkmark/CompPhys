#!/usr/bin/env python
###############################################################################
# Modified version of E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

filename = "KPE_0.1"
folder = '4/'
array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)
modeNums = range(1,33)
maxtIndex = -2

dt = 100.0
avgEs = np.zeros((10001,32))
avgEs[0,:] = array[0,3::3]
for ti,t in enumerate(array[1:,0]):
    avgEs[ti+1,:] = np.mean(array[:ti+2,3::3], axis=0)

fig, ax = plt.subplots(figsize=(13,9))
for k in modeNums:
    ax.loglog(array[:maxtIndex + 1, 0], avgEs[:maxtIndex + 1, k-1], label = "Total energy, mode " + str(k))

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Energy (eV)')
ax.set_title('$\Delta t = 0.1$, $\\alpha = 0.1$')
ax.legend(loc='best', prop={'size': 8})
ax.grid()

fig.savefig(folder + filename + '_equip.pdf')
'''
maxtIndex = 2500
array1 = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)
array2 = np.genfromtxt('3/' + filename + '.csv', delimiter=',', skip_header=1)
k = 1
diff = array1[:maxtIndex + 1, 1 + 3 * (k - 1)]-array2[:maxtIndex + 1, 1 + 3 * (k - 1)]
print(str(diff))
'''
