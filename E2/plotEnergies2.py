#!/usr/bin/env python
###############################################################################
# Modified version of E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import sys
np.set_printoptions(threshold=sys.maxsize)
# skip_header skips the first
# row in data.csv

filename = "KPE_0.01"
folder = '3/'
array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)
modeNums = range(1,33)
maxtIndex = -2

fig, ax = plt.subplots(figsize=(11, 7))
sum = 0
for k in modeNums:
    ax.plot(array[:, 0].flatten(), array[:, 3 + 3*(k-1)].flatten(), label = "$E_k$, mode " + str(k))
    print('Mode ' + str(k))
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Energy (eV)')
ax.set_title('Total energy, $\\alpha = 0.01$')
#ax.legend(loc='best', prop={'size': 8})
#ax.set_yscale('log')
ax.grid()
plt.show()
fig.savefig(folder + filename + '.pdf')
'''
maxtIndex = 2500
array1 = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)
array2 = np.genfromtxt('3/' + filename + '.csv', delimiter=',', skip_header=1)
k = 1
diff = array1[:maxtIndex + 1, 1 + 3 * (k - 1)]-array2[:maxtIndex + 1, 1 + 3 * (k - 1)]
print(str(diff))
'''