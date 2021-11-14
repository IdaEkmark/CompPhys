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
folder = '4/'
array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)
modeNums = range(1,33)
maxtIndex = -2

fig, ax = plt.subplots(figsize=(11, 7))
sum = 0
for k in modeNums:
    #ax.plot(array[:, 0], array[:, 1 + 3*(k-1)], label = "Kinetic energy")
    #ax.plot(array[:, 0], array[:, 2 + 3*(k-1)], label = "Potential energy")
    ax.plot(array[:, 0].flatten(), array[:, 3 + 3*(k-1)].flatten())#, label = "Total energy, mode " + str(k))
    sum += array[:, 3 + 3 * (k - 1)]
    #print('Mode ' + str(k))
    #print('Array:\n' + str(array[:, 3 + 3*(k-1)]))
    #print('Sum:\n' + str(sum) + '\n\n')
ax.plot(array[:, 0].flatten(), sum.flatten())#, label = "Total energy, mode " + str(k))
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Energy (eV)')
ax.set_title('Total energy, $\\alpha = 0.01$')
#ax.legend(loc='best')
#ax.set_yscale('log')
ax.grid()
plt.show()
fig.savefig(folder + 'TotalEnergy_alpha0.01.pdf')
'''
maxtIndex = 2500
array1 = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)
array2 = np.genfromtxt('3/' + filename + '.csv', delimiter=',', skip_header=1)
k = 1
diff = array1[:maxtIndex + 1, 1 + 3 * (k - 1)]-array2[:maxtIndex + 1, 1 + 3 * (k - 1)]
print(str(diff))
'''