#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
filename = "block"
array = np.genfromtxt('4/' + filename + '.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots(figsize=(11,7))
ax.plot(array[:, 0], array[:, 1])

ax.set_xlabel('B')
ax.set_ylabel('s')
ax.grid()

j = 0
s_mean = 0
for i in range(48):
    if array[i, 0] >= 1000 and array[i, 0] <= 100000:
        j += 1
        s_mean += array[i, 1]
s_mean /= j

ax.set_xscale('log')
ax.set_title('$s\\approx$' + str(np.rint(s_mean)))

fig.savefig('4/' + filename + '.pdf')

filename = "phi"
array = np.genfromtxt('4/' + filename + '.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots(figsize=(11,7))
ax.plot(np.linspace(1,181,181), array[:, 0])

ax.set_ylabel('$\phi$')
ax.set_xlabel('$s$')
ax.grid()

#ax.set_xscale('log')

fig.savefig('4/' + filename + '.pdf')