#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
array = np.genfromtxt('powerspectrum26.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots()
ax.plot(array[:, 0], array[:, 1])

ax.set_xlabel('frequency (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()

fig.savefig('powerspectrum26.pdf')
