#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
filename = "E1u5_PS_q3"
array = np.genfromtxt('./' + filename + '.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots(figsize=(11,7))
ax.plot(array[:, 0], array[:, 1])

ax.set_xlabel('Frequency (arb. unit)')
ax.set_ylabel('Signal (arb. unit)')
ax.grid()

ax.set_xlim([-100,100])

fig.savefig('./' + filename + '.pdf')