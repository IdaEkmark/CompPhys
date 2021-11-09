#!/usr/bin/env python
###############################################################################
# Modified version of E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
filename = "E1u5"
array = np.genfromtxt('./' + filename + '.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots(figsize=(11,7))
ax.plot(array[:, 0], array[:, 1])
ax.plot(array[:, 0], array[:, 3])
ax.plot(array[:, 0], array[:, 5])

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Displacement (Å)')
ax.grid()

fig.savefig('./' + filename + '.pdf')