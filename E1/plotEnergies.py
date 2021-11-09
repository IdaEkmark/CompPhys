#!/usr/bin/env python
###############################################################################
# Modified version of E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
filename = "E1u5_Energies_dt2.5e4"
array = np.genfromtxt('5/' + filename + '.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots(figsize=(11,7))
ax.plot(array[:, 0], array[:, 1], label = "Total energy")
ax.plot(array[:, 0], array[:, 2], label = "Kinetic energy")
ax.plot(array[:, 0], array[:, 3], label = "Potential energy")

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Energy (eV)')
ax.set_title('$dt=2.5e4$')
ax.legend(loc='best')
ax.grid()

fig.savefig('5/' + filename + '.pdf')