#!/usr/bin/env python
###############################################################################
# Modified version of E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
filename = "E1u6_dt2.5e4"
array = np.genfromtxt('6/' + filename + '.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots(figsize=(11,7))
ax.plot(array[:, 0], array[:, 1], label = "$q_1=$C")
ax.plot(array[:, 0], array[:, 3], label = "$q_2=$O")
ax.plot(array[:, 0], array[:, 5], label = "$q_3=$C")

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Displacement (Å)')
ax.set_title('$dt=2.5e4$')
ax.legend(loc='best')
ax.grid()

fig.savefig('6/' + filename + '.pdf')