#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
filename = "E1u5_PS_q3_dt2.5e4"
array = np.genfromtxt('5/' + filename + '.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots(figsize=(11,7))
ax.plot(array[:, 0], array[:, 1])

ax.set_xlabel('Frequency (ps$^{-1}$)')
ax.set_ylabel('Signal (arb. unit)')
ax.set_title('$dt=2.5e4$')
ax.grid()

ax.set_xlim([-100,100])

fig.savefig('5/' + filename + '.pdf')
