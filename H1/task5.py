import numpy as np
import matplotlib.pyplot as plt

dt = 1e-3
filename = "MSD_dt" + str(dt) + '_solid'
folder = "5/"
N = 256

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t = array[:, 0]
MSD = array[:, 1]

fig, ax = plt.subplots(figsize=(15,7))
ax.plot(t, MSD)
ax.set_xlabel('$t$ (ps)')
ax.set_ylabel('$\Delta_{MSD}$ for solid (A)')
ax.grid()
plt.show()
fig.savefig(folder + filename + '.pdf')

filename = "MSD_dt" + str(dt) + '_liquid'

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t = array[:, 0]
MSD = array[:, 1]

fig, ax = plt.subplots(figsize=(15,7))
ax.plot(t, MSD)
ax.set_xlabel('$t$ (ps)')
ax.set_ylabel('$\Delta_{MSD}$ for liquid (A)')
ax.grid()
plt.show()
fig.savefig(folder + filename + '.pdf')
