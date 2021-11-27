import numpy as np
import matplotlib.pyplot as plt

axisfontsize = 14
dt = 1e-3
filename = "MSD_dt" + str(dt) + '_solid'
folder = "5/"
N = 256

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t = array[:, 0]
MSD = array[:, 1]

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(t, MSD)
ax.set_xlabel('$t$ (ps)', fontsize=axisfontsize)
ax.set_ylabel('$\Delta_\mathrm{MSD}$ ($\mathrm{\AA}^2$)', fontsize=axisfontsize)
ax.set_title('Solid aluminium', fontsize=axisfontsize)
ax.grid()
plt.show()
fig.savefig(folder + filename + '.pdf')

filename = "MSD_dt" + str(dt) + '_liquid'

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t = array[:, 0]
MSD = array[:, 1]

Ds = np.mean(MSD[-100:]/(6*t[-100:]))

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(t, MSD)
ax.set_xlabel('$t$ (ps)', fontsize=axisfontsize)
ax.set_ylabel('$\Delta_\mathrm{MSD}$ ($\mathrm{\AA}^2$)', fontsize=axisfontsize)
ax.set_title('Liquid aluminium, $D_\mathrm{s}=$' + str(np.round(Ds,3)) + ' $\mathrm{\AA}^2$/ps', fontsize=axisfontsize)
ax.grid()
plt.show()
fig.savefig(folder + filename + '.pdf')
