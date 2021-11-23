import numpy as np
import matplotlib.pyplot as plt

dt = 1e-3
filename = "VCF_dt" + str(dt) + '_standard'
folder = "6/"
N = 256

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t = array[:, 0]
VCF = array[:, 1]

fig, ax = plt.subplots(figsize=(15,7))
ax.plot(t, VCF)
ax.set_xlabel('$t$ (ps)')
ax.set_ylabel('$\Phi$ for standard algortihm')
ax.set_ylim([-20,20])
ax.grid()
plt.show()
fig.savefig(folder + filename + '.pdf')

filename = "VCF_dt" + str(dt) + '_fast'

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t = array[:, 0]
VCF = array[:, 1]

fig, ax = plt.subplots(figsize=(15,7))
ax.plot(t, VCF)
ax.set_xlabel('$t$ (ps)')
ax.set_ylabel('$\Delta_{VCF}$ for liquid')
ax.grid()
plt.show()
fig.savefig(folder + filename + '.pdf')
