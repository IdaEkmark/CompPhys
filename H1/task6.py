import numpy as np
import matplotlib.pyplot as plt

dt = 1e-3
folder = "6/"
filename = "VCF_dt" + str(dt) + '_standard'

N = 256

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t_s = array[:, 0]
VCF_s = array[:, 1]
'''
fig, ax = plt.subplots(figsize=(15,7))
ax.plot(t, VCF)
ax.set_xlabel('$t$ (ps)')
ax.set_ylabel('$\Phi$ for standard algortihm')
ax.set_ylim([-20,20])
ax.grid()
plt.show()
fig.savefig(folder + filename + '.pdf')
'''
filename = "VCF_dt" + str(dt) + '_standard'

filename = "VCF_dt" + str(dt) + '_fast'

array2 = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t_f = array2[:, 0]
VCF_f = array2[:, 1]

fig, ax = plt.subplots(figsize=(15,7))
ax.plot(t_s, VCF_s)
ax.plot(t_f, VCF_f)
ax.legend('Standard', 'Fast')
ax.set_xlabel('$t$ (ps)')
ax.set_ylabel('$\Phi$ for standard algortihm')
ax.set_ylim([-20,20])
ax.set_xlim([0,2])
ax.grid()
plt.show()
fig.savefig(folder + filename + '.pdf')
