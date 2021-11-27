import numpy as np
import matplotlib.pyplot as plt

axisfontsize = 14
dt = 5e-4
folder = "6/"
filename = "VCF_dt" + str(dt) + '_standard'

N = 256

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t_s = array[:, 0]
VCF_s = array[:, 1]

filename = "VCF_dt" + str(dt) + '_fast'

array2 = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t_f = array2[:, 0]
VCF_f = array2[:, 1]

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(t_s, VCF_s)
ax.plot(t_f, VCF_f)
ax.legend('Standard Correlation Algorithm', 'Fast Correlation Algorithm', loc = 'best', fontsize=axisfontsize)
ax.set_xlabel('$t$ (ps)', fontsize=axisfontsize)
ax.set_ylabel('$\Phi(t)$ ($\mathrm{\AA}^2$/ps$^2$)', fontsize=axisfontsize)
ax.set_ylim([-20,20])
ax.set_xlim([0,1])
ax.grid()
plt.show()
filename = "VCF_dt" + str(dt)
fig.savefig(folder + filename + '.pdf')
