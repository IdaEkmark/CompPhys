import numpy as np
import matplotlib.pyplot as plt

tau_list  = [48.5e-6, 147.3e-6]
color = ['r', 'k']
fig, ax = plt.subplots(1, 1, figsize=(7,5))
for tau, c in zip(tau_list, color):
    filename = 'velcorr_' + str(tau*1e6) + 'e-6'

    data = np.genfromtxt('3/' + filename + '.csv', delimiter=',', skip_header=1)

    t = data[:, 0]*1e3
    vc = data[:, 1]/data[0, 1]

    ax.plot(t, vc, c, label = '$\\tau = $' + str(tau*1e6) + ' $\mu$s')
ax.set_xlabel('t [ms]')
ax.set_ylabel('VCF normalized')
ax.legend()
ax.set_xlim([0,1])
fig.savefig('3/velcorr.pdf')
plt.show()