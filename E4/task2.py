import numpy as np
import matplotlib.pyplot as plt

dtau_list = ['25']
tau_list  = [48.5e-6, 147.3e-6]
color = ['k']
fig, ax = plt.subplots(1, 2, figsize=(15,7))
for dtau, c in zip(dtau_list, color):
    i=0
    for tau in tau_list:
        filename = 'power_' + str(tau*1e6) + 'e-6_dtau' + dtau + 'dt_Long'

        data = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

        f = data[:, 0]
        P = data[:, 1]

        ax[i].plot(f, P, c, label = 'd$\\tau = $' + dtau + 'd$t$')
        ax[i].set_xlabel('f [Hz]')
        ax[i].set_ylabel('Signal')
        ax[i].legend()
        ax[i].set_title('$\\tau = $' + str(tau*1e6) + ' $\mu$s')
        ax[i].set_xlim([0,10e3])
        i += 1
fig.savefig('2/power.pdf')
plt.show()