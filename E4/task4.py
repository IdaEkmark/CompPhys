import numpy as np
import matplotlib.pyplot as plt

tau_list  = [48.5e-6, 147.3e-6]
color = ['k', 'r--']
fig, ax = plt.subplots(1, 2, figsize=(15,7))
i=0
for tau in tau_list:
    filename = 'power_' + str(tau*1e6) + 'e-6'

    data = np.genfromtxt('4/' + filename + '.csv', delimiter=',', skip_header=1)

    f = data[:, 0]/(2*np.pi)
    P = data[:, 1]

    ax[i].plot(f, P, 'k')
    ax[i].set_xlabel('f [Hz]')
    ax[i].set_ylabel('Signal')
    #ax[i].legend()
    ax[i].set_title('$\\tau = $' + str(tau*1e6) + ' $\mu$s')
    #ax[i].set_xlim([0,10e3])
    i += 1
fig.savefig('4/power.pdf')
plt.show()