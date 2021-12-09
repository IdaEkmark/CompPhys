import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

tau_list  = [48.5e-6, 147.3e-6]
fig, ax = plt.subplots(2, 2, figsize=(12,12))
i=0
t = 2
for tau in tau_list:
    filenameV = 'rhoV_tau' + str(tau*1e6) + 'e-6_t' + str(round(t,2))
    filenameX = 'rhoX_tau' + str(tau*1e6) + 'e-6_t' + str(round(t,2))

    dataV = np.genfromtxt('3/' + filenameV + '.csv', delimiter=',', skip_header=1)
    dataX = np.genfromtxt('3/' + filenameX + '.csv', delimiter=',', skip_header=1)

    v = dataV[:, 0]*1e3
    rhoV  = dataV[:, 1]
    x = dataX[:, 0]*1e9
    rhoX  = dataX[:, 1]

    ax[0, i].plot(x, rhoX)
    ax[0, i].set_xlabel('$x$ [nm]')
    ax[0, i].set_ylabel('$\\rho_x$')
    #ax[0, i].set_xlim([0, 2])
    ax[1, i].plot(v, rhoV)
    ax[1, i].set_xlabel('$v$ [mm/s]')
    ax[1, i].set_ylabel('$\\rho_v$')
    #ax[1, i].set_xlim([0, 2])
    ax[0, i].set_title('$\\tau = $' + str(tau*1e6) + ' $\mu$s')
    i += 1
plt.suptitle('t=' + str(round(t,2)))
fig.savefig('3/density_t' + str(round(t,2)) + '.pdf')
plt.show()