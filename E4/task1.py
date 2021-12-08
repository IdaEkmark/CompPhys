import numpy as np
import matplotlib.pyplot as plt

dt_list = [5e-6, 1e-6]
tau_list  = [147.3e-6, 48.5e-6]
color = ['r', 'k']
for dt in dt_list:
    i=0
    fig, ax = plt.subplots(2, 2, figsize=(12,12))
    for tau, c in zip(tau_list, color):
        filenameV = 'velocity_tau' + str(tau*1e6) + 'e-6_dt' + str(int(dt*1e6)) + 'e-6'
        filenameX = 'position_tau' + str(tau*1e6) + 'e-6_dt' + str(int(dt*1e6)) + 'e-6'

        dataV = np.genfromtxt('1/' + filenameV + '.csv', delimiter=',', skip_header=1)
        dataX = np.genfromtxt('1/' + filenameX + '.csv', delimiter=',', skip_header=1)

        tV = dataV[:, 0]*1e3
        V  = dataV[:, 1]*1e3
        tX = dataX[:, 0]*1e3
        X  = dataX[:, 1]*1e9

        ax[0, i].plot(tX, X, color=c, marker='o', markersize=2, linestyle='-') # (tX, X, '-o'+c, markersize=2)
        ax[0, i].set_xlabel('Time [ms]')
        ax[0, i].set_ylabel('Position [nm]')
        ax[0, i].set_xlim([0, 2])
        ax[1, i].plot(tV, V, color=c, marker='o', markersize=2, linestyle='-') # (tV, V, '-o'+c, markersize=2)
        ax[1, i].set_xlabel('Time [ms]')
        ax[1, i].set_ylabel('Velocity [mm/s]')
        ax[1, i].set_xlim([0, 2])
        ax[0,i].set_title('$\\tau = $' + str(tau*1e6) + ' $\mu$s')
        i += 1
    plt.suptitle('d$t = $' + str(dt*1e6) + ' $\mu$s')
    fig.savefig('1/VP_dt'+str(dt)+'.pdf')
    plt.show()