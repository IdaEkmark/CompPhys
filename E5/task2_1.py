import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

tau_list  = [48.5e-6, 147.3e-6]
fig, ax = plt.subplots(2, 2, figsize=(12,12))
for j in range(1,6):
    i=0
    for tau in tau_list:
        filenameV = 'velocity_tau' + str(tau*1e6) + 'e-6_' + str(int(j))
        filenameX = 'position_tau' + str(tau*1e6) + 'e-6_' + str(int(j))

        dataV = np.genfromtxt('2/' + filenameV + '.csv', delimiter=',', skip_header=1)
        dataX = np.genfromtxt('2/' + filenameX + '.csv', delimiter=',', skip_header=1)

        tV = dataV[:, 0]*1e3
        V  = dataV[:, 1]*1e3
        tX = dataX[:, 0]*1e3
        X  = dataX[:, 1]*1e9

        col = pl.cm.Blues([(j-1)/8])
        ax[0, i].plot(tX, X, color=col)
        ax[0, i].set_xlabel('Time [ms]')
        ax[0, i].set_ylabel('Position [nm]')
        #ax[0, i].set_xlim([0, 2])
        ax[1, i].plot(tV, V, color=col)
        ax[1, i].set_xlabel('Time [ms]')
        ax[1, i].set_ylabel('Velocity [mm/s]')
        #ax[1, i].set_xlim([0, 2])
        ax[0, i].set_title('$\\tau = $' + str(tau*1e6) + ' $\mu$s')
        i += 1
fig.savefig('2/VP.pdf')
plt.show()

tau_list  = [48.5e-6, 147.3e-6]
fig, ax = plt.subplots(2, 2, figsize=(20,12))

for j in range(1,6):
    i=0
    for tau in tau_list:
        filenameV = 'velocity_tau' + str(tau*1e6) + 'e-6_' + str(int(j))
        filenameX = 'position_tau' + str(tau*1e6) + 'e-6_' + str(int(j))

        dataV = np.genfromtxt('2/' + filenameV + '.csv', delimiter=',', skip_header=1)
        dataX = np.genfromtxt('2/' + filenameX + '.csv', delimiter=',', skip_header=1)

        tV = dataV[:, 0]*1e3
        V  = dataV[:, 1]*1e3
        tX = dataX[:, 0]*1e3
        X  = dataX[:, 1]*1e9

        col = pl.cm.Blues([(j-1)/8])
        ax[0, i].plot(tX, X, color=col)
        ax[0, i].set_xlabel('Time [ms]')
        ax[0, i].set_ylabel('Position [nm]')
        #ax[0, i].set_xlim([0, 2])
        ax[1, i].plot(tV, V, color=col)
        ax[1, i].set_xlabel('Time [ms]')
        ax[1, i].set_ylabel('Velocity [mm/s]')
        #ax[1, i].set_xlim([0, 2])
        ax[0, i].set_title('$\\tau = $' + str(tau*1e6) + ' $\mu$s')
        i += 1
for shift in [-1, 0, 1]:
    i=0
    for tau in tau_list:
        filenameVmean = 'meanvelocity_tau' + str(tau*1e6) + 'e-6'
        filenameXmean = 'meanposition_tau' + str(tau*1e6) + 'e-6'
        filenameVvar = 'varvelocity_tau' + str(tau*1e6) + 'e-6'
        filenameXvar = 'varposition_tau' + str(tau*1e6) + 'e-6'

        meanDataV = np.genfromtxt('2/' + filenameVmean + '.csv', delimiter=',', skip_header=1)
        meanDataX = np.genfromtxt('2/' + filenameXmean + '.csv', delimiter=',', skip_header=1)
        varDataV = np.genfromtxt('2/' + filenameVvar + '.csv', delimiter=',', skip_header=1)
        varDataX = np.genfromtxt('2/' + filenameXvar + '.csv', delimiter=',', skip_header=1)

        tV = meanDataV[:, 0]*1e3
        mV = meanDataV[:, 1]*1e3
        stdV = np.sqrt(varDataV[:, 1])*1e3
        tX = meanDataX[:, 0]*1e3
        mX = meanDataX[:, 1]*1e9
        stdX = np.sqrt(varDataX[:, 1]) * 1e9

        col = pl.cm.Blues([(7+shift)/8])
        if shift == 0:
            labX = '$\mu_x$'
            labV = '$\mu_v$'
        elif shift == -1:
            labX = '$\mu_x-\sigma_x$'
            labV = '$\mu_v-\sigma_v$'
        elif shift == 1:
            labX = '$\mu_x+\sigma_x$'
            labV = '$\mu_v+\sigma_v$'
        ax[0, i].plot(tX, mX+shift*stdX, color=col, label=labX)
        ax[0, i].set_xlabel('Time [ms]')
        ax[0, i].set_ylabel('Position [nm]')
        ax[1, i].plot(tV, mV+shift*stdV, color=col, label=labV)
        ax[1, i].set_xlabel('Time [ms]')
        ax[1, i].set_ylabel('Velocity [mm/s]')
        ax[0, i].set_title('$\\tau = $' + str(tau*1e6) + ' $\mu$s')
        i += 1
ax[0,0].legend()
ax[1,0].legend()
fig.savefig('2/VP_withmean.pdf')
plt.show()