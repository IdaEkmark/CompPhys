#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

# skip_header skips the first
# row in data.csv
dt = 1e-3
filename = "QP_dt" + str(dt)
folder = "4/"
N = 256
axisfontsize = 14
nequi= 20000

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t = array[:, 0]
Q1x = array[:, 1]
Q1y = array[:, 2]
Q1z = array[:, 3]
Q2x = array[:, 7]
Q2y = array[:, 8]
Q2z = array[:, 9]
Q3x = array[:, 13]
Q3y = array[:, 14]
Q3z = array[:, 15]
Q4x = array[:, 19]
Q4y = array[:, 20]
Q4z = array[:, 21]
Q_list = [[Q1x, Q1y, Q1z], [Q2x, Q2y, Q2z], [Q3x, Q3y, Q3z], [Q4x, Q4y, Q4z]]
i = 1
for Q in Q_list:
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(5,7))
    ax1.plot(t, Q[0])
    ax1.set_xlabel('$t$ (ps)', fontsize=axisfontsize)
    ax1.set_ylabel('$x$ ($\mathrm{\AA}$)', fontsize=axisfontsize)
    ax1.grid()

    ax2.plot(t, Q[1])
    ax2.set_xlabel('$t$ (ps)', fontsize=axisfontsize)
    ax2.set_ylabel('$y$ ($\mathrm{\AA}$)', fontsize=axisfontsize)
    ax2.grid()

    ax3.plot(t, Q[2])
    ax3.set_xlabel('$t$ (ps)', fontsize=axisfontsize)
    ax3.set_ylabel('$z$ ($\mathrm{\AA}$)', fontsize=axisfontsize)
    ax3.grid()

    plt.suptitle('Position of particle ' + str(i), fontsize=axisfontsize)
    plt.tight_layout()
    plt.show()
    figurefilename = 'Position_particle_' + str(int(i))
    fig.savefig(folder + figurefilename + '_liquid.pdf')
    i = i + 85

filename_T = "T_dt" + str(dt)
filename_P = "P_dt" + str(dt)

array_T = np.genfromtxt(folder + filename_T + '.csv', delimiter=',', skip_header=1)
array_P = np.genfromtxt(folder + filename_P + '.csv', delimiter=',', skip_header=1)

t_T = array_T[:, 0]
T_inst = array_T[:, 1] - 273.15
t_P = array_P[:, 0]
P_inst = array_P[:, 1] / (6.24e-7)

print('T = ' + str(np.mean(T_inst[2*nequi:])))
print('P = ' + str(np.mean(P_inst[2*nequi:])))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
ax1.plot(np.array([t_T[nequi],t_T[nequi]]), np.array([0,1300]), '--', color = 'k')
ax1.plot(np.array([t_T[2*nequi],t_T[2*nequi]]), np.array([0,1300]), ':', color = 'k')
ax1.plot(t_T, T_inst)
ax1.set_ylabel('$T_\mathrm{instantaneous}$ ($^\mathrm{o}$C)', fontsize=axisfontsize)
ax1.set_xlabel('$t$ (ps)', fontsize=axisfontsize)
ax1.set_ylim([0,1300])
ax1.grid()
#ax1.set_xlim([-0.1,20])
ax2.plot(np.array([t_P[nequi],t_P[nequi]]), np.array([-3*10000, 2*10000]), '--', color = 'k')
ax2.plot(np.array([t_P[2*nequi],t_P[2*nequi]]), np.array([-3*10000, 2*10000]), ':', color = 'k')
ax2.plot(t_P, P_inst)
ax2.set_xlabel('$t$ (ps)', fontsize=axisfontsize)
ax2.set_ylabel('$P_\mathrm{instantaneous}$ (bar)', fontsize=axisfontsize)
ax2.set_ylim([-1.5*10000,2*10000])
ax2.grid()

ax2.legend(['Equilibration with $T=1000$ K', 'Equilibration with $T=700$ K'], loc='upper right', fontsize=axisfontsize)
#ax2.set_xlim([-0.1,20])
plt.suptitle('$T=$' + str(int(np.rint(np.mean(T_inst[2*nequi:])))) + ' $^\mathrm{o}$C and $P=$' +
             str(int(np.rint(np.mean(P_inst[2*nequi:])))) + ' bar during constant energy and volume simulation',
             fontsize=axisfontsize)
plt.tight_layout()
plt.show()
fig.savefig(folder + 'InstantaneousTemperaturePressure' + '_liquid.pdf')


filename = "a0_dt" + str(dt)

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

t1 = array[0:20002, 0]
t2 = array[20002:, 0]
a01 = array[0:20002, 1]
a02 = array[20002:, 1]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
ax1.plot(t1, a01)
ax1.set_ylabel('$a_0$ ($\mathrm{\AA}$)', fontsize=axisfontsize)
ax1.set_xlabel('$t$ (ps)', fontsize=axisfontsize)
ax1.set_title('Equilibration for $T=1000$ K', fontsize=axisfontsize)
ax1.grid()
ax2.plot(t2, a02)
ax2.set_ylabel('$a_0$ ($\mathrm{\AA}$)', fontsize=axisfontsize)
ax2.set_xlabel('$t$ (ps)', fontsize=axisfontsize)
ax2.set_title('Equilibration for $T=700$ K', fontsize=axisfontsize)
ax2.grid()
plt.suptitle('Lattice parameter during equilibration',  fontsize=axisfontsize)
plt.tight_layout()
plt.show()
fig.savefig(folder + 'latticeparam' + '_liquid.pdf')

