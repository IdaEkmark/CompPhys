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

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
dt = 1e-3
filename = "QP_dt" + str(dt)
folder = "4/"
N = 256

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
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(15,7))
    ax1.plot(t, Q[0])
    ax1.set_xlabel('$t$ (ps)')
    ax1.set_ylabel('$x$ (A)')
    ax1.grid()

    ax2.plot(t, Q[1])
    ax2.set_xlabel('$t$ (ps)')
    ax2.set_ylabel('$y$ (A)')
    ax2.grid()

    ax3.plot(t, Q[2])
    ax3.set_xlabel('$t$ (ps)')
    ax3.set_ylabel('$z$ (A)')
    ax3.grid()

    plt.suptitle('Position, $dt=$' + str(dt) + ', particle ' + str(i))
    plt.show()
    figurefilename = 'Position_particle_' + str(int(i))
    fig.savefig(folder + figurefilename + '.pdf')
    i = i + 85

filename_T = "T_dt" + str(dt)
filename_P = "P_dt" + str(dt)

array_T = np.genfromtxt(folder + filename_T + '.csv', delimiter=',', skip_header=1)
array_P = np.genfromtxt(folder + filename_P + '.csv', delimiter=',', skip_header=1)

t_T = array_T[:, 1]
T_inst   = array_T[:, 0] - 273.15
t_P = array_P[:, 1]
P_inst   = array_P[:, 0] / (6.24e-7)

print('T = ' + str(np.mean(T_inst)))
print('P = ' + str(np.mean(P_inst)))
