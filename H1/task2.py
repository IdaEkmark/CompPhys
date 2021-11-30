#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
folder = "2/"
N = 256
d = .01
axisfontsize = 14

# skip_header skips the first
# row in data.csv
dt = 1e-3
filename = "E_dt" + str(dt)
filename_p = "Epot_dt" + str(dt)
filename_k = "Ekin_dt" + str(dt)


array_p = np.genfromtxt(folder + filename_p + '.csv', delimiter=',', skip_header=1)
array_k = np.genfromtxt(folder + filename_k + '.csv', delimiter=',', skip_header=1)

t_p = array_p[:, 0]
E_p = array_p[:, 1]
t_k = array_k[:, 0]
E_k = array_k[:, 1]
T = 2*np.mean(E_k)/(3*N) / (8.6e-5)

fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 5))
ax.plot(t_k, E_k+E_p, label = '$E_\mathrm{total}$')
ax.plot(t_p, E_p, label = '$E_\mathrm{potential}$')
ax.plot(t_k, E_k, label = '$E_\mathrm{kinetic}$')
ax2.plot(t_k, E_k+E_p, label = '$E_\mathrm{total}$')
ax2.plot(t_p, E_p, label = '$E_\mathrm{potential}$')

ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax.set_ylim([0, 60])
ax.set_yticks([10, 30, 50])
ax.set_xlim([0, 2.5])
ax2.set_ylim([-860, -800])
ax2.set_yticks([-810, -830, -850])
ax2.set_xlim([0, 2.5])
ax2.set_xticks([0, 0.5, 1, 1.5, 2, 2.5])

# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

fig.supxlabel('$t$ (ps)', fontsize=axisfontsize)
fig.supylabel('$E$ (eV)', fontsize=axisfontsize)
ax.grid()
ax2.grid()
ax.legend(loc='upper right', fontsize=axisfontsize)
plt.suptitle('$T=$' + str(int(round(T,0))) + ' K, $dt=$' + str(dt) + ' ps', fontsize=axisfontsize)
plt.tight_layout()
fig.savefig(folder + filename + '.pdf')
plt.close()
#'''

dt = 1e-2
filename = "E_dt" + str(dt)
filename_p = "Epot_dt" + str(dt)
filename_k = "Ekin_dt" + str(dt)

array_p = np.genfromtxt(folder + filename_p + '.csv', delimiter=',', skip_header=1)
array_k = np.genfromtxt(folder + filename_k + '.csv', delimiter=',', skip_header=1)

t_p = array_p[:, 0]
E_p = array_p[:, 1]
t_k = array_k[:, 0]
E_k = array_k[:, 1]
T = 2*np.mean(E_k)/(3*N) / (8.6e-5)

fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 5))
ax.plot(t_k, E_k+E_p, label = '$E_\mathrm{total}$')
ax.plot(t_p, E_p, label = '$E_\mathrm{potential}$')
ax.plot(t_k, E_k, label = '$E_\mathrm{kinetic}$')
ax2.plot(t_k, E_k+E_p, label = '$E_\mathrm{total}$')
ax2.plot(t_p, E_p, label = '$E_\mathrm{potential}$')

ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax.set_ylim([0, 60])
ax.set_xlim([0, 2.5])
ax.set_yticks([10, 30, 50])
ax2.set_ylim([-860, -800])
ax2.set_xlim([0, 2.5])
ax2.set_yticks([-810, -830, -850])
ax2.set_xticks([0, 0.5, 1, 1.5, 2, 2.5])

# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

fig.supxlabel('$t$ (ps)', fontsize=axisfontsize)
fig.supylabel('$E$ (eV)', fontsize=axisfontsize)
ax.grid()
ax2.grid()
ax.legend(loc='upper right', fontsize=axisfontsize)
plt.suptitle('$T=$' + str(int(round(T,0))) + ' K, $dt=$' + str(dt) + ' ps', fontsize=axisfontsize)

plt.tight_layout()
fig.savefig(folder + filename + '.pdf')
plt.close()

#'''

dt = 2e-2
filename = "E_dt" + str(dt)
filename_p = "Epot_dt" + str(dt)
filename_k = "Ekin_dt" + str(dt)

array_p = np.genfromtxt(folder + filename_p + '.csv', delimiter=',', skip_header=1)
array_k = np.genfromtxt(folder + filename_k + '.csv', delimiter=',', skip_header=1)

t_p = array_p[:, 0]
E_p = array_p[:, 1]
t_k = array_k[:, 0]
E_k = array_k[:, 1]
T = 2*np.mean(E_k)/(3*N) / (8.6e-5)

fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 5))
ax.plot(t_k, E_k+E_p, label = '$E_\mathrm{total}$')
ax.plot(t_p, E_p, label = '$E_\mathrm{potential}$')
ax.plot(t_k, E_k, label = '$E_\mathrm{kinetic}$')
ax2.plot(t_k, E_k+E_p, label = '$E_\mathrm{total}$')
ax2.plot(t_p, E_p, label = '$E_\mathrm{potential}$')

ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax.set_ylim([0, 100])
ax.set_xlim([0, 2.5])
ax.set_yticks([20, 60, 100])
ax2.set_ylim([-860, -760])
ax2.set_xlim([0, 2.5])
ax2.set_yticks([-780, -820, -860])
ax2.set_xticks([0, 0.5, 1, 1.5, 2, 2.5])

# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

fig.supxlabel('$t$ (ps)', fontsize=axisfontsize)
fig.supylabel('$E$ (eV)', fontsize=axisfontsize)
ax.grid()
ax2.grid()
ax.legend(loc='upper left', fontsize=axisfontsize)
plt.suptitle('$T=$' + str((round(T/1e9,1))) + '$\cdot 10^{9}$ K, $dt=$' + str(dt) + ' ps', fontsize=axisfontsize)

plt.tight_layout()
fig.savefig(folder + filename + '.pdf')
plt.close()
#'''