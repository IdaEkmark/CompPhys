#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
axisfontsize = 14
N = 4e6

s_P_corr = np.genfromtxt('2/s_P_corr.csv', delimiter=',', skip_header=1)
s_P_block = np.genfromtxt('2/s_P_block.csv', delimiter=',', skip_header=1)
varP = np.genfromtxt('2/var_P.csv', delimiter=',', skip_header=1)

varP_corr = s_P_corr[:, 1]/N*varP[:, 1]
varP_block = s_P_block[:, 1]/N*varP[:, 1]

filename = "TPmean"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
P = array[:, 1]

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,5))
ax1.errorbar(T, P, yerr=np.sqrt(varP_corr))
ax2.errorbar(T, P, yerr=np.sqrt(varP_block))

ax1.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax1.set_ylabel('$P$', fontsize=axisfontsize)
ax1.set_title('Error evaluated from correlation')
ax1.grid()

ax2.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax2.set_ylabel('$P$', fontsize=axisfontsize)
ax2.set_title('Error evaluated from block averaging')
ax2.grid()

plt.show()
fig.savefig('2/TP_error.pdf')


s_r_corr = np.genfromtxt('2/s_r_corr.csv', delimiter=',', skip_header=1)
s_r_block = np.genfromtxt('2/s_r_block.csv', delimiter=',', skip_header=1)
varr = np.genfromtxt('2/var_r.csv', delimiter=',', skip_header=1)

varr_corr = s_r_corr[:, 1]/N*varr[:, 1]
varr_block = s_r_block[:, 1]/N*varr[:, 1]

filename = "Trmean"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
r = array[:, 1]

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,5))
ax1.errorbar(T, r, yerr=np.sqrt(varr_corr))
ax2.errorbar(T, r, yerr=np.sqrt(varr_block))

ax1.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax1.set_ylabel('$r$', fontsize=axisfontsize)
ax1.set_title('Error evaluated from correlation')
ax1.grid()

ax2.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax2.set_ylabel('$r$', fontsize=axisfontsize)
ax2.set_title('Error evaluated from block averaging')
ax2.grid()

plt.show()
fig.savefig('2/Tr_error.pdf')


s_Em_corr = np.genfromtxt('2/s_E_m_corr.csv', delimiter=',', skip_header=1)
s_Em_block = np.genfromtxt('2/s_E_m_block.csv', delimiter=',', skip_header=1)
varEm = np.genfromtxt('2/var_E_m.csv', delimiter=',', skip_header=1)

varU_corr = s_Em_corr[:, 1]/N*varEm[:, 1]
varU_block = s_Em_block[:, 1]/N*varEm[:, 1]

negInd = np.where(varU_corr<0)[0]
varU_corr[negInd] = 0

filename = "TU"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
U = array[:, 1]

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,5))
ax1.errorbar(T, U, yerr=np.sqrt(varU_corr))
ax2.errorbar(T, U, yerr=np.sqrt(varU_block))

ax1.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax1.set_ylabel('$U$ (eV)', fontsize=axisfontsize)
ax1.set_title('Error evaluated from correlation')
ax1.grid()

ax2.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax2.set_ylabel('$U$ (eV)', fontsize=axisfontsize)
ax2.set_title('Error evaluated from block averaging')
ax2.grid()

plt.show()
fig.savefig('2/TU_error.pdf')
