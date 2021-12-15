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
T = s_P_corr[:,0]
s_P_corr = s_P_corr[:,1]
s_P_block = s_P_block[:, 5]
negInd = np.where(s_P_corr<0)[0]
s_P_corr[negInd] = N
negInd = np.where(s_P_block<0)[0]
s_P_block[negInd] = 32000

varP = np.genfromtxt('2/var_P.csv', delimiter=',', skip_header=1)
varP = varP[:, 1]

varP_corr = s_P_corr/N*varP
varP_block = s_P_block/N*varP

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,5))
ax1.plot(T, s_P_corr, label='Correlation')
ax1.plot(T, s_P_block, label='block')

ax2.plot(T, varP_corr, label='Correlation')
ax2.plot(T, varP_block, label='block')

ax1.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax1.set_ylabel('$n_s$', fontsize=axisfontsize)
ax1.set_title('Statistical inefficiency', fontsize=axisfontsize)
ax1.grid()
ax1.set_yscale('log')
ax1.legend(fontsize=axisfontsize)

ax2.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax2.set_ylabel('Var$(P)$', fontsize=axisfontsize)
ax2.set_title('Variance', fontsize=axisfontsize)
ax2.grid()
ax2.set_yscale('log')
ax2.legend(fontsize=axisfontsize)

plt.show()
fig.savefig('2/TP_se.pdf')

filename = "TPmean"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
P = array[:, 1]

skip=2
fig1, (ax1) = plt.subplots(1, 1,figsize=(12,5))
fig2, (ax2) = plt.subplots(1, 1,figsize=(12,5))
ax1.errorbar(T[::skip], P[::skip], yerr=3*np.sqrt(varP_corr)[::skip], capsize=2)
ax2.errorbar(T[::skip], P[::skip], yerr=3*np.sqrt(varP_block)[::skip], capsize=2)

ax1.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax1.set_ylabel('$P$', fontsize=axisfontsize)
ax1.set_title('Error evaluated from correlation', fontsize=axisfontsize)
ax1.legend(['$3\sigma_P$'], fontsize=axisfontsize)
ax1.grid()

ax2.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax2.set_ylabel('$P$', fontsize=axisfontsize)
ax2.set_title('Error evaluated from block averaging', fontsize=axisfontsize)
ax2.legend(['$3\sigma_P$'], fontsize=axisfontsize)
ax2.grid()

fig1.savefig('2/TP_error_corr.pdf')
fig2.savefig('2/TP_error_block.pdf')
plt.show()


s_r_corr = np.genfromtxt('2/s_r_corr.csv', delimiter=',', skip_header=1)
s_r_block = np.genfromtxt('2/s_r_block.csv', delimiter=',', skip_header=1)
T = s_r_corr[:,0]
s_r_corr = s_r_corr[:,1]
s_r_block = s_r_block[:, 5]
negInd = np.where(s_r_corr<0)[0]
s_r_corr[negInd] = N
negInd = np.where(s_r_block<0)[0]
s_r_block[negInd] = 32000

varr = np.genfromtxt('2/var_r.csv', delimiter=',', skip_header=1)
varr = varr[:, 1]

varr_corr = s_r_corr/N*varr
varr_block = s_r_block/N*varr

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,5))
ax1.plot(T, s_r_corr, label='Correlation')
ax1.plot(T, s_r_block, label='block')

ax2.plot(T, varr_corr, label='Correlation')
ax2.plot(T, varr_block, label='block')

ax1.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax1.set_ylabel('$n_s$', fontsize=axisfontsize)
ax1.set_title('Statistical inefficiency', fontsize=axisfontsize)
ax1.grid()
ax1.set_yscale('log')
ax1.legend(fontsize=axisfontsize)

ax2.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax2.set_ylabel('Var$(r)$', fontsize=axisfontsize)
ax2.set_title('Variance', fontsize=axisfontsize)
ax2.grid()
ax2.set_yscale('log')
ax2.legend(fontsize=axisfontsize)

plt.show()
fig.savefig('2/Tr_se.pdf')

filename = "Trmean"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
r = array[:, 1]

fig1, (ax1) = plt.subplots(1, 1,figsize=(12,5))
fig2, (ax2) = plt.subplots(1, 1,figsize=(12,5))
ax1.errorbar(T[::skip], r[::skip], yerr=3*np.sqrt(varr_corr)[::skip], capsize=2)
ax2.errorbar(T[::skip], r[::skip], yerr=3*np.sqrt(varr_block)[::skip], capsize=2)

ax1.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax1.set_ylabel('$r$', fontsize=axisfontsize)
ax1.set_title('Error evaluated from correlation', fontsize=axisfontsize)
ax1.legend(['$3\sigma_R$'], fontsize=axisfontsize)
ax1.grid()

ax2.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax2.set_ylabel('$r$', fontsize=axisfontsize)
ax2.set_title('Error evaluated from block averaging', fontsize=axisfontsize)
ax2.legend(['$3\sigma_r$'], fontsize=axisfontsize)
ax2.grid()

fig1.savefig('2/Tr_error_corr.pdf')
fig2.savefig('2/Tr_error_block.pdf')
plt.show()



s_Em_corr = np.genfromtxt('2/s_E_m_corr.csv', delimiter=',', skip_header=1)
s_Em_block = np.genfromtxt('2/s_E_m_block.csv', delimiter=',', skip_header=1)
T = s_Em_corr[:,0]
s_Em_corr = s_Em_corr[:,1]
s_Em_block = s_Em_block[:, 5]
negInd = np.where(s_Em_corr<0)[0]
s_Em_corr[negInd] = N
negInd = np.where(s_Em_block<0)[0]
s_Em_block[negInd] = 32000

varEm = np.genfromtxt('2/var_E_m.csv', delimiter=',', skip_header=1)
varEm = varEm[:, 1]

varU_corr = s_Em_corr/N*varEm
varU_block = s_Em_block/N*varEm

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,5))
ax1.plot(T, s_Em_corr, label='Correlation')
ax1.plot(T, s_Em_block, label='block')

ax2.plot(T, varU_corr, label='Correlation')
ax2.plot(T, varU_block, label='block')

ax1.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax1.set_ylabel('$n_s$', fontsize=axisfontsize)
ax1.set_title('Statistical inefficiency', fontsize=axisfontsize)
ax1.grid()
ax1.set_yscale('log')
ax1.legend(fontsize=axisfontsize)

ax2.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax2.set_ylabel('Var$(U)$', fontsize=axisfontsize)
ax2.set_title('Variance', fontsize=axisfontsize)
ax2.grid()
ax2.set_yscale('log')
ax2.legend(fontsize=axisfontsize)

plt.show()
fig.savefig('2/TEm_se.pdf')

filename = "TU"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

T = array[:, 0]
U = array[:, 1]

fig1, (ax1) = plt.subplots(1, 1,figsize=(12,5))
fig2, (ax2) = plt.subplots(1, 1,figsize=(12,5))
ax1.errorbar(T[::skip], U[::skip], yerr=3*np.sqrt(varU_corr)[::skip], capsize=2)
ax2.errorbar(T[::skip], U[::skip], yerr=3*np.sqrt(varU_block)[::skip], capsize=2)

ax1.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax1.set_ylabel('$U$ (eV)', fontsize=axisfontsize)
ax1.set_title('Error evaluated from correlation', fontsize=axisfontsize)
ax1.legend(['$3\sigma_U$'], fontsize=axisfontsize)
ax1.grid()

ax2.set_xlabel('$T$ (K)', fontsize=axisfontsize)
ax2.set_ylabel('$U$ (eV)', fontsize=axisfontsize)
ax2.set_title('Error evaluated from block averaging', fontsize=axisfontsize)
ax2.legend(['$3\sigma_U$'], fontsize=axisfontsize)
ax2.grid()

fig1.savefig('2/TU_error_corr.pdf')
fig2.savefig('2/TU_error_block.pdf')
plt.show()
