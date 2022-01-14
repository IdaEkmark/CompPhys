#!/usr/bin/env python
###############################################################################
# E1code2 (slightly modified)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace

# skip_header skips the first
# row in data.csv

HBAR = 6.582119569e-1
HMASS = 104.453702

# Position space
x = np.zeros((6,401))
pDens = np.zeros((6,401))

xAvgs = np.zeros(6)
x2Avgs = np.zeros(6)
ts = np.arange(0,51,10)

for i in range(6):
    filename = "pDensPos" + str(i) + "0fs"
    array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)

    x[i,:] = array[:, 0]
    pDens[i,:] = array[:, 1]

    dx = x[i,1]-x[i,0]

    xList = [x[i,j]*pDens[i,j] for j in range(401)]
    xAvgs[i] = sum(xList)*dx
    #print("xAvg = " + str(xAvgs[i]))

    x2List = [x[i,j]**2 * pDens[i,j] for j in range(401)]
    x2Avgs[i] = sum(x2List)*dx
    #print("x2Avg = " + str(x2Avgs[i]))

xVars = x2Avgs - np.square(xAvgs)

# Momentum space
filename = "pDensMom00fs"
array = np.genfromtxt('2/' + filename + '.csv', delimiter=',', skip_header=1)
ps = array[:, 0]
pDensMom = array[:, 1]

dp = ps[1]-ps[0]

pList = [ps[j]*pDensMom[j] for j in range(401)]
pAvg = sum(pList)*dp
#print("xAvg = " + str(pAvgs[i]))

p2List = [ps[j]**2 * pDensMom[j] for j in range(401)]
p2Avg = sum(p2List)*dp
#print("x2Avg = " + str(p2Avgs[i]))

pVar = p2Avg - np.square(pAvg)
print("--------")
print("Momentum uncertainty is: " + str(pVar))
print("Momentum uncertainty should be: " + str(HBAR**2/(4*xVars[0])))
print("--------")

tLin = linspace(-5,55)

def sigSq(sigSq0, sigSqMom, t):
    return sigSq0 + sigSqMom/(HMASS**2) * t**2

# Plotting
axisfontsize = 14
fig, ax = plt.subplots(figsize=(9,5))
ax.plot(tLin, sigSq(xVars[0], pVar, tLin), '--', label=("$\sigma_x^2(t) = \sigma_0^2 + \\frac{\sigma_p^2}{m^2} t^2$, Theory"))
ax.plot(ts, xVars, 'ks', label=("$\sigma_x^2(t)$, Data"))

ax.set_xlabel('$t$ (fs)', fontsize=axisfontsize)
ax.set_ylabel('$\sigma_x^2$ (Å$^2$)', fontsize=axisfontsize)
ax.set_title('Variance in $x$ as a function of time', fontsize=axisfontsize)
ax.legend()

ax.grid()

plt.show()
fig.savefig('2/' + "sigSq" + '.pdf')