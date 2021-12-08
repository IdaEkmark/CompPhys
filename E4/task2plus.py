import numpy as np
import matplotlib.pyplot as plt

tau_list  = [48.5e-6, 147.3e-6]

def theoryPowSpec(fs, tau):
    kB = 1.380649e-23
    T = 297
    radius = 2.79e-6/2
    rho = 2.65e3
    m = (4.0/3.0 * np.pi * radius**3) * rho
    mplConst = kB*T/m
    f0 = 3.1e3
    eta = 1.0/tau

    powspec = mplConst * 2*eta*fs**2 / ( 4*np.pi**2 * (fs**2 - f0**2)**2 + (eta*fs)**2 )
    return powspec

fs = np.linspace(0,1e4,1000)
color = ['k']

fig, ax = plt.subplots(1, 2, figsize=(15,7))
for c in color:
    i=0
    for tau in tau_list:

        ax[i].plot(fs, theoryPowSpec(fs, tau), c)
        ax[i].set_xlabel('f [Hz]')
        ax[i].set_ylabel('Signal')
        ax[i].legend()
        ax[i].set_title('$\\tau = $' + str(tau*1e6) + ' $\mu$s')
        ax[i].set_xlim([0,10e3])
        i += 1
fig.savefig('2/powerTheory.pdf')
plt.show()