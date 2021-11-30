import numpy as np
import matplotlib.pyplot as plt

axisfontsize = 14
filename = "powerspectrum"
folder = "7/"

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

f = array[:, 0]
powerspectrum = array[:, 1]
fig, ax = plt.subplots(figsize=(7,5))
ax.plot(f, powerspectrum)
ax.set_xlabel('$\omega$ (rad/ps)', fontsize=axisfontsize)
ax.set_ylabel('$\Phi(\omega)$ ($\mathrm{\AA}^2$/ps$^2$)', fontsize=axisfontsize)
ax.set_title('$D_\mathrm{s}=$' + str(np.round(powerspectrum[0]/6, 3)) + ' $\mathrm{\AA}^2$/ps', fontsize=axisfontsize)
ax.grid()
#ax.set_xlim([-0.1, 2])
#ax.set_ylim([-5.5, 4])
plt.show()
fig.savefig(folder + filename + '.pdf')

