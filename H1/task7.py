import numpy as np
import matplotlib.pyplot as plt

filename = "powerspectrum"
folder = "7/"

array = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

f = array[:, 0]
powerspectrum = array[:, 1]
fig, ax = plt.subplots(figsize=(15,7))
ax.plot(f, powerspectrum)
ax.set_xlabel('$\omega$ (ps$^{-1}$)')
ax.set_title('$D_s=$' + str(np.round(powerspectrum[0]/6, 3)))
ax.grid()
plt.show()
fig.savefig(folder + filename + '.pdf')

