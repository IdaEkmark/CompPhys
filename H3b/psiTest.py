import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

# Position space
filename = "pDensPosCalc"
axisfontsize = 14
array = np.genfromtxt('1/' + filename + '.csv', delimiter=',', skip_header=1)

x = array[:, 0]
pDens = array[:, 1]

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(x, pDens)

ax.set_xlabel('$x$ (Å)', fontsize=axisfontsize)
ax.set_ylabel('Probability density', fontsize=axisfontsize)
ax.set_title('Probability density as a function of position x', fontsize=axisfontsize)
ax.grid()

plt.show()
fig.savefig('1/' + filename + '.pdf')
