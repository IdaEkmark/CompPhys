import numpy as np
import matplotlib.pyplot as plt

N_vec = [int(10), int(100), int(1000), int(10000)]
I_vec = [0.1489, 0.1807, 0.1665, 0.1666]
for N, I in zip(N_vec, I_vec):
    filename = "f_x_i_" + str(N)
    folder = "1/"

    arr = np.genfromtxt(folder + filename + '.csv', delimiter=',', skip_header=1)

    x_i = arr[:, 0]
    f_i = arr[:, 1]

    f_var = np.mean(f_i**2) - np.mean(f_i)**2
    err = np.sqrt(f_var / (float(N)-1))



    plt.hist(x_i, bins = 'auto')
    plt.title('I=' + str(I) + ', N=' + str(N) + ', error is ' + str(np.round(err, 5)))
    plt.show()

