import numpy as np
import matplotlib.pyplot as pp

x = np.loadtxt("../data/mu_dos.csv", delimiter=",")
pp.plot(x)
pp.show()
