import matplotlib.pyplot as pp
import scipy.integrate as integrate
import numpy as np


x = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/Us_Plot.csv", delimiter=',')[0:-1]
y = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/Ms_Plot.csv", delimiter=',')[0:-1]

pp.plot(x, y)
pp.xlabel("Interaction strength U")
pp.ylabel("M_z")
pp.errorbar(x, y, yerr=0.0001)
pp.show()