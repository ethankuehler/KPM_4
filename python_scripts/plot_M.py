import matplotlib.pyplot as pp
import scipy.integrate as integrate
import numpy as np


x = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/Us_Plot.csv", delimiter=',')
avg = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/AVG_Ms_Plot.csv", delimiter=',')
rms = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/RMS_Ms_Plot.csv", delimiter=',')

pp.plot(x, avg)
pp.plot(x, rms)
pp.xlabel("Interaction strength U")
pp.ylabel("average M_z")
pp.title("Average M_z vs interaction strength")
pp.errorbar(x, avg, yerr=0.001)
pp.xticks(x)
pp.savefig("new_flip.png")
pp.show()