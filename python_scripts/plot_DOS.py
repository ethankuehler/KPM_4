import matplotlib.pyplot as pp
import scipy.integrate as integrate
import numpy as np

x = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_Plot_Points.csv", delimiter=',')[0:-1]
y = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_Plot.csv", delimiter=',')[0:-1]
up = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/UP_Plot.csv", delimiter=',')[0:-1]
down = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOWN_Plot.csv", delimiter=',')[0:-1]

local = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOSlocal_Plot.csv", delimiter=',')[0:-1]
localover = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOSlocalover_Plot.csv", delimiter=',')[0:-1]

print(x)
print(down.min())
print("N = ")
print(np.trapz(x, y*np.heaviside(-(x), 0.5)))

print("E = ")
print(np.trapz(x, y*np.heaviside(x, 0.5)*x))

nup = np.trapz(x, up*np.heaviside(-(x), 0.5))
print("<N_up> = {}".format(nup))
ndown = np.trapz(x, down*np.heaviside((-x), 0.5))
print("<N_down> = {}".format(ndown))
print("M = {}".format(nup-ndown))

pp.plot(x, y)
pp.plot(x, up)
pp.plot(x, down)
pp.vlines(0, ymin=0, ymax=y.max(), colors='red')

pp.show()

pp.plot(x, local)
pp.plot(x, localover)
pp.show()

