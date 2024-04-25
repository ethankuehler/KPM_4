import matplotlib.pyplot as pp
import scipy.integrate as integrate
import numpy as np

x = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_Plot_Points.csv", delimiter=',')
y = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_Plot.csv", delimiter=',')
up = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/UP_Plot.csv", delimiter=',')
down = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOWN_Plot.csv", delimiter=',')

local = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_local_Plot.csv", delimiter=',')
local_mag = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_local_mag_Plot.csv", delimiter=',')

localover = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_local_over_Plot.csv", delimiter=',')
local_over_mag = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_local_over_mag_Plot.csv", delimiter=',')
edge = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_edge.csv", delimiter=',')

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
pp.vlines(0, ymin=0, ymax=y.max(), colors='red')
pp.show()

pp.plot(x,local_mag)
pp.plot(x, local_over_mag)
pp.vlines(0, ymin=0, ymax=y.max(), colors='red')
pp.show()

pp.plot(x, y)
pp.plot(x, edge)
pp.plot(x, local)
pp.vlines(0, ymin=0, ymax=y.max(), colors='red')
pp.show()

