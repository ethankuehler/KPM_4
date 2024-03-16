import numpy as np
import matplotlib.pyplot as pp
import scipy.integrate as integrate

DOS = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS.csv")
mag = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/mag.csv")
x = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/points.csv")

def Z(B):
    return integrate.trapz(DOS*np.exp(-x/B), x=x)

def A(B):
    return integrate.trapz(mag*np.exp(-x/B)/Z(B), x=x)


B = np.linspace(0.01, 10, 100)

DOS[DOS==np.inf] = 0
DOS[DOS==-np.inf] = 0
mag[mag==np.inf] = 0
mag[mag==-np.inf] = 0

f = np.vectorize(A)
print(f(B))
print(mag.min())
pp.plot(B, f(B))
pp.show()
