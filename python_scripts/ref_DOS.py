import numpy as np
import matplotlib.pyplot as pp


t = 1
tp = 1
m = 1.5
N = 10**7
eps = 0.1
scale = (9)/(2 - eps)
Q = np.arccos(m/t-2)

def Energy(Kx, Ky, Kz):
    return np.sqrt(tp**2 * ((np.sin(Kx))**2 + np.sin(Ky)**2) + (t*(np.cos(Kx) + np.cos(Ky) + np.cos(Kz)) - m)**2)

def EngeryLeft(Kx, Ky, Kz):
    return np.sqrt((t*np.sin(Q)*(Kz - Q))**2 + tp**2*(Kx**2 + Ky**2))

def EngeryRight(Kx, Ky, Kz):
    return np.sqrt((t*np.sin(-Q)*(Kz + Q))**2 + tp**2*(Kx**2 + Ky**2))

def Single(Kx, Ky, Kz):
    return np.sqrt(Kx**2 + Ky**2 + Kz**2)

def SingleDOS(e):
    return 4*np.pi*e**2/((2*np.pi)**3)

def LowDOS(e):
    return np.pi*e**2/np.sin(Q)





Kx = (2*np.random.random(N) - 1)*2*np.pi
Ky = (2*np.random.random(N) - 1)*2*np.pi
Kz = (2*np.random.random(N) - 1)*2*np.pi

E = Energy(Kx, Ky, Kz)
E = np.append(E, -E)
E = 2*E/(E.max() - E.min())

Kx = (2*np.random.random(N) - 1)*2*np.pi
Ky = (2*np.random.random(N) - 1)*2*np.pi
Kz = (2*np.random.random(N) - 1)*2*np.pi

leftE = EngeryLeft(Kx, Ky, Kz)
leftE = np.append(leftE, -leftE)
RightE = EngeryRight(Kx, Ky, Kz)
RightE = np.append(RightE, -RightE)
LowE = np.append(leftE, RightE)
LowE = 2*LowE/(LowE.max() - LowE.min())


x = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_Plot_Points.csv", delimiter=',')[0:-1]
y = np.loadtxt("/Users/bruvisland/CLionProjects/KPM_4/data/DOS_Plot.csv", delimiter=',')[0:-1]


#pp.hist(E, bins=1000, density=True)

pp.plot(x[np.abs(x) < 0.25], LowDOS(x[np.abs(x) < 0.25]))
n,x2,_ = pp.hist(E,bins=1000, histtype=u'step',density=True )
bin_centers = 0.5*(x2[1:]+x2[:-1])
pp.plot(bin_centers, n) ## using bin_centers rather than edges
pp.hist(LowE,bins=1000,density=True )

pp.show()


singleE = Single(Kx, Ky, Kz)
singleE = np.append(singleE, -singleE)
singleE = singleE[np.abs(singleE) < 5]
pp.hist(singleE,bins=500,density=True )
x = np.linspace(singleE.min(), singleE.max(), 1000)
pp.plot(x, SingleDOS(x))
pp.show()