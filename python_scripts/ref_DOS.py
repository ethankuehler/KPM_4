import numpy as np
import matplotlib.pyplot as pp

def Energy(Kx, Ky, Kz):
    return np.sqrt(tp**2 * ((np.sin(Kx))**2 + np.cos(Ky)**2) + (t*(np.cos(Kx) + np.cos(Ky) + np.cos(Kz)) - m)**2)


t = 1
tp = 1
m = 0.5
N = 10**8


Kx = (2*np.random.random(N) - 1)*2*np.pi
Ky = (2*np.random.random(N) - 1)*2*np.pi
Kz = (2*np.random.random(N) - 1)*2*np.pi

E = Energy(Kx, Ky, Kz)
K = np.sqrt(Kx**2 + Ky**2 + Kz**2)

#print(np.array([Kx, Ky, Kz]))
#print(Energy(Kx[0], Ky[0], Kz[0]))
#print(E)

print(np.sort(E)[0])

pp.hist(E, bins=1000, density=True)
pp.show()