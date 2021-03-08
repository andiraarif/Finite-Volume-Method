import numpy as np
from matplotlib import pyplot as plt


# Input parameters
L = 1.0 # domain length (in m)
rho = 1.0 # fluid density in kg/m3
gamma = 0.1 # fluid diffusivity in kg/ms
phi0 = 1 # a scalar quantity (non-dimensional)
phi1 = 0 # a scalar quantity (non-dimensional)
U = 1 # fluid velocity in m/s
nx = 100 # number of grid

# Generate mesh
dx = L/nx

# Prepare the matrix
a = np.zeros((nx, nx))
b = np.zeros((nx, 1))

# Compute fluxes
F = rho*U
D = gamma/dx

for i in range(nx):
    if i == 0:
        aW = 0
        aE = D - F/2
        Sp = -(2*D + F)
        Su = -Sp*phi0
        a[i, i + 1] = -aE
    elif i > 0 and i < nx - 1:
        aW = D + F/2
        aE = D - F/2
        Sp = 0
        Su = 0
        a[i, i - 1] = -aW
        a[i, i + 1] = -aE
    elif i == nx - 1:
        aW = D + F/2
        aE = 0
        Sp = -(2*D - F)
        Su = -Sp*phi1
        a[i, i - 1] = -aW

    aP = aW + aE - Sp
    a[i, i] = aP
    b[i] = Su

# Solve the matrix
phi = np.linalg.solve(a, b)
phi = np.insert(phi, 0, phi0)
phi = np.append(phi,phi1)

# Post Processing
x = np.linspace(0 + dx/2, L - dx/2, nx)
x = np.insert(x, 0, 0)
x = np.append(x, L)
y = x
x_analytic =np.linspace(0, L, 50)

xx, yy = np.meshgrid(x, y)
phi_ = np.zeros((nx + 2, nx + 2))
for i in range(nx + 2):
    phi_[:, i] = phi[i]

#phi_analytic = [(2.7183 - np.exp(x_))/1.7183  for x_ in x_analytic] # u = 0.1 m/s
phi_analytic = [1 + (1 - np.exp(25*x_))/7.20e10  for x_ in x_analytic] # u = 2.5 m/s

#plt.plot(x, phi, '-sk')
#plt.plot(x_analytic, phi_analytic, '-k')
plt.contourf(xx, yy, phi_, 256, cmap = 'turbo')
plt.show()
