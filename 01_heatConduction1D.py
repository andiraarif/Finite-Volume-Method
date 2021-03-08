import numpy as np
from matplotlib import pyplot as plt


# Input parameters
L = 0.5 # in m
A = 1e-2 # in m2
k = 1000
TA = 100 # in deg celcius
TB = 500 # in deg celcius
nx = 5

# Generate mesh
dx = L/nx

# Prepare the matrix
a = np.zeros((nx, nx))
b = np.zeros((nx, 1))

# Start calculation
for i in range(nx):
    if i == 0:
        aW = 0
        aE = k*A/dx
        Sp = -2*k*A/dx
        Su = -Sp*TA
        a[i, i + 1] = -aE
    elif i > 0 and i < nx - 1:
        aW = k*A/dx
        aE = k*A/dx
        Sp = 0
        Su = 0
        a[i, i - 1] = -aW
        a[i, i + 1] = -aE
    elif i == nx - 1:
        aW = k*A/dx
        aE = 0
        Sp = -2*k*A/dx
        Su = -Sp*TB
        a[i, i - 1] = -aW

    aP = aW + aE - Sp
    a[i, i] = aP
    b[i] = Su

# Solve the matrix
T = np.linalg.solve(a, b)
T = np.insert(T, 0, TA)
T = np.append(T, TB)

# Post Processing
x = np.linspace(0 + dx/2, L - dx/2, nx)
x = np.insert(x, 0, 0)
x = np.append(x, L)
y = x

xx, yy = np.meshgrid(x, y)
TT = np.zeros((nx + 2, nx + 2))
for i in range(nx + 2):
    TT[:, i] = T[i]

T_analytic = [800*x_ + 100 for x_ in x]

fig, ax = plt.subplots()
#plt.plot(x, T, 'sk')
#plt.plot(x, T_analytic, '-k')
ax.set_aspect(1e-1/0.5)
plt.contourf(xx, yy, TT, 256, cmap = 'turbo')
plt.show()

print(a)
print(b)
print(T)
print(x)
