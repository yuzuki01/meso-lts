from numpy import *

# molecular collision model constant
alpha = 1.0
# omega_HS = 0.5
# omega_VHS = 0.81
omega = 0.74
# Boltzmann constant
k = 1.38066e-23
# Molecular mass, kg
m = 4.65e-26
# Molecular hard sphere diameter, m
d = 4.17e-10
# Gas specific gas constant
R = k / m
# Characteristic length, m
L = 20.0e-06
# Knudsen number
Kn = 0.72
# Characteristic temperature, K
T = 273.0
# dimensionless viscosity
muBar = 5.0 * (alpha + 1.0) * (alpha + 2.0) * sqrt(pi) * Kn \
        / (4.0 * alpha * (5.0 - 2 * omega) * (7.0 - 2.0 * omega))
# mean free path, m
mfp = L * Kn
# initial density field, Kg/m^3
rho = m / (mfp * sqrt(2.0) * pi * d * d)
# Most probable speed of molecular, m/2
C = sqrt(2 * R * T)
# viscosity
mu = rho * C * L * muBar
# Adiabatic constant
g = 1.4

print("Refer Temperature  (T )  = ", T)
print("Knudsen number     (Kn)  = ", Kn)
print("Viscosity          (mu)  = ", mu)
print("Density           (rho)  = ", rho)
print("Spefic gas constant (R)  = ", R)
print("Most probable speed (C)  = ", C)
