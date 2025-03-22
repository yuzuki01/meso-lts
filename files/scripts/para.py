from numpy import *

# molecular collision model constant
alpha = 1.0
# omega_HS = 0.5
# omega_VHS = 0.81
omega = 0.74
# Boltzmann constant
k = 1.38066e-23
# Molecular mass, kg
m = 46.5e-27
# Molecular hard sphere diameter, m
d = 4.17e-10
# Gas specific gas constant
R = k / m
# Characteristic length, m
L = 4.0e-7
# Knudsen number
Kn = 0.1
# Characteristic temperature, K
T = 273.0
# mean free path, m
mfp = L * Kn
# initial density field, Kg/m^3
rho = m / (mfp * sqrt(2.0) * pi * d * d)
# Most probable speed of molecular, m/2
C = sqrt(2 * R * T)
# viscosity
mu = (15 * rho * mfp) / (2 * (5 - 2 * omega) * (7 - 2 * omega)) \
     * sqrt(2 * pi * R * T)
# adiabatic constant
g = (5 + k) / (3 + k)

print("Refer Density      (rho)  = ", rho)
print("Refer Length         (L)  = ", L)
print("Refer Temperature    (T)  = ", T)
print("Knudsen number      (Kn)  = ", Kn)
print("Viscosity           (mu)  = ", mu)
print("Spefic gas constant  (R)  = ", R)
print("Most probable speed  (C)  = ", C)
