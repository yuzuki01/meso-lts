from numpy import *

# molecular collision model constant
omega = 0.50
# initial density field, Kg/m^3
rho = 1.0
# Gas specific gas constant
R = 0.5
# Characteristic length, m
L = 1.0
# Knudsen number
Kn = 0.00018054066673528204
# Characteristic temperature, K
T = 1.0
# Inner freedom degree
k = 0
# Characteristic velocity, m/s
U = 0.1
# mean free path, m
mfp = L * Kn
# Most probable speed of molecular, m/2
C = sqrt(2 * R * T)
# viscosity
mu = (15 * rho * mfp) / (2 * (5 - 2 * omega) * (7 - 2 * omega)) \
     * sqrt(2 * pi * R * T)
# adiabatic constant
g = (5 + k) / (3 + k)
# Mach number
Ma = U / sqrt(g * R * T)
# Reynolds number
Re = rho * U * L / mu

print("Refer Density      (rho)  = ", rho)
print("Refer Length         (L)  = ", L)
print("Refer Temperature    (T)  = ", T)
print("Knudsen number      (Kn)  = ", Kn)
print("Mach number         (Ma)  = ", Ma)
print("Reynolds number     (Re)  = ", Re)
print("Viscosity           (mu)  = ", mu)
print("Spefic gas constant  (R)  = ", R)
print("Most probable speed  (C)  = ", C)
