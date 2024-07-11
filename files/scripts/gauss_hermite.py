import numpy as np

n = 3
x, w = np.polynomial.hermite.hermgauss(n)

with open(f"gh-{n}.dvs", "w") as f:
    f.write("Gauss-Hermite\n")
    for i, xi in enumerate(x):
        f.write(f"{xi:.18f} {w[i]:.18f}\n")
