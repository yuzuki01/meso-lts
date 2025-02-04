import sys
import numpy as np


def print_help():
    print('''
Usage exampels:
    Using half-range Gasuu-Hermit quaderature points:
        setDV GH 28
    where 28 is the number of discrete velocities in each velocity space directions.
''')


def save_file(Xis, weights):
    print("Writing Xis and weights...\n")
    data = np.column_stack((Xis, weights))
    with open('../../build/hGH.dvs', 'w') as f_Xi:
        f_Xi.write("half-range-Gauss-Hermite\n")
        for xi, weight in data:
            f_Xi.write(f"{xi:18.15e} {weight:18.15e}\n")
    print("Writing Done!\n")


def generate_dv(N2):
    N = N2 // 2

    a = np.zeros(N)
    b = np.zeros(N)
    a[0] = 1.0 / np.sqrt(np.pi)
    a[1] = 2.0 / np.sqrt(np.pi) / (np.pi - 2.0)
    b[1] = a[0] / (a[0] + a[1]) / 2.0

    for i in range(2, N):
        b[i] = (i - 1) + 1.0 / 2.0 - b[i - 1] - a[i - 1] ** 2
        a[i] = (i ** 2 / 4.0 / b[i] - b[i - 1] - 1.0 / 2) / a[i - 1] - a[i - 1]

    J = np.diag(a) + np.diag(np.sqrt(b[1:N]), 1) \
        + np.diag(np.sqrt(b[1:N]), -1)

    v, V = np.linalg.eig(J)

    w = V[0, :] * V[0, :] * np.sqrt(np.pi) / 2.0

    vw = np.transpose(np.vstack((v, w)))
    vw = vw[vw[:, 0].argsort()]
    v = vw[:, 0]
    w = vw[:, 1]

    Xis = np.hstack((-np.flipud(v), v))
    weights = np.hstack((np.flipud(w), w))
    weights = weights * np.exp(Xis ** 2)
    Xis = Xis

    return Xis, weights


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print_help()
        exit()
    elif sys.argv[1] == 'GH':
        if int(sys.argv[2]) % 2 != 0 or int(sys.argv[2]) < 7:
            print("ERROR!  Number of discrete velocity should be even, and at least 6.")
            print_help()
            exit()

        Xis, weights = generate_dv(int(sys.argv[2]))
        save_file(Xis, weights)
    else:
        print_help()
        exit()
