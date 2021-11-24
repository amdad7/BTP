from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


# def hit_spring(t, y):
#     return y[0] + 1
#
#
# hit_spring.terminal = True


def F1(t, y, m1, m2, k1, k2, c1, c2, f1, f2):
##########################################################################################################################################################################################
    m1 = 3
    m2 = 5

    k1 = 7
    k2 = 9

    c1 = 1
    c2 = 2

    f1 = 40*np.cos(3*t)
    f2 = 0
##########################################################################################################################################################################################

    M = np.array([
        [m1, 0],
        [0, m2]
    ])

    C = np.array([
        [c1 + c2, -c2],
        [-c2, c2]
    ])

    K = np.array([
        [k1 + k2, -k2],
        [-k2, k2]
    ])

    KK = np.array([
        [2*k1 + k2, -k2],
        [-k2, k2]
    ])

    A = np.vstack([
        np.hstack([
            np.zeros((2, 2)),
            np.eye(2, 2)
        ]),
        np.hstack([
            -np.matmul(np.linalg.inv(M), K),
            -np.matmul(np.linalg.inv(M), C)
        ])
    ])

    AA = np.vstack([
        np.hstack([
            np.zeros((2, 2)),
            np.eye(2, 2)
        ]),
        np.hstack([
            -np.matmul(np.linalg.inv(M), KK),
            -np.matmul(np.linalg.inv(M), C)
        ])
    ])


    B = np.vstack([
        np.zeros((2, 2)),
        np.linalg.inv(M)
    ])

    F = np.array([
        [f1],
        [f2]
    ])

    yvec = np.array([[y[i] for i in range(4)]]).T

    ydot1 = np.matmul(A, yvec) + np.matmul(B, F)    # k1 = 7
    ydot2 = np.matmul(AA, yvec) + np.matmul(B, F)   # k1 = 14
    print(y[0])

    return ydot1


def main():
    pass


if __name__ == 'main':
    main()
