from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def make_stiffness_matrix(n, mod_elas, mom_inertia, length):
    k_matrix = np.array([
        [
            12*mod_elas*mom_inertia/(length/n)**3,
            6*mod_elas*mom_inertia/(length/n)**2,
            -12*mod_elas*mom_inertia/(length/n)**3,
            6*mod_elas*mom_inertia/(length/n)**2
        ],
        [
            6*mod_elas*mom_inertia/(length/n)**2,
            4*mod_elas*mom_inertia/(length/n),
            -6*mod_elas*mom_inertia/(length/n)**2,
            2*mod_elas*mom_inertia/(length/n)
        ],
        [
            -12*mod_elas*mom_inertia/(length/n)**3,
            -6*mod_elas*mom_inertia/(length/n)**2,
            12*mod_elas*mom_inertia/(length/n)**3,
            -6*mod_elas*mom_inertia/(length/n)**2
        ],
        [
            6*mod_elas*mom_inertia/(length/n)**2,
            2*mod_elas*mom_inertia/(length/n),
            -6*mod_elas*mom_inertia/(length/n)**2,
            4*mod_elas*mom_inertia/(length/n)
        ]
    ])

    stiffness_matrix = np.zeros((2*n + 2, 2*n + 2), dtype='float')
    i, j = 0, 0
    counti, countj = 0, 0

    while i < 2*n + 2:
        while j < 2*n + 2:
            stiffness_matrix[i][j] += k_matrix[counti][countj]
            j += 1
            countj += 1
            if countj == 4:
                if counti != 3:
                    i += 1
                    countj = 0
                    j -= 4
                    counti += 1
                    continue
                if counti == 3:
                    if i == 2*n + 1:
                        break
                    j -= 2
                    i -= 1
                    counti, countj = 0, 0
        break
    return stiffness_matrix


def make_mass_matrix(n, area, density, length):
    m_matrix = np.array([
        [
            13*area*(length/n)*density/35,
            11*area*((length/n)**2)*density/210,
            9*area*(length/n)*density/70,
            -13*area*((length/n)**2)*density/420
        ],
        [
            11*area*((length/n)**2)*density/210,
            area*((length/n)**3)*density/105,
            13*area*((length/n)**2)*density/420,
            -area*((length/n)**3)*density/140
        ],
        [
            9*area*(length/n)*density/70,
            13*area*((length/n)**2)*density/420,
            13*area*(length/n)*density/35,
            -11*area*((length/n)**2)*density/210
        ],
        [
            -13*area*((length/n)**2)*density/420,
            -area*((length/n)**3)*density/140,
            -11*area*((length/n)**2)*density/210,
            area*((length/n)**3)*density/105
        ]
    ])

    mass_matrix = np.zeros((2*n + 2, 2*n + 2), dtype='float')
    i, j = 0, 0
    counti, countj = 0, 0

    while i < 2*n + 2:
        while j < 2*n + 2:
            mass_matrix[i][j] += m_matrix[counti][countj]
            j += 1
            countj += 1
            if countj == 4:
                if counti != 3:
                    i += 1
                    countj = 0
                    j -= 4
                    counti += 1
                    continue
                if counti == 3:
                    if i == 2*n + 1:
                        break
                    j -= 2
                    i -= 1
                    counti, countj = 0, 0
        break
    return mass_matrix


def main():

    number = 10  # number of elements taken

    global_stiffness_matrix = make_stiffness_matrix(number, )
    # make_csv(STIFFNESS_MATRIX, "K")

    global_mass_matrix = make_mass_matrix(number, )
    # make_csv(MASS_MATRIX, "M")

    STIFFNESS_MATRIX = np.delete(STIFFNESS_MATRIX, np.s_[0:2], 0)
    STIFFNESS_MATRIX = np.delete(STIFFNESS_MATRIX, np.s_[0:2], 1)

    MASS_MATRIX = np.delete(MASS_MATRIX, np.s_[0:2], 0)
    MASS_MATRIX = np.delete(MASS_MATRIX, np.s_[0:2], 1)

    evals, evecs = la.eigh(STIFFNESS_MATRIX, MASS_MATRIX)

    # make_csv(evals, "evals")
    # make_csv(evecs, "evecs")

    frequencies = np.sqrt(evals) / (2 * np.pi)
    print(frequencies)


if __name__ == 'main':
    main()
