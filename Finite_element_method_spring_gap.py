from scipy.integrate import solve_ivp
import scipy.linalg as la
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
    number = 10                     # number of elements taken
    length_of_beam = 5              # in meters
    breadth_of_beam = 0.05          # in meters
    height_of_beam = 0.05           # in meters
    modulus_of_elasticity = 7E10
    density_of_beam = 2700  

    cross_sect_area = breadth_of_beam*height_of_beam

    area_mom_inertia = (breadth_of_beam*height_of_beam**3)/12

    global_stiffness_matrix = make_stiffness_matrix(number, modulus_of_elasticity, area_mom_inertia, length_of_beam)

    global_mass_matrix = make_mass_matrix(number, cross_sect_area, density_of_beam, length_of_beam)

    global_stiffness_matrix = np.delete(global_stiffness_matrix, np.s_[0:2], 0)
    global_stiffness_matrix = np.delete(global_stiffness_matrix, np.s_[0:2], 1)

    global_mass_matrix = np.delete(global_mass_matrix, np.s_[0:2], 0)
    global_mass_matrix = np.delete(global_mass_matrix, np.s_[0:2], 1)

    evals, evecs = la.eigh(global_stiffness_matrix, global_mass_matrix)

    frequencies = np.sqrt(evals)/(2*np.pi)
    print(frequencies)

    return 



if __name__ == '__main__':
    main()
