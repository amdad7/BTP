from scipy.integrate import solve_ivp
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt


def plot_graph(x_value, y_values, x_axis_label, y_axis_label, graph_title, fig_length=20, fig_height=9, gap=None):
    plt.figure(figsize=(fig_length, fig_height))
    plt.plot(x_value, y_values, linewidth=0.8, color = 'black')
    if gap is not None:
        plt.plot(x_value, np.full((len(y_values)), gap))
    plt.title(graph_title)
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    plt.show()


def make_stiffness_matrix(n, length, breadth, height, mod_elas):
    
    mom_inertia = (breadth*height**3)/12
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
    
    stiffness_matrix = np.delete(stiffness_matrix, np.s_[0:2], 0)
    stiffness_matrix = np.delete(stiffness_matrix, np.s_[0:2], 1)

    return stiffness_matrix


def make_mass_matrix(n, length, breadth, height, density):

    area = breadth*height
    mass_matrix = np.zeros((2*n + 2, 2*n + 2), dtype='float')

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

    mass_matrix = np.delete(mass_matrix, np.s_[0:2], 0)
    mass_matrix = np.delete(mass_matrix, np.s_[0:2], 1)

    return mass_matrix


def find_frequencies(mass_matrix, stiffness_matrix):

    evals, evecs = la.eigh(stiffness_matrix, mass_matrix)

    return np.sqrt(evals)/(2*np.pi)


def ode(t, y, n, a, b, forces, freq, amplitude, type, spring_vec, gap=0):

    forces[-2] = amplitude*type(2*np.pi*freq*t)

    yvec = np.array([[y[i] for i in range(4*n)]]).T

    yvec1 = np.matmul(a, yvec) + np.matmul(b, forces)
    yvec2 = np.matmul(a, yvec) + np.matmul(b, forces)


    if y[2*n - 2] < gap:
        return yvec2
    else:
        return yvec1
    return yvec2


def state_space_formulation(n, mass_matrix, stiffness_matrix, start_t, end_t, delta_t, spring_k):

    t_vector = np.arange(start_t, end_t, delta_t)
    t_interval = np.array([start_t, end_t])
    init_conditions = np.zeros(4*n)
    spring_vec = np.zeros(4*n)

    spring_vec[n - 1] = spring_k

    damping = np.zeros((2*n, 2*n))

    a = np.vstack([
        np.hstack([
            np.zeros((2*n, 2*n)), 
            np.eye(2*n, 2*n)
        ]), 
        np.hstack([
            -np.matmul(np.linalg.inv(mass_matrix), stiffness_matrix), 
            -np.matmul(np.linalg.inv(mass_matrix), damping)
        ])
    ])

    b = np.vstack([
        np.zeros((2*n, 2*n)), 
        np.linalg.inv(mass_matrix)
    ])

    vertical_forces = np.zeros((2*n, 1))

    return [init_conditions, damping, a, b, vertical_forces, t_vector, t_interval, spring_vec]



def main():
    number = 4                     # number of elements taken
    length_of_beam = 5              # in meters
    breadth_of_beam = 0.05          # in meters
    height_of_beam = 0.05           # in meters
    modulus_of_elasticity = 7E10    
    density_of_beam = 2700  

    global_mass_matrix = make_mass_matrix(
        number, 
        length_of_beam, 
        breadth_of_beam, 
        height_of_beam, 
        density_of_beam
    )

    global_stiffness_matrix = make_stiffness_matrix(
        number, 
        length_of_beam, 
        breadth_of_beam, 
        height_of_beam, 
        modulus_of_elasticity
    )

    frequencies = find_frequencies(
        global_mass_matrix, 
        global_stiffness_matrix
    )

    print("The first 5 frequencies are : ", frequencies[:5])

    start_time = 0
    end_time = 4
    time_step = 1/2000

    force_frequency = frequencies[0]
    force_amplitude = 25
    force_type = np.sin

    spring_stiffness = 0
    beam_spring_gap = 0.1

    s_s_f_r = state_space_formulation(
        number, 
        global_mass_matrix, 
        global_stiffness_matrix, 
        start_time, 
        end_time, 
        time_step, 
        spring_stiffness
    )

    initial_conditions = s_s_f_r[0] 
    damping_matrix = s_s_f_r[1]
    a_matrix = s_s_f_r[2]
    b_matrix = s_s_f_r[3]
    force_vector = s_s_f_r[4]
    time_vector = s_s_f_r[5]
    time_interval = s_s_f_r[6]
    spring_vector = s_s_f_r[7]

    sol = solve_ivp(
        ode,
        time_interval, 
        initial_conditions, 
        t_eval=time_vector, 
        vectorized=True, 
        args=(
            number, 
            a_matrix, 
            b_matrix, 
            force_vector, 
            force_frequency, 
            force_amplitude, 
            force_type, 
            beam_spring_gap,
            spring_vector
        )
    )

    time = sol.t
    total_beam_profile = sol.y

    plot_graph(
        time, 
        total_beam_profile[2*number - 2], 
        'time in s', 
        'displacement in m', 
        'displacement at end',
        gap=-beam_spring_gap
    )



if __name__ == '__main__':
    main()

