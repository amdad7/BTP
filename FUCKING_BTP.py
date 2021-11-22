from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


def hit_spring(t, y):
    return y[0] + 1


hit_spring.terminal = True


def plot_all_graphs(time_series, solution_of_ode):
    legend_names = [
        'Displacement of m1 vs time',
        'Displacement of m2 vs time',
        'Velocity of m1 vs time',
        'Velocity of m2 vs time'
    ]

    minimum_y_value = min(solution_of_ode)
    maximum_y_value = max(solution_of_ode)
    plt.plot(T, dis_m_1)
    plt.plot(T, np.full((len(dis_m_1)), -1, dtype='int32'))
    plt.vlines(x=TT, ymin=-9, ymax=9)
    plt.title("Displacement of m1 vs time")
    plt.xlabel("Time")
    plt.grid()
    plt.show()

    plt.plot(T, dis_m_2)
    plt.plot(T, np.full((len(dis_m_1)), -1, dtype='int32'))
    plt.vlines(x=TT, ymin=-9, ymax=9)
    plt.title("Displacement of m2 vs time")
    plt.xlabel("Time")
    plt.grid()
    plt.show()

    plt.plot(T, vel_m_1)
    plt.plot(T, np.full((len(dis_m_1)), -1, dtype='int32'))
    plt.vlines(x=TT, ymin=-9, ymax=9)
    plt.title("Velocity of m1 vs time")
    plt.xlabel("Time")
    plt.grid()
    plt.show()


def main():
    pass


if __name__ == '__main__':
    main()
