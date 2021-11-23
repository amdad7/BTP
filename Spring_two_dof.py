from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


def hit_spring(t, y):
    return y[0] + 1


hit_spring.terminal = True


def plot_all_graphs(time_series, solution_of_ode, legends):
    for i in range(len(solution_of_ode)):
        plt.plot(time_series, solution_of_ode[i])
        plt.plot(time_series, np.full((len(solution_of_ode[i])), -1, dtype='int32'))
        plt.title(legends[i])
        plt.xlabel("Time")
        plt.grid()
        plt.show()


def main():
    pass


if __name__ == '__main__':
    main()
