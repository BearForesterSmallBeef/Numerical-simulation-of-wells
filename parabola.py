import numpy as np

from BasicSolver import *


class ParabolicWell(QuantumWellSolver):
    def __init__(self, z_min, z_max, N, a=2.0, b=2.0, A=1.0):
        self.A = A
        super().__init__(z_min, z_max, N, a, b)

    def _define_potential(self):
        return self.A * self.z ** 2


class HalfParabolicWell(QuantumWellSolver):
    def __init__(self, z_max, N, a=2.0, b=2.0, z_min=0, A=1.0):
        if z_min != 0:
            z_min = 0
        self.A = A
        super().__init__(z_min, z_max, N, a, b)

    def _define_potential(self):
        potential = self.A * self.z ** 2
        return potential


def check_dependence(energies: np.ndarray, gipotesis=lambda x: 2 * x + 1, path="results/pivo_gipotesis.png"):
    fig, axs = plt.subplots()
    plt.tick_params(axis='both', labelsize=16)
    fig.set_size_inches(12, 9)
    font1 = {'family': 'serif',
             'color': 'black',
             'weight': 'normal',
             'size': 20,
             }
    font2 = {'family': 'serif',
             'color': 'black',
             'weight': 'normal',
             'size': 24,
             }
    c = ["r", "g", "b", "black", "yellow"]
    axs.plot(range(len(energies)), energies, label='Моделирование')
    axs.plot(range(len(energies)), [gipotesis(i) for i in range(len(energies))], label='Модель')
    plt.title("Собственные функции и уровни энергии")
    plt.xlabel(r"n")
    plt.ylabel(r"$\frac{E}{\frac{hw}{2}}$")
    plt.legend()
    plt.savefig(path)


parabola = ParabolicWell(-50, 50, 20000)
parabola.solve(num_levels=100)
parabola.plot_wavefunctions(path="results/parabola/1_wave.png")
check_dependence(parabola.get_energy(), lambda x: 2*x + 1, "results/parabola/1_gip.png")

half_parabola = HalfParabolicWell(50, 20000)
half_parabola.solve(num_levels=50)
half_parabola.plot_wavefunctions(path="results/half_parabola/1_wave.png")
check_dependence(half_parabola.get_energy(), lambda x: 4*x + 3, "results/half_parabola/1_gip.png")
