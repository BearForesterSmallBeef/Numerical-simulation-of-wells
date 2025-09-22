import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal
# беспонятия как это работает, но это решает ангем для трехдиагональной матрицы


class QuantumWellSolver:
    def __init__(self, z_min, z_max, N, a=2.0, b=2.0):
        self.z_min = z_min
        self.z_max = z_max
        self.N = N
        self.a = a
        self.b = b

        self.z, self.dz = np.linspace(z_min, z_max, N, retstep=True)

        self.potential = self._define_potential()

        self.eigenvalues = None
        self.eigenvectors = None

    def _define_potential(self):
        pass

    def solve(self, num_levels=10):

        # элементы главной диагонали матрицы, здесь же прописываем начальные условия
        diagonal = ([self.a] + [2] * (self.N - 4) + [self.b]) / self.dz ** 2 + self.potential[1:-1]

        # элементы побочных диагоналей
        off_diagonal = -1 / self.dz ** 2 * np.ones(self.N - 2)

        self.eigenvalues, self.eigenvectors = eigh_tridiagonal(
            diagonal,
            off_diagonal[1:],
            select='i',
            select_range=(0, num_levels - 1)
        )
        # подгон размерности - тут тоже играет роль a и b, но тут надо как-то переписывать
        full_eigenvectors = np.zeros((self.N, num_levels))
        full_eigenvectors[1:-1, :] = self.eigenvectors
        self.eigenvectors = full_eigenvectors

        for i in range(num_levels):
            norm = np.sqrt(np.sum(self.eigenvectors[:, i] ** 2) * self.dz)
            self.eigenvectors[:, i] /= norm

        print([f"E_{i} = {E}" for i, E in enumerate(self.eigenvalues)], sep=";\t")

    def plot_potential(self, path="results/pivo_pot.png"):
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

        axs.plot(self.z, self.potential, color="gray", label=r"$U(z)$")
        plt.title("Профиль потенциальной ямы")
        plt.xlabel(r"$z$")
        plt.ylabel(r"$U(z)$")
        plt.legend()
        plt.savefig(path)

    def plot_wavefunctions(self, num_levels_to_plot=10, path="results/pivo_wave.png"):
        if self.eigenvalues is None:
            print("Ничо не посчитано")
            return
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

        plt.plot(self.z, self.potential, color='black', label='$U(z)$', lw=2)

        num_levels = min(num_levels_to_plot, len(self.eigenvalues))

        for i in range(num_levels):
            plt.plot(self.z, self.eigenvalues[i] + self.eigenvectors[:, i], label=f"E_{i}={round(self.eigenvalues[i], 2)}")
            plt.axhline(self.eigenvalues[i], color="gray", linestyle='--')

        plt.title("Собственные функции и уровни энергии")
        plt.xlabel(r"$z$")
        plt.ylabel(r"$\psi(z)$")
        plt.legend()
        plt.savefig(path)


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


parabola = ParabolicWell(-5, 5, 2000)
parabola.solve()
parabola.plot_wavefunctions(path="results/parabola/1_wave.png")

half_parabola = HalfParabolicWell(5, 2000)
half_parabola.solve(num_levels=5)
half_parabola.plot_wavefunctions(path="results/half_parabola/1_wave.png")
