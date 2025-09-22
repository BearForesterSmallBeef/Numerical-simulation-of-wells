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


parabola = ParabolicWell(-5, 5, 2000)
parabola.solve()
parabola.plot_wavefunctions(path="results/parabola/1_wave.png")

half_parabola = HalfParabolicWell(5, 2000)
half_parabola.solve(num_levels=5)
half_parabola.plot_wavefunctions(path="results/half_parabola/1_wave.png")
