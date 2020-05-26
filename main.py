import numpy as np

class Protein():
    def __init__(self):
        pass

    def amino(self, amino):

        np.zeros((4, 1))  # position + atom_number

    def get_sequence(self, amino=-1):
        word = "DAEFRHDSGYEVHHQKLVFFAEDVGSNK"
        return [char for char in word]  # get 1AMC Peptide smallest part of AMYLOID Beta

    def get_position(self, sequence, n):
        pass

    def set_position(self, sequence, n):
        pass

    def initialize_positions(self):

        sequence = self.get_sequence()

        for i in range(len(sequence)):
            if i == 0:
                amino = sequence[i]
                self.amino(amino)



class Force():
    def __init__(self):
        pass

    def bond(self, k_b, b, b_zero):
        return 2 * k_b * (b - b_zero)

    def angle(self, k_phi, a, a_zero):
        return 2 * k_phi * (a - a_zero)

    def torsion(self, k_psi, n, psi, delta):
        return -k_psi * n * np.sin(n * psi + delta)

    def elec_lenard(self, r, qi, qj, a, c):
        return -qi * qj / (r ^ 2) - a/(r ^ 13) + c(r ^ 7)


Force_Field = Force()












