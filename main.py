import numpy as np

class Protein():
    def __init__(self, params, presequence="XXXX"):
        self.params = params
        self.presequence = presequence
        self.sequence = self.get_sequence()


    def force(self,  type=None):
        k_b, b, b_zero, k_a, a, a_zero, k_phi, n_phi, phi, delta, r, qi, qj, a, c  = self.params
        if type == "bond":
            return 2 * k_b * (b - b_zero)

        if type == "angle":
            return 2 * k_phi * (a - a_zero)

        if type == "torsion":
            return -k_phi * n_phi * np.sin(n_phi * phi + delta)

        if type == "elec_lenard":
            return -qi * qj / (r ^ 2) - a / (r ^ 13) + c(r ^ 7)

        if type == None:
            pass


    def amino(self, n):
        amino = self.sequence[n]
        if amino == "X":
            return np.zeros((3, 3))
        else:
            pass


    def get_sequence(self):
        # use X for tri-chain of carbon c-c-c
        word = "XXXXX"  # "DAEFRHDSGYEVHHQKLVFFAEDVGSNK"
        return [char for char in self.presequence]  # get 1AMC Peptide smallest part of AMYLOID Beta


    def get_position(self, n):
        if n-1 < 0:
            pass
        else:
            pass


    def set_position(self, n):
        k_b, b, b_zero, k_a, a, a_zero, k_phi, n_phi, phi, delta, r, qi, qj, a, c = self.params
        pos = self.amino(n)
        if n > 0:
            pos_bef = self.amino(n-1)
        if n%2 == 0:
            zickzack = -1
        else:
            zickzack = 1
        for i in range(len(pos[:, 0])):
            b = np.random.randint(130, 160)
            a = np.radians(np.random.randint(90, 180))
            if i == 0:
                if n > 0:
                    pos[0, i] = pos_bef[0, 2] + (b * np.sin(-a * zickzack))
                    pos[1, i] = pos_bef[1, 2] + (b * np.cos(-a * zickzack))
            if i == 1:
                pos[0, i] = pos[0, i-1] + (b * np.sin(a * zickzack))
                pos[1, i] = pos[1, i-1] + (b * np.cos(a * zickzack))
            if i == 2:
                pos[0, i] = pos[0, i-1] + (b * np.sin(-a * zickzack))
                pos[1, i] = pos[1, i-1] + (b * np.cos(-a * zickzack))
        return pos

        # else:
        #     pos_bond = self.amino(n - 1)
        #     # bond_dist = sqrt((pos[0, i] - pos[0, i]) ^ 2 + (ya - yb) ^ 2 + (za - zb) ^ 2)
        #     bond = self.force(type = "bond")
        #     angle = self.force(type = "angle")
        #     pass


    def initialize_positions(self, n):
        k_b, b, b_zero, k_a, a, a_zero, k_phi, n_phi, phi, delta, r, qi, qj, a, c = self.params
        pos = self.amino(n)
        if n > 0:
            pos_bef = self.amino(n-1)
        if n % 2 == 0:
            zickzack = -1
        else:
            zickzack = 1
        for i in range(len(pos[:, 0])):
            b = np.random.randint(130, 160)
            a = np.radians(np.random.randint(90, 180))
            if i == 0:
                if n > 0:
                    pos[0, i] = pos_bef[0, 2] + (b * np.cos(-a * zickzack))
                    pos[1, i] = pos_bef[1, 2] + (b * np.sin(-a * zickzack))
            if i == 1:
                pos[0, i] = pos[0, i-1] + (b * np.cos(a * zickzack))
                pos[1, i] = pos[1, i-1] + (b * np.sin(a * zickzack))
            if i == 2:
                pos[0, i] = pos[0, i-1] + (b * np.cos(-a * zickzack))
                pos[1, i] = pos[1, i-1] + (b * np.sin(-a * zickzack))
        return pos


    def structure(self):
        structure = []
        for n in range(len(self.sequence)):
            print(self.initialize_positions(n))
            structure.append(self.initialize_positions(n))
        return structure


k_b = 1
b_zero = 147  # in pm
k_a = 1
a_zero = np.radians(120)  # in deg
b = np.random.randint(130, 160)
a = np.radians(np.random.randint(90, 180))
k_phi, n_phi, phi, delta, r, qi, qj, a, c = 1, 1, 1, 1, 1, 1, 1, 1, 1


params = [k_b, b, b_zero, k_a, a, a_zero, k_phi, n_phi, phi, delta, r, qi, qj, a, c]

carbonchain = Protein(params, presequence="XXXXX")

print(carbonchain.structure())
