import numpy as np

class Protein():
    def __init__(self, params, presequence="XXXX"):
        self.params = params
        self.presequence = presequence
        self.sequence = self.get_sequence()


    def force(self, pos, n, atom, type=None):  # pos must be array [3, 3, n]

        k_b, b, b_zero, k_a, a, a_zero, k_phi, n_phi, phi, delta, r, qi, qj, a, c = self.params
        if type == "bond":
            if atom != 0:
                b_x = pos[0, atom, n] - pos[0, atom - 1, n]
                b_y = pos[1, atom, n] - pos[1, atom - 1, n]
            elif n != 0:  # previous amino
                b_x = pos[0, atom, n] - pos[0, 2, n-1]
                b_y = pos[1, atom, n] - pos[1, 2, n-1]
            else:  # redundant, n=0, atom=0 needs to be forbidden
                print("Oups this schould not be used!")
                b_x = pos[0, atom, n] - pos[0, 0, 0]  # should return zero
                b_y = pos[1, atom, n] - pos[1, 0, 0]  # should return zero
            # calculation of force angle
            angle_b_x_y = np.arctan(b_x / b_y)
            b_zero_x = b_zero * np.cos(angle_b_x_y)
            b_zero_y = b_zero * np.sin(angle_b_x_y)
            f_bound_x, f_bound_y = 2 * k_b * (b_x - b_zero_x), 2 * k_b * (b_y - b_zero_y)
            return f_bound_x, f_bound_y


        if type == "angle":
           if n != 0:
                if atom == 0:

                    d = np.sqrt((pos[0, atom, n] - pos[0, 2, n - 1])^2 + (pos[1, atom, n] - pos[1, 2, n - 1])^2)       # distance between first and second atom

                    e = np.sqrt((pos[0, atom+1, n] - pos[0, atom, n])^2 + (pos[1, atom +1, n] - pos[1, atom, n ])^2)   # distance between second and third atom

                    f = np.sqrt((pos[0, atom+1, n] - pos[0, 2, n - 1])^2 + (pos[1, atom + 1, n] - pos[1, 2, n - 1])^2) # distance between first and third atom


                elif atom == 1:

                    d = np.sqrt((pos[0, atom, n] - pos[0, atom - 1, n - 1]) ^ 2 + (pos[1, atom, n] - pos[1, atom - 1, n]) ^ 2)  # distance between first and second atom

                    e = np.sqrt((pos[0, atom + 1, n] - pos[0, atom, n]) ^ 2 + (pos[1, atom + 1, n] - pos[1, atom, n]) ^ 2)  # distance between second and third atom

                    f = np.sqrt((pos[0, atom + 1, n] - pos[0, atom - 1, n]) ^ 2 + (pos[1, atom + 1, n] - pos[1, atom - 1, n]) ^ 2)  # distance between first and third atom

                elif atom == 2:
                    d = np.sqrt((pos[0, atom, n] - pos[0, atom - 1, n - 1]) ^ 2 + (pos[1, atom, n] - pos[1, atom - 1, n]) ^ 2)  # distance between first and second atom

                    e = np.sqrt((pos[0, 0, n + 1] - pos[0, atom, n]) ^ 2 + (pos[1, 0, n + 1] - pos[1, atom, n]) ^ 2)  # distance between second and third atom

                    f = np.sqrt((pos[0, 0, n + 1] - pos[0, atom - 1, n]) ^ 2 + (pos[1, 0, n + 1] - pos[1, atom - 1, n]) ^ 2)  # distance between first and third atom

        a = np.arccos(d ^ 2 + e ^ 2 - f ^ 2 / (2 * d * e))
        return 2 * k_phi * (a - a_zero), 2 * k_phi * (a - a_zero)  # torque

        if type == "torsion":  # not included yet, 2d impossible
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
        return pos(n)



    def update_position(self, pos, t):
        k_b, b, b_zero, k_a, a, a_zero, k_phi, n_phi, phi, delta, r, qi, qj, a, c = self.params
        mass = 1  # relative mass, if masses of atoms are different
        for i in range(len(self.sequence)):
            for j in range(len(pos[:, 0, 0])):
                if i == 0 and j == 0:
                    pos[:, 0, 0]
                else:
                    acc_x, acc_y = self.force(pos, i, j, type="bond") / mass
                    distance_x = acc_x * t
                    distance_y = acc_y * t
                    pos[0, i, j] += distance_x
                    pos[1, i, j] += distance_y
        for i in range(len(self.sequence)):
            for j in range(len(pos[:, 0, 0])):
                acc_rad = self.force(pos, i, j, type="angle")
                # winkelberechnung 
        return pos

    def initialize_positions(self):
        structure = []
        for n in range(len(self.sequence)):
            k_b, b, b_zero, k_a, a, a_zero, k_phi, n_phi, phi, delta, r, qi, qj, a, c = self.params
            pos = self.amino(n)
            if n > 0:
                pos_bef = self.amino(n - 1)
            if n % 2 == 0:
                zickzack = -1  # changes direction of rotation
            else:
                zickzack = 1
            for i in range(len(pos[:, 0])):
                b = np.random.randint(130, 160)  # create random distance between atoms
                a = np.radians(np.random.randint(90, 180))  # create random angle between atom (in radiant)
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
            np.array(structure.append(pos))
        return structure





k_b = 1
b_zero = 147  # in pm
k_a = 1
a_zero = np.radians(120)  # in deg
b = np.random.randint(130, 160)
a = np.radians(np.random.randint(90, 180))
k_phi, n_phi, phi, delta, r, qi, qj, a, c = 1, 1, 1, 1, 1, 1, 1, 1, 1


params = [k_b, b, b_zero, k_a, a, a_zero, k_phi, n_phi, phi, delta, r, qi, qj, a, c]

carbonchain = Protein(params, presequence="XXXXXXX")

pos = carbonchain.initialize_positions()
print(pos)

#np.dot(b_x,b_y)/ (np.linalg.norm(b_x) * np.linalg.norm(b_y))