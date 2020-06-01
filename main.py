import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class Protein():
    def __init__(self, params, presequence="XXXX"):
        self.params = params
        self.presequence = presequence
        self.sequence = self.get_sequence()

    def force(self, pos, n, atom, type=None):  # pos must be array [3, 3, n]

        k_b, b, b_zero, k_a, a, a_zero, k_phi, n_phi, phi, delta, r, qi, qj, a, c = self.params
        x_ax = 0
        y_ax = 1

        if type == "bond":
            if atom != 0:
                b_x = pos[n][x_ax, atom] - pos[n][x_ax, atom - 1]
                b_y = pos[n][y_ax, atom] - pos[n][y_ax, atom - 1]
            elif n != 0:  # previous amino
                b_x = pos[n][x_ax, atom] - pos[n - 1][x_ax, 2]
                b_y = pos[n][y_ax, atom] - pos[n - 1][y_ax, 2]
            else:  # redundant, n=0, atom=0 needs to be forbidden
                print("Oups this schould not be used!")
                b_x = pos[n][x_ax, atom] - pos[0][x_ax, 0]  # should return zero
                b_y = pos[n][y_ax, atom] - pos[0][y_ax, 0]  # should return zero
            # calculation of force angle
            angle_b_x_y = np.arctan(b_x / b_y)
            b_zero_x = b_zero * np.cos(angle_b_x_y)
            b_zero_y = b_zero * np.sin(angle_b_x_y)
            f_bound_x, f_bound_y = 2 * k_b * (b_x - b_zero_x), 2 * k_b * (b_y - b_zero_y)
            return f_bound_x, f_bound_y

        if type == "angle":
            if n == 0 and atom == 0:
                return 0
            else:
                if atom == 0:

                    d = np.sqrt((pos[n][x_ax, atom] - pos[n - 1][x_ax, 2]) ** 2 + (
                                pos[n][y_ax, atom] - pos[n - 1][y_ax, 2]) ** 2)  # distance between first and second atom

                    e = np.sqrt((pos[n][x_ax, atom + 1] - pos[n][x_ax, atom]) ** 2 + (pos[n][y_ax, atom + 1] - pos[n][
                        y_ax, atom]) ** 2)  # distance between second and third atom

                    f = np.sqrt((pos[n][x_ax, atom + 1] - pos[n - 1][x_ax, 2]) ** 2 + (
                                pos[n][y_ax, atom + 1] - pos[n - 1][
                            y_ax, 2]) ** 2)  # distance between first and third atom


                elif atom == 1:

                    d = np.sqrt((pos[n][x_ax, atom] - pos[n][x_ax, atom - 1]) ** 2 + (pos[n][y_ax, atom] - pos[n][
                        y_ax, atom - 1]) ** 2)  # distance between first and second atom

                    e = np.sqrt((pos[n][x_ax, atom + 1] - pos[n][x_ax, atom]) ** 2 + (pos[n][y_ax, atom + 1] - pos[n][y_ax, atom]) ** 2)  # distance between second and third atom

                    f = np.sqrt((pos[n][x_ax, atom + 1] - pos[n][x_ax, atom - 1]) ** 2 + (pos[n][y_ax, atom + 1] - pos[n][y_ax, atom - 1]) ** 2)  # distance between first and third atom

                elif atom == 2:
                    d = np.sqrt((pos[n][x_ax, atom] - pos[n][x_ax, atom - 1]) ** 2 + (pos[n][y_ax, atom] - pos[n][
                        y_ax, atom - 1]) ** 2)  # distance between first and second atom

                    e = np.sqrt((pos[n + 1][x_ax, 0] - pos[n][x_ax, atom]) ** 2 + (
                                pos[n + 1][y_ax, 0] - pos[n][y_ax, atom]) ** 2)  # distance between second and third atom

                    f = np.sqrt((pos[n + 1][x_ax, 0] - pos[n + 1][x_ax, atom - 1]) ** 2 + (pos[n + 1][y_ax, 0] - pos[n][
                        y_ax, atom - 1]) ** 2)  # distance between first and third atom

                a = np.arccos((d ** 2 + e ** 2 - f ** 2) / (2 * d * e))
                return 2 * k_phi * (a - a_zero)  # torque
        if type == "torsion":  # not included yet, 2d impossible
            return -k_phi * n_phi * np.sin(n_phi * phi + delta)

        if type == "elec_lenard":
            return -qi * qj / (r ** 2) - a / (r ** 13) + c(r ** 7)

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

        x_ax = 0
        y_ax = 1

        for i in range(len(self.sequence)):  # i = n
            for j in range(len(pos[0])-1):  # j = atom
                if i == 0 and j == 0:
                    pass
                elif i == len(pos)-1 and j == 2:
                    pass
                else:
                    acc_x, acc_y = self.force(pos, i, j, type="bond")  # calc boundforce
                    distance_x = acc_x * t
                    distance_y = acc_y * t
                    pos[i][x_ax, j] += distance_x
                    pos[i][y_ax, j] += distance_y

        for i in range(len(self.sequence)):  # i = n
            for j in range(len(pos[0])-1):  # j = atom
                acc_rad = self.force(pos, i, j, type="angle") # calc angle force
                angle_to_y_ax = (np.pi / 2) - (acc_rad/mass) * t

                if i == 0 and j == 0:
                    pass
                elif i == len(pos)-1 and j == 2:
                    pass
                else:
                    if j == 0:
                        e_dist = np.sqrt(
                            (pos[i][x_ax, j + 1] - pos[i][x_ax, j]) ** 2 + (pos[i][y_ax, j + 1] - pos[i][y_ax, j]) ** 2)
                        pos[i][x_ax, j + 1] = pos[i - 1][x_ax, 2] + e_dist * np.cos(angle_to_y_ax)
                        pos[i][y_ax, j + 1] = pos[i - 1][y_ax, 2] + e_dist * np.sin(angle_to_y_ax)
                    elif j == 1:
                        e_dist = np.sqrt((pos[i][x_ax, j + 1] - pos[i][x_ax, j]) ** 2 + (pos[i][y_ax, j + 1] - pos[i][
                            y_ax, j]) ** 2)  # distance between second and third atom
                        pos[i][x_ax, j + 1] = pos[i][x_ax, j - 1] + e_dist * np.cos(angle_to_y_ax)
                        pos[i][y_ax, j + 1] = pos[i][y_ax, j - 1] + e_dist * np.sin(angle_to_y_ax)
                    elif j == 2:
                        e_dist = np.sqrt((pos[i + 1][x_ax, 0] - pos[i][x_ax, j]) ** 2 + (pos[i + 1][y_ax, 0] - pos[i][
                            y_ax, j]) ** 2)  # distance between second and third atom
                        pos[i+1][x_ax, 0] = pos[i + 1][x_ax, 0] + e_dist * np.cos(angle_to_y_ax)
                        pos[i+1][y_ax, 0] = pos[i + 1][y_ax, 0] + e_dist * np.sin(angle_to_y_ax)
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
                    pos[0, i] = pos[0, i - 1] + (b * np.cos(a * zickzack))
                    pos[1, i] = pos[1, i - 1] + (b * np.sin(a * zickzack))
                if i == 2:
                    pos[0, i] = pos[0, i - 1] + (b * np.cos(-a * zickzack))
                    pos[1, i] = pos[1, i - 1] + (b * np.sin(-a * zickzack))
            structure.append(pos)  # np.array()
        return structure


k_b = 1
b_zero = 147  # in pm
k_a = 1
a_zero = np.radians(120)  # in deg
b = np.random.randint(130, 160)
a = np.radians(np.random.randint(70, 110))
k_phi, n_phi, phi, delta, r, qi, qj, a, c = 1, 1, 1, 1, 1, 1, 1, 1, 1
time_step = 0.01

params = [k_b, b, b_zero, k_a, a, a_zero, k_phi, n_phi, phi, delta, r, qi, qj, a, c]

carbonchain = Protein(params, presequence="XXXXXXX")

pos = carbonchain.initialize_positions()

print(
    f" pos Array: {len(pos)} \n Pos[0] Array: {pos[0]} \n pos[0][0] Array: {pos[0][0]} \n pos[0][0,1] Array: {pos[0][0, 1]}")




# visualization using Matplotlib:



fig = plt.figure()
ax = plt.axes(xlim=(0, 4), ylim=(-2, 2))
line, = ax.plot([], [], lw=3)

def init():
    line.set_data([], [])
    return line,

def animate(i,pos):
    pos = carbonchain.update_position(pos, time_step)

    x_data = []
    y_data = []

    for h in range(len(pos)): # count n
        for k in range(len(pos[0])):
            x = pos[h][0,k]
            y = pos[h][1,k]
            x_data.append(x)
            y_data.append(y)
    line.set_data(x_data, y_data)
    return line,
    #return x_data,y_data, pos

"""
for test in range(100):
   x,y,pos = animate(test,pos)
   fig = plt.figure() 
   ax = fig.add_subplot(111)
   #ax.set_xlim(0,1e9)
   ax.plot(x,y,"o-")
   plt.savefig("C:\\Users\Jacob\Desktop\Testfolding/FIG" + str(test) +".png")
"""

anim = animation.FuncAnimation(fig, animate, init_func=init,
							frames=500, interval=20, blit=True)


anim.save('coil.gif',writer='imagemagick')