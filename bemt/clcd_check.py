from lookup_table import create_table
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt


airfoil = 'SDA1075'
this_Re = 90000
alpha, Re, CL, CD = create_table(airfoil)
print "Re from table = " + str(Re)
print "Alpha from table = " + str(alpha)

Cl_table = ((alpha, Re), CL)
Cd_table = ((alpha, Re), CD)


def getCl(aoa, R):
    return griddata(Cl_table[0], Cl_table[1], (aoa, R))


def getCd(aoa, R):
    return griddata(Cd_table[0], Cd_table[1], (aoa, R))

# max_alpha = 0
# min_alpha = 100
# for pair in zip(alpha, Re):
#     if abs(pair[1] - this_Re) < 0.05:
#         if pair[0] > max_alpha:
#             max_alpha = pair[0]
#         elif pair[0] < min_alpha:
#             min_alpha = pair[0]


this_alphas = [Cl_table[0][0][i] for i in xrange(len(Cl_table[0][0])) if abs(Cl_table[0][1][i] - this_Re) < 0.05]
max_alpha = max(this_alphas)
min_alpha = min(this_alphas)
print "max_alpha = " + str(max_alpha)
print "min_alpha = " + str(min_alpha)

print getCl(15.9, 90000)
print getCd(15.9, 90000)

alpha_vec = np.linspace(0., 12., 15)
print alpha_vec

Cl = [getCl(a, this_Re) for a in alpha_vec]

# plt.figure(1)
# plt.plot(alpha_vec, Cl)
# plt.xlabel('Angle of attack')
# plt.ylabel('Cl')
# plt.show()
