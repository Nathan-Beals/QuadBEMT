from lookup_table import create_table, interpolate, interpolate_griddata, interp_weights
from scipy.interpolate import griddata, interp2d
import numpy as np
import matplotlib.pyplot as plt


airfoil = 'SDA1075_494p'
this_Re = 90000
alpha, Re, CL, CD = create_table(airfoil)

Cl_table = ((alpha, Re), CL)
Cd_table = ((alpha, Re), CD)

# Cl_table = interp2d(alpha, Re, CL)
# Cd_table = interp2d(alpha, Re, CD)

this_Re = [10940.86555083, 13419.20972376, 21080.74104361, 29969.83072693, 36438.89917782, 41373.73264073,
           45842.63697827, 50839.13217244, 55877.23014428, 60612.78159125, 65482.99369756, 70332.56842598,
           75197.56024198, 80340.67175443, 85377.27859881, 90330.64760795, 89890.63129127, 8576.42658423]

this_alpha = [18.64167177, 24.07602684, 22.66423279, 17.7874695, 15.73267543, 14.17419872, 13.62414577, 13.25507582,
              12.86662159, 12.49683968, 12.0619762, 11.66907308, 11.28712304, 10.99017599, 10.7326172, 10.53891917,
              8.73624857, 1.69642724]

# Cl_vec = [getCl(this_alpha[i], this_Re[i]) for i in xrange(len(this_Re))]
# Cd_vec = [getCd(this_alpha[i], this_Re[i]) for i in xrange(len(this_Re))]
# print "Cl = " + str(Cl_vec)
# print "Cd = " + str(Cd_vec)


# this_alphas = [alpha[i] for i in xrange(len(alpha)) if abs(Re[i] - this_Re) < 0.05]
# max_alpha = max(this_alphas)
# min_alpha = min(this_alphas)
# print "max_alpha = " + str(max_alpha)
# print "min_alpha = " + str(min_alpha)

print "griddata Cl = " + str(interpolate_griddata(Cl_table, 5.0, 80000))
print "griddata Cd = " + str(interpolate_griddata(Cd_table, 5.0, 80000))

vtx_Cl, wts_Cl = interp_weights(Cl_table[0])

alpha_vec = np.linspace(0., 12., 15)
print alpha_vec

# plt.figure(1)
# plt.plot(alpha_vec, Cl)
# plt.xlabel('Angle of attack')
# plt.ylabel('Cl')
# plt.show()
