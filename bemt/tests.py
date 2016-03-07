import numpy as np
import matplotlib.pyplot as plt
import propeller
import bemt


def ideal_twist(tip_theta, r):
    return tip_theta/r


def linear_twist(theta_0, total_twist, r):
    return theta_0 + r * total_twist


def ideal_taper(tip_chord, r):
    return tip_chord / r


def linear_taper(root_chord, tip_chord, r):
    return root_chord - (root_chord-tip_chord)*r

n_elements = 40
n_azi_elements = 60
n_blades = 2
blade_radius = 6
root_cutout = 0.1
dy = float(blade_radius-root_cutout)/n_elements
dr = float(1)/n_elements
y = root_cutout + dy*np.arange(1, n_elements+1)
r = y/blade_radius

twist = 8 * 2 * np.pi / 360 * np.ones(n_elements)
chord = 0.4 * np.ones(n_elements)
omega_rpm = 400
omega_radps = omega_rpm * 2 * np.pi / 60
alpha_deg = 8
alpha_rad = alpha_deg * 2 * np.pi / 360
Clalpha = 2 * np.pi

# Ghazirah test cases
pitch = 0
CT0 = 0.002
speeds = np.linspace(1, 100, 50)
prop = propeller.Propeller(twist, chord, blade_radius, n_blades, r, y, dr, dy, Clalpha)
prop_CTs = []
prop_CPs = []
num_iter = []
CT_azifun_mat = []
power_req = []
mu = speeds * np.cos(alpha_rad) / (omega_radps*blade_radius)
for speed in speeds:
    CT_old = CT0
    converged = False
    n_iter = 0
    while not converged:
        n_iter += 1
        CT, CP, dCt, Ptot, local_inflow, rel_inflow_angle, dCl, CT_azifun = bemt.bemt(prop, pitch, omega_radps, alpha_rad,
                                                                           v_inf=speed, CT0=CT_old, mach_corr=False)
        converged = abs((CT - CT_old) / CT) < 0.0005
        if converged:
            prop_CTs.append(CT)
            prop_CPs.append(CP)
            num_iter.append(n_iter)
            CT_azifun_mat.append(CT_azifun)
            power_req.append(Ptot)
        CT_old = CT

for i in xrange(len(speeds)):
    print "(v_inf, mu, CT, iter) = (%f, %f, %f, %d)" % (speeds[i], mu[i], prop_CTs[i], num_iter[i])

plt.figure(1)
plt.plot(speeds*1.944, [745.699872*pwr for pwr in power_req])
plt.xlabel("Airspeed, kts")
plt.ylabel("Power required, watts")



plt.show()

# # Leishman 3.2
# solidity = 0.1
# prop = propeller.Propeller(twist_zero, chord, blade_radius, n_blades, r, y, dr, dy, Clalpha, solidity=solidity)
# CT_target = 0.008
# start_pitch = 6*CT_target/(solidity*Clalpha)+3*np.sqrt(2)/4*np.sqrt(CT_target)
# dCT_list = []
# CT_list = []
# inflow_list = []
# for tip_loss in [True, False]:
#     next_pitch = start_pitch
#     converged = False
#     while not converged:
#         CT, CP, dCT, Ptot, local_inflow, rel_inflow_angle, dCl = bemt.bemt(prop, next_pitch, omega_rad, alpha_rad,
#                                                                       tip_loss=tip_loss)
#         next_pitch += (6*(CT_target - CT)/(solidity*Clalpha) + 3*np.sqrt(2)/4*(np.sqrt(CT_target)-np.sqrt(CT)))
#         converged = abs((CT - CT_target)/CT) < 0.0005
#         if converged:
#             CT_list.append(CT)
#             dCT_list.append(dCT)
#             inflow_list.append(local_inflow)
#
#
# plt.figure(1)
# plt.plot(r, dCT_list[0]/dr, 'b', r, dCT_list[1]/dr, 'r')
# plt.xlabel("Radial position, r")
# plt.ylabel("Thrust gradient, dCT/dr")
# plt.legend(["Tip losses", "No tip losses"], loc="upper left")
# plt.grid(True)
#
# plt.figure(2)
# plt.plot(r, inflow_list[0], 'b', r, inflow_list[1], 'r')
# plt.xlabel("Radial distribution, r")
# plt.ylabel("Induced inflow ratio")
# plt.legend(["Tip losses", "No tip losses"], loc="upper left")
# plt.grid(True)

#plt.show()

# # Leishman twist verification
# solidity = 0.1
# CT_target = 0.008
# start_pitch = 6*CT_target/(solidity*Clalpha)+3*np.sqrt(2)/4*np.sqrt(CT_target)
# root_twist = 12*2*np.pi/360
# tip_theta = 4*CT_target/solidity/Clalpha+np.sqrt(CT_target/2)
# print "Tip theta: " + str(tip_theta)
# inflow_list = []
# dCT_list = []
# dCl_list = []
#
# twist_neg10 = linear_twist(root_twist, -10*2*np.pi/360, r)
# twist_neg20 = linear_twist(root_twist, -20*2*np.pi/360, r)
# twist_zero = np.zeros(n_elements)
# twist_ideal = ideal_twist(4*CT_target/solidity/Clalpha+np.sqrt(CT_target/2), r)
# twists = [twist_ideal]
#     #, twist_zero, twist_neg10, twist_neg20]
#
#
# prop = propeller.Propeller(twist_zero, chord, blade_radius, n_blades, r, y, dr, dy, Clalpha, solidity=solidity)
# for twist in twists:
#     prop.twist = twist
#     next_pitch = start_pitch
#     converged = False
#     while not converged:
#         CT, CP, dCT, Ptot, local_inflow, rel_inflow_angle, dCl = bemt.bemt(prop, next_pitch, omega_rad, alpha_rad,
#                                                                            tip_loss=False)
#         next_pitch += (6*(CT_target - CT)/(solidity*Clalpha) + 3*np.sqrt(2)/4*(np.sqrt(CT_target)-np.sqrt(CT)))
#         converged = abs((CT - CT_target)/CT) < 0.0005
#         if converged:
#             inflow_list.append(local_inflow)
#             dCT_list.append(dCT)
#             dCl_list.append(dCl)
#
# plt.figure(1)
# for x in inflow_list:
#     plt.plot(r, x)
# plt.xlabel("Radial distribution, r")
# plt.ylabel("Inflow ratio")
# plt.legend(["Ideal", "Zero", "-10", "-20"], loc='upper left')
# plt.ylim([0, 0.1])
# plt.grid(True)
#
# plt.figure(2)
# for x in dCT_list:
#     plt.plot(r, x/dr)
# plt.xlabel("Radial distribution, r")
# plt.ylabel("Thrust gradient, dCT/dr")
# plt.legend(["Ideal", "Zero", "-10", "-20"], loc='upper left')
# plt.ylim([0, 0.03])
# plt.grid(True)
#
# plt.figure(3)
# for x in dCl_list:
#     plt.plot(r, x)
# plt.xlabel("Radial distribution, r")
# plt.ylabel("Sectional lift coefficient")
# plt.legend(["Ideal", "Zero", "-10", "-20"], loc='upper left')
# plt.ylim([0, 2])
# plt.grid(True)
#
# plt.show()








