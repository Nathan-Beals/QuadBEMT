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


# # Geometry for 9x4.7 apc slow flyer
# n_azi_elements = 60
# n_blades = 2
# blade_radius = 8.90 * 0.0254 / 2
# root_cutout = 0.1 * blade_radius
# chord = np.array([0.127, 0.135, 0.158, 0.178, 0.195, 0.209, 0.219, 0.225, 0.227, 0.226, 0.221, 0.212, 0.199, 0.182,
#                   0.161, 0.135, 0.097, 0.058]) * blade_radius
# twist = np.array([27.54, 25.28, 26.08, 25.47, 24.07, 22.18, 20.00, 18.18, 16.38, 14.83, 13.63, 12.56, 11.56, 10.65,
#                   9.68, 8.51, 6.72, 4.89]) * 2 * np.pi / 360
# n_elements = chord.size
# pitch = 0

# Geometry for 10x5 apc thin electric
n_blades = 2
blade_radius = 10 * 0.0254 / 2
root_cutout = 0.1 * blade_radius
chord = np.array([0.130, 0.149, 0.173, 0.189, 0.197, 0.201, 0.200, 0.194, 0.186, 0.174, 0.160, 0.145, 0.128, 0.112,
                  0.096, 0.081, 0.061, 0.041]) * blade_radius
twist = np.array([32.76, 37.19, 33.54, 29.25, 25.64, 22.54, 20.27, 18.46, 17.05, 15.97, 14.87, 14.09, 13.39, 12.84,
                  12.25, 11.37, 10.19, 8.99]) * 2 * np.pi / 360
n_elements = chord.size
pitch = 0

dy = float(blade_radius-root_cutout)/n_elements
dr = float(1)/n_elements
y = root_cutout + dy*np.arange(1, n_elements+1)
r = y/blade_radius
Clalpha = 2 * np.pi

thin_electric_10 = propeller.Propeller(twist, chord, blade_radius, n_blades, r, y, dr, dy, Clalpha,
                                       airfoils=(('NACA4412', 0, 1),))

# omegas = np.array([2763, 3062, 3310, 3622, 3874, 4153, 4422, 4687, 4942, 5226, 5473, 5736, 6026, 6285, 6554, 6768]) \
#          * 2 * np.pi / 60
omegas = np.array([2508, 2795, 3093, 3344, 3630, 3919, 4176, 4446, 4743, 5025, 5314, 5596, 5869, 6146, 6434, 6708]) \
         * 2 * np.pi / 60
CT_array = np.empty([omegas.size])
CP_array = np.empty([omegas.size])
CQ_array = np.empty([omegas.size])
high_rpm_Cl = 0
high_rpm_dT = 0
prop_CT = 0
prop_CP = 0
prop_Re = 0
prop_aoa = 0
for i in xrange(omegas.size):
    CT, CP, CQ, Cl, dT, pCT, pCP, Re, aoa = bemt.bemt_axial(thin_electric_10, pitch, omegas[i])
    CT_array[i] = CT
    CP_array[i] = CP
    CQ_array[i] = CQ
    if i == omegas.size-1:
        high_rpm_Cl = Cl
        high_rpm_dT = dT
        prop_CT = pCT
        prop_CP = pCP
        prop_Re = Re
        prop_aoa = aoa
print "CT(RPM) = " + str(CT_array)
print "CP(RPM) = " + str(CP_array)
print "CQ(RPM) = " + str(CQ_array)
print "propCT = " + str(prop_CT)
print "propCP = " + str(prop_CP)

# plt.figure(1)
# plt.plot(omegas * 60 / 2 / np.pi, CT_array)
# plt.xlabel("\Omega, RPM")
# plt.ylabel("Thrust Coefficient")
# plt.xlim([1500, 7500])
# plt.ylim([0, 0.15])
#
# plt.figure(2)
# plt.plot(omegas * 60 / 2 / np.pi, CP_array)
# plt.xlabel("\Omega, RPM")
# plt.ylabel("Power Coefficient")
# plt.xlim([1500, 7500])
# plt.ylim([0, 0.10])
#
# plt.figure(3)
# plt.plot(r, high_rpm_Cl)
# plt.xlabel("Radial location")
# plt.ylabel("Lift coefficient")
#
# plt.figure(4)
# plt.plot(r, high_rpm_dT)
# plt.xlabel("Radial location")
# plt.ylabel("Thrust distribution")

plt.figure(5)
plt.plot(r, prop_Re)
plt.xlabel("Radial location")
plt.ylabel("Reynolds number")

plt.figure(6)
plt.plot(r, prop_aoa*360/2/np.pi)
plt.xlabel("Radial location")
plt.ylabel("Angle of attack, deg")


plt.show()


# pitch = 0
# CT0 = 0.002
# speeds = np.linspace(1, 100, 50)
# alt = 1585
# prop = propeller.Propeller(twist, chord, blade_radius, n_blades, r, y, dr, dy, Clalpha)
# prop_CTs = []
# prop_Pos = []
# prop_Pis = []
# num_iter = []
# CT_azifun_mat = []
# power_req = []
# mu = speeds * np.cos(alpha_rad) / (omega_radps*blade_radius)
# for speed in speeds:
#     CT_old = CT0
#     converged = False
#     n_iter = 0
#     while not converged:
#         n_iter += 1
#         CT, Pi, Po, dCt, Ptot, local_inflow, rel_inflow_angle, dCl, CT_azifun = bemt.bemt(prop, pitch, omega_radps, alpha_rad,
#                                                                            v_inf=speed, CT0=CT_old, mach_corr=False)
#         converged = abs((CT - CT_old) / CT) < 0.0005
#         if converged:
#             prop_CTs.append(CT)
#             prop_Pis.append(Pi)
#             prop_Pos.append(Po)
#             num_iter.append(n_iter)
#             CT_azifun_mat.append(CT_azifun)
#             power_req.append(Ptot)
#         CT_old = CT
#
# for i in xrange(len(speeds)):
#     print "(v_inf, mu, CT, iter) = (%f, %f, %f, %d)" % (speeds[i], mu[i], prop_CTs[i], num_iter[i])
#
# plt.figure(1)
# plt.plot(speeds*1.944, [pwr/745.699872 for pwr in power_req])
# plt.xlabel("Airspeed, kts")
# plt.ylabel("Power required, hp")
#
# plt.figure(2)
# plt.plot(speeds*1.944, [p/745.699872 for p in prop_Pos])
# plt.plot(speeds*1.944, [p/745.699872 for p in prop_Pis])
# plt.xlabel("Airspeed, kts")
# plt.ylabel("Induced and profile power, hp")
# plt.legend(["Profile", "Induced"], loc="upper left")
#
# plt.show()

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








