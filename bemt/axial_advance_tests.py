import numpy as np
import matplotlib.pyplot as plt
import propeller
import bemt


# Geometry for DA4002 9x6.75
n_blades = 2
diameter = 9 * 0.0254
blade_radius = diameter / 2
root_cutout = 0.1 * blade_radius
chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
                  0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154]) * blade_radius
twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
                  18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
n_elements = chord.size
pitch = 0

dy = float(blade_radius-root_cutout)/n_elements
dr = float(1)/n_elements
y = root_cutout + dy*np.arange(1, n_elements+1)
r = y/blade_radius
Clalpha = 2 * np.pi

thin_electric_10 = propeller.Propeller(twist, chord, blade_radius, n_blades, r, y, dr, dy, Clalpha,
                                       airfoils=(('NACA4412', 0, 1),))
DA4002 = propeller.Propeller(twist, chord, blade_radius, n_blades, r, y, dr, dy, Clalpha, airfoils=(('SDA1075', 0, 1),))

omega_rpm = 5000
omega_radps = omega_rpm * 2 * np.pi / 60
n = omega_rpm / 60
adv_ratios = np.arange(0.15, 0.85, 0.01)
print len(adv_ratios)

CT_array = np.empty([adv_ratios.size])
CP_array = np.empty([adv_ratios.size])
CQ_array = np.empty([adv_ratios.size])
Re_array = np.empty([adv_ratios.size])
T_array = np.empty([adv_ratios.size])
high_rpm_Cl = 0
high_rpm_dT = 0
prop_CT = 0
prop_CP = 0
prop_Re = 0
prop_aoa = 0
for i in xrange(adv_ratios.size):
    v_climb = adv_ratios[i] * n * diameter
    CT, CP, CQ, Cl, dT, pCT, pCP, Re, aoa = bemt.bemt_axial(DA4002, pitch, omega_radps, v_climb=v_climb)
    T_array[i] = sum(dT)
    CT_array[i] = pCT
    CP_array[i] = pCP
    CQ_array[i] = CQ
    Re_array[i] = Re[3*len(Re)/4]

# Experimental data of CT and CP vs J for the DA4002 at nominally 5000 RPM
exp_CT = [0.130239, 0.127207, 0.124940, 0.122901, 0.119574, 0.115206, 0.111528, 0.106872, 0.102531, 0.098328, 0.093013,
          0.088433, 0.082563, 0.076371, 0.071175, 0.063821, 0.057435, 0.052166, 0.045150, 0.038378, 0.032477, 0.025451,
          0.019905, 0.010939, 0.003606]
exp_CP = [0.076611, 0.075522, 0.074664, 0.074731, 0.074178, 0.073309, 0.072770, 0.071610, 0.070558, 0.069457, 0.067796,
          0.066534, 0.064268, 0.061482, 0.059543, 0.055666, 0.052418, 0.049923, 0.045885, 0.041691, 0.038005, 0.032932,
          0.028840, 0.021637, 0.015514]
exp_J = [0.123969, 0.155062, 0.186357, 0.218201, 0.250515, 0.281467, 0.312083, 0.343845, 0.374369, 0.404401, 0.436725,
         0.468579, 0.500784, 0.531197, 0.565824, 0.593765, 0.625372, 0.659234, 0.689432, 0.721506, 0.751023, 0.783654,
         0.810712, 0.847761, 0.876248]

print Re_array

plt.figure(1)
plt.plot(adv_ratios, CT_array, 'ro', exp_J, exp_CT, 'go')
plt.xlabel("J")
plt.ylabel("Thrust Coefficient")
plt.grid(b='on')
plt.xlim(0, 1.0)

plt.show()