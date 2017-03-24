import numpy as np
import matplotlib.pyplot as plt
import propeller
import bemt
import unit_conversion
import lookup_table
import aero_coeffs


####################
# Experimental Data
####################

# DA4002 Experimental Static (J = 0 i.e., hover)
rpm_static = np.array([1546.667, 1700.000, 2053.333, 2253.333, 2503.333, 2746.667, 2946.667, 3280.000, 3526.667,
                       3743.333, 3940.000, 4213.333, 4443.333, 4676.667, 4930.000, 5186.667, 5430.000, 5710.000,
                       5943.333])
CT_exp_static = np.array([0.125931, 0.127314, 0.129033, 0.127408, 0.131518, 0.130653, 0.132682, 0.130863, 0.132955,
                          0.134311, 0.134936, 0.136319, 0.135364, 0.138118, 0.138094, 0.138747, 0.139023, 0.139617,
                          0.140450])
CP_exp_static = np.array([0.088432, 0.088785, 0.087836, 0.085943, 0.087917, 0.086259, 0.086557, 0.083443, 0.083443,
                          0.083174, 0.082488, 0.079168, 0.078000, 0.079092, 0.078738, 0.079026, 0.079163, 0.079595,
                          0.080081])

# DA4002 Experimental Axial RPM = 2013
J_2013rpm = np.array([0.315488, 0.392566, 0.465086, 0.543643, 0.642669, 0.688683, 0.779962, 0.894262])
CT_exp_2013rpm = np.array([0.103248, 0.095053, 0.083917, 0.068049, 0.046900, 0.034373, 0.014050, -0.006909])
CP_exp_2013rpm = np.array([0.073956, 0.071617, 0.067243, 0.059680, 0.047706, 0.039326, 0.025041, 0.008631])

# DA4002 Experimental Axial RPM = 3030
J_3030rpm = np.array([0.202688, 0.254532, 0.318318, 0.368036, 0.417834, 0.465380, 0.518250, 0.565931, 0.635838,
                      0.682688, 0.734338, 0.787576, 0.836754, 0.887498])
CT_exp_3030rpm = np.array([0.116273, 0.111396, 0.104832, 0.097600, 0.091369, 0.083231, 0.073491, 0.065815, 0.053303,
                           0.042818, 0.030541, 0.018225, 0.005563, -0.004863])
CP_exp_3030rpm = np.array([0.077403, 0.074837, 0.072432, 0.069647, 0.067753, 0.064561, 0.060397, 0.057236, 0.051087,
                           0.044609, 0.036364, 0.027660, 0.017716, 0.008957])

# DA4002 Experimental Axial RPM = 4054
J_4054rpm = np.array([0.151501, 0.191256, 0.232175, 0.269752, 0.307040, 0.346660, 0.387705, 0.426859, 0.465805, 0.507860,
                      0.547924, 0.589943, 0.625818, 0.665307, 0.702598, 0.749396, 0.785937])
CT_exp_4054rpm = np.array([0.123138, 0.119479, 0.116798, 0.112796, 0.107361, 0.103260, 0.096600, 0.091014, 0.085481,
                           0.077791, 0.070923, 0.062089, 0.056236, 0.047958, 0.039769, 0.030508, 0.022348])
CP_exp_4054rpm = np.array([0.074794, 0.074956, 0.074264, 0.072859, 0.071121, 0.070301, 0.068196, 0.066546, 0.064837,
                           0.061803, 0.058943, 0.054657, 0.051932, 0.047361, 0.042256, 0.036304, 0.030356])

# DA4002 Experimental Axial RPM = 4049
J_4049rpm = np.array([0.660673, 0.705412, 0.746229, 0.788221, 0.825365, 0.865364])
CT_exp_4049rpm = np.array([0.047624, 0.040267, 0.031661, 0.022344, 0.012597, 0.003434])
CP_exp_4049rpm = np.array([0.046930, 0.042808, 0.037144, 0.030502, 0.022874, 0.015445])

# DA4002 Experimental Axial RPM = 5028
J_5028rpm = np.array([0.123969, 0.155062, 0.186357, 0.218201, 0.250515, 0.281467, 0.312083, 0.343845, 0.374369, 0.404401,
                      0.436725, 0.468579, 0.500784, 0.531197, 0.565824, 0.593765, 0.625372])
CT_exp_5028rpm = np.array([0.130239, 0.127207, 0.124940, 0.122901, 0.119574, 0.115206, 0.111528, 0.106872, 0.102531,
                           0.098328, 0.093013, 0.088433, 0.082563, 0.076371, 0.071175, 0.063821, 0.057435])
CP_exp_5028rpm = np.array([0.076611, 0.075522, 0.074664, 0.074731, 0.074178, 0.073309, 0.072770, 0.071610, 0.070558,
                           0.069457, 0.067796, 0.066534, 0.064268, 0.061482, 0.059543, 0.055666, 0.052418])

# DA4002 Experimental Axial RPM = 5064
J_5064rpm = np.array([0.530584, 0.564462, 0.596597, 0.626289, 0.659234, 0.689432, 0.721506, 0.751023, 0.783654, 0.810712,
                      0.847761, 0.876248, 0.914534])
CT_exp_5064rpm = np.array([0.077156, 0.070439, 0.064012, 0.058029, 0.052166, 0.045150, 0.038378, 0.032477, 0.025451,
                           0.019905, 0.010939, 0.003606, -0.005644])
CP_exp_5064rpm = np.array([0.062018, 0.059045, 0.055967, 0.052999, 0.049923, 0.045885, 0.041691, 0.038005, 0.032932,
                           0.028840, 0.021637, 0.015514, 0.007598])

# Combined experimental results for the two cases around 5000 RPM
J_5000rpm = np.array([0.123969, 0.155062, 0.186357, 0.218201, 0.250515, 0.281467, 0.312083, 0.343845, 0.374369,
                      0.404401, 0.436725, 0.468579, 0.500784, 0.531197, 0.565824, 0.593765, 0.625372, 0.659234,
                      0.689432, 0.721506, 0.751023, 0.783654, 0.810712, 0.847761, 0.876248, 0.914534])
rpm_5000rpm = np.array([5028.]*len(J_5028rpm)+[5064.]*(len(J_5000rpm)-len(J_5028rpm)))
CT_exp_5000rpm = np.array([0.130239, 0.127207, 0.124940, 0.122901, 0.119574, 0.115206, 0.111528, 0.106872, 0.102531,
                           0.098328, 0.093013, 0.088433, 0.082563, 0.076371, 0.071175, 0.063821, 0.057435, 0.052166,
                           0.045150, 0.038378, 0.032477, 0.025451, 0.019905, 0.010939, 0.003606, -0.005644])
CP_exp_5000rpm = np.array([0.076611, 0.075522, 0.074664, 0.074731, 0.074178, 0.073309, 0.072770, 0.071610, 0.070558,
                           0.069457, 0.067796, 0.066534, 0.064268, 0.061482, 0.059543, 0.055666, 0.052418, 0.049923,
                           0.045885, 0.041691, 0.038005, 0.032932, 0.028840, 0.021637, 0.015514, 0.007598])

#######################
# End Experimental Data
#######################
kine_visc = 1.460 * 10**-5
n_blades = 2
n_elements = 18
radius = unit_conversion.in2m(9.0)/2
root_cutout = 0.1 * radius
dy = float(radius-root_cutout)/n_elements
dr = float(1)/n_elements
y = root_cutout + dy*np.arange(1, n_elements+1)
r = y/radius
pitch = 0.0
airfoils = (('SDA1075_494p', 0.0, 1.0),)
#allowable_Re = [1000000., 500000., 250000., 100000., 90000., 80000., 70000., 60000., 50000., 40000., 30000., 20000., 10000.]
allowable_Re = []
tip_loss = True
mach_corr = False

Cl_tables = {}
Cd_tables = {}
Clmax = {}
# Get lookup tables
if any(airfoil[0] != 'simple' for airfoil in airfoils):
    for airfoil in airfoils:
        Cl_table, Cd_table, Clmax = aero_coeffs.create_Cl_Cd_table(airfoil[0])

        Cl_tables[airfoil[0]] = Cl_table
        Cd_tables[airfoil[0]] = Cd_table
        Clmax[airfoil[0]] = Clmax

# Create list of Cl functions. One for each Reynolds number. Cl_tables (and Cd_tables) will be empty for the
# 'simple' case, therefore this will be skipped for the simple case. For the full table lookup case this will be
# skipped because allowable_Re will be empty.
Cl_funs = {}
Cd_funs = {}
if Cl_tables and allowable_Re:
    Cl_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cl_fun(Re, Cl_tables[airfoils[0][0]], Clmax[airfoils[0][0]][Re]) for Re in allowable_Re]))
    Cd_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cd_fun(Re, Cd_tables[airfoils[0][0]]) for Re in allowable_Re]))


def get_performance(o, c, t, v_climb):
        chord_meters = c * radius
        prop = propeller.Propeller(t, chord_meters, radius, n_blades, r, y, dr, dy, airfoils=airfoils,
                                   Cl_tables=Cl_tables, Cd_tables=Cd_tables)

        perf = bemt.bemt_axial(prop, pitch, o, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                               tip_loss=tip_loss, mach_corr=mach_corr, output='long', v_climb=v_climb)
        return perf[-2], perf[-1]

# DA4002 Propeller with 18 elements
chord_18e = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
                      0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
twist_18e = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
                      18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360

# Construct the CT vs RPM and CP vs RPM for Static Hover
r75_ind = min(range(len(r)), key=lambda i: abs(r[i]-0.75))
r75 = r[r75_ind]
chord75 = chord_18e[r75_ind] * radius
Re75 = [o*2*np.pi*r75*radius*chord75/kine_visc for o in rpm_static*2*np.pi/60]
CT_static = []
CP_static = []
for omega in rpm_static*2*np.pi/60:
    this_CT, this_CP = get_performance(omega, chord_18e, twist_18e, 0)
    CT_static.append(this_CT)
    CP_static.append(this_CP)
CT_static = np.array(CT_static)
CP_static = np.array(CP_static)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(rpm_static, CT_exp_static, 'ko-', markerfacecolor='white')
ax1.plot(rpm_static, CT_static, 'ks-', markerfacecolor='white')
ax1.plot(rpm_static, CP_exp_static, 'ko-')
ax1.plot(rpm_static, CP_static, 'ks-')
ax1.set_xlim([0, 6500])
ax1.set_ylim([0.0, 0.20])
ax1.set_xlabel(r'$\Omega\,\mathrm{in}\,\mathrm{RPM}$', fontsize=18)
ax1.set_ylabel(r'$\mathrm{C}_\mathrm{T},\, \mathrm{ C}_\mathrm{P}$', fontsize=18)
ax1.grid()
ax1.legend([r'$\mathrm{C}_\mathrm{T}\,\mathrm{UIUC}$', r'$\mathrm{C}_\mathrm{T}\,\mathrm{BEMT}\,\mathrm{XFOIL}$',
            r'$\mathrm{C}_\mathrm{P}\,\mathrm{UIUC}$', r'$\mathrm{C}_\mathrm{P}\,\mathrm{BEMT}\,\mathrm{XFOIL}$'], loc=2, fontsize=14)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.tick_params(axis='both', which='minor', labelsize=14)
ax2 = ax1.twiny()
Re6500 = 6500*2*np.pi/60*r75*radius*chord75/kine_visc
ax2.plot(range(int(round(Re6500)/1000)), np.ones(int(round(Re6500)/1000))*0.2)
ax2.set_ylim([0.0, 0.2])
ax2.set_xlabel(r'$\mathrm{Re}_{3/4}\,\mathrm{(10}^\mathrm{3}\mathrm{)}$', fontsize=18)
ax2.tick_params(axis='both', which='major', labelsize=14)
ax2.tick_params(axis='both', which='minor', labelsize=14)
plt.show()

# # Construct the CT vs J and CP vs J for 5000 RPM
# chord75 = chord_18e[int(round(len(r)*0.75))]*radius
# r75 = r[int(round(len(r)*0.75))]
#
# omegas_5000 = rpm_5000rpm * 2*np.pi/60
# n_5000 = rpm_5000rpm / 60
# v_climbs_5000 = np.array([j*n*2*radius for j, n in zip(J_5000rpm, n_5000)])
# Re75_5000 = 5000.*2*np.pi/60**r75*radius*chord75/kine_visc
#
# CT_5000rpm = []
# CP_5000rpm = []
# for omega, v_climb in zip(omegas_5000, v_climbs_5000):
#     this_CT, this_CP = get_performance(omega, chord_18e, twist_18e, v_climb)
#     CT_5000rpm.append(this_CT)
#     CP_5000rpm.append(this_CP)
# CT_5000rpm = np.array(CT_5000rpm)
# CP_5000rpm = np.array(CP_5000rpm)
#
# omega_2000 = 2013. * 2*np.pi/60
# n_2000 = 2013. / 60
# v_climbs_2000 = np.array([j*n_2000*2*radius for j in J_2013rpm])
# Re75_2000 = omega_2000*r75*radius*chord75/kine_visc
# CT_2000rpm = []
# CP_2000rpm = []
# for v_climb in v_climbs_2000:
#     this_CT, this_CP = get_performance(omega_2000, chord_18e, twist_18e, v_climb)
#     CT_2000rpm.append(this_CT)
#     CP_2000rpm.append(this_CP)
# CT_2000rpm = np.array(CT_2000rpm)
# CP_2000rpm = np.array(CP_2000rpm)
#
# Go
# CT_2000rpm = CT_2000rpm[CT_2000rpm > 0]
# CP_2000rpm = CP_2000rpm[CP_2000rpm > 0]
#
# omega_3000 = 3030. * 2*np.pi/60
# n_3000 = 3030. / 60
# v_climbs_3000 = np.array([j*n_3000*2*radius for j in J_3030rpm])
# Re75_3000 = omega_3000*r75*radius*chord75/kine_visc
# print Re75_3000
# CT_3000rpm = []
# CP_3000rpm = []
# for v_climb in v_climbs_3000:
#     this_CT, this_CP = get_performance(omega_3000, chord_18e, twist_18e, v_climb)
#     CT_3000rpm.append(this_CT)
#     CP_3000rpm.append(this_CP)
# CT_3000rpm = np.array(CT_3000rpm)
# CP_3000rpm = np.array(CP_3000rpm)
#
#
# plt.figure(1)
# plt.plot(J_5000rpm, CT_exp_5000rpm, 'ko-', markerfacecolor='white')
# plt.plot(J_5000rpm, CT_5000rpm, 'ko-')
# plt.plot(J_3030rpm, CT_exp_3030rpm, 'kv-', markerfacecolor='white')
# plt.plot(J_3030rpm, CT_3000rpm, 'kv-')
# plt.plot(J_2013rpm, CT_exp_2013rpm, 'ks-', markerfacecolor='white')
# plt.plot(J_2013rpm, CT_2000rpm, 'ks-')
# #plt.ylim([0.0, 0.15])
# plt.xlim([0.0, 1.0])
# plt.xlabel(r'$\mathrm{J}$', fontsize=18)
# plt.ylabel(r'$\mathrm{C}_\mathrm{T}$', fontsize=18)
# plt.grid()
# plt.legend(['UIUC 5000rpm', 'BEMT 5000rpm', 'UIUC 3000rpm', 'BEMT 3000rpm', 'UIUC 2000rpm', 'BEMT 2000rpm'],
#            loc='lower left')
# plt.tick_params(axis='both', which='major', labelsize=14)
# plt.tick_params(axis='both', which='minor', labelsize=14)
#
# plt.figure(2)
# plt.plot(J_5000rpm, CP_exp_5000rpm, 'ko-', markerfacecolor='white')
# plt.plot(J_5000rpm, CP_5000rpm, 'ko-')
# plt.plot(J_3030rpm, CT_exp_3030rpm, 'kv-', markerfacecolor='white')
# plt.plot(J_3030rpm, CT_3000rpm, 'kv-')
# plt.plot(J_2013rpm, CT_exp_2013rpm, 'ks-', markerfacecolor='white')
# plt.plot(J_2013rpm, CT_2000rpm, 'ks-')
# #plt.ylim([0.0, 0.10])
# plt.xlim([0.0, 1.0])
# plt.xlabel(r'$\mathrm{J}$', fontsize=18)
# plt.ylabel(r'$\mathrm{C}_\mathrm{P}$', fontsize=18)
# plt.grid()
# plt.legend(['UIUC 5000rpm', 'BEMT 5000rpm', 'UIUC 3000rpm', 'BEMT 3000rpm', 'UIUC 2000rpm', 'BEMT 2000rpm'])
# plt.tick_params(axis='both', which='major', labelsize=14)
# plt.tick_params(axis='both', which='minor', labelsize=14)
#
# plt.show()


























