import numpy as np
from unit_conversion import in2m
import propeller
import quadrotor
import aero_coeffs
import unit_conversion
import bemt
import matplotlib.pyplot as plt
import trim


def calc_twist_dist(t0, dt_vec):
    t = np.array([0.0]*(len(dt_vec)+1))
    t[0] = t0
    for i, dt in enumerate(dt_vec):
        t[i+1] = t[i] + dt
    return t


def calc_chord_dist(c0, dc_vec):
    c = np.array([0.0]*(len(dc_vec)+1))
    c[0] = c0
    for i, dc in enumerate(dc_vec):
        c[i+1] = c[i] + dc
    return c

n_blades = 2
n_elements = 22
radius = unit_conversion.in2m(15.)/2
root_cutout = 0.08 * radius
dy = float(radius-root_cutout)/n_elements
dr = float(1)/n_elements
y = root_cutout + dy*np.arange(1, n_elements+1)
r = y/radius
pitch = 0.0
airfoils = (('SDA1075_494p', 0.0, 1.0),)
allowable_Re = [1000000., 500000., 250000., 100000., 90000., 80000., 70000., 60000., 50000., 40000., 30000., 20000., 10000.]
vehicle_weight = 26.7
alt = 0
tip_loss = True
mach_corr = False

# Forward flight parameters
v_inf = 4.    # m/s
alpha0 = 0.0454  # Starting guess for trimmed alpha in radians
n_azi_elements = 5

# Mission times
time_in_hover = 5. * 60     # Time in seconds
time_in_ff = 10. * 60
mission_time = [time_in_hover, time_in_ff]


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
lift_curve_info_dict = {}
if any(Cl_tables) and allowable_Re:
    Cl_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cl_fun(Re, Cl_tables[airfoils[0][0]], Clmax[airfoils[0][0]][Re]) for Re in allowable_Re]))
    Cd_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cd_fun(Re, Cd_tables[airfoils[0][0]]) for Re in allowable_Re]))
    lift_curve_info_dict = aero_coeffs.create_liftCurveInfoDict(allowable_Re, Cl_tables[airfoils[0][0]])


radius = in2m(15./2)
twist = np.array([6.97, 15.97, 21.77, 21.72, 19.91, 18.14, 16.55, 15.28, 14.01, 13., 12.18, 11.39, 10.76, 10.24, 9.85,
                  9.4, 9.07, 8.7, 8.46, 8.29, 8.19, 8.17]) * 2 * np.pi / 360
chord = in2m(np.array([0.77, 0.97, 1.18, 1.34, 1.44, 1.49, 1.50, 1.49, 1.46, 1.42, 1.38, 1.33, 1.27, 1.21, 1.14, 1.07, 1.,
                       0.93, 0.86, 0.77, 0.67, 0.44]))

prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy, airfoils=airfoils, Cl_tables=Cl_tables,
                           Cd_tables=Cd_tables)

quad = quadrotor.Quadrotor(prop, vehicle_weight)

ff_kwargs = {'propeller': prop, 'pitch': pitch, 'n_azi_elements': n_azi_elements, 'allowable_Re': allowable_Re,
             'Cl_funs': Cl_funs, 'Cd_funs': Cd_funs, 'tip_loss': tip_loss, 'mach_corr': mach_corr, 'alt': alt,
             'lift_curve_info_dict': lift_curve_info_dict}

omega = 2800 * 2*np.pi/60
Th, Ph = bemt.bemt_axial(prop, pitch, omega, allowable_Re=allowable_Re, Cd_funs=Cd_funs, Cl_funs=Cl_funs)
print "(T, P) for %d RPM = (%f, %f)" % (omega*60/2/np.pi, sum(Th)*0.2248*4, Ph)

omega = 4200 * 2*np.pi/60
Th, Ph = bemt.bemt_axial(prop, pitch, omega, allowable_Re=allowable_Re, Cd_funs=Cd_funs, Cl_funs=Cl_funs)
print "(T, P) for %d RPM = (%f, %f)" % (omega*60/2/np.pi, sum(Th)*0.2248*4, Ph)

## Hover performance of SUI Endurance for a range of rotor speeds
# hthrust = []
# hpower = []
# o_start = 2000.
# o_end = 4600.
# o_range = np.linspace(o_start, o_end, 10)
# for o in o_range*2*np.pi/60:
#     dT_h, P_h = bemt.bemt_axial(prop, pitch, o, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
#                                 tip_loss=True, mach_corr=mach_corr, alt=alt)
#     hthrust.append(dT_h.sum())
#     hpower.append(P_h)
#     print "(Omega, Thrust, Power) = (%f, %f, %f)" % (o*60/2/np.pi, dT_h.sum(), P_h)
#
# plt.figure(1)
# plt.plot(o_range, hthrust)
# plt.xlabel("RPM")
# plt.ylabel("Thrust (N)")
#
# plt.figure(2)
# plt.plot(o_range, hpower)
# plt.xlabel("RPM")
# plt.ylabel("Thrust (N)")
#
# plt.show()

lift_exp = np.array([[3.2, 3.85, 4.1, 4.15, 4.2], [5.25, 6.05, 6.25, 6.35, 6.5], [7.75, 8.65, 8.9, 9.2, 9.25]])
lift_exp = np.array([lift_exp[0], lift_exp[2]])
drag_exp = np.array([[-0.8, -0.25, 0., 0.25, 0.3], [-1.25, -0.5, -0.1, 0.2, 0.3], [-2.2, -0.9, -0.25, 0.2, 0.45]])
drag_exp = np.array([drag_exp[0], drag_exp[2]])
alphas_exp = np.array([-20., -10., -5., -2.5, 0.])[::-1]

q = 22.98252   # dynamic pressure in N/m**2
dens = 1.225
v_inf = np.sqrt(2*q/dens)
frame_drag = 0.27 * 0.092903 * q
omegas = np.array([2800., 4200.])
alphas = np.linspace(20., 0., 10)[::-1]
lift = np.empty([len(omegas), len(alphas)], dtype=float)
drag = np.empty([len(omegas), len(alphas)], dtype=float)
for j, o in enumerate(omegas*2*np.pi/60):
    for k, a in enumerate(alphas*2*np.pi/360):
        T, H, P = bemt.bemt_forward_flight(quad, pitch, o, a, v_inf, n_azi_elements, alt=alt,
                                           tip_loss=tip_loss, mach_corr=mach_corr, allowable_Re=allowable_Re,
                                           Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                           lift_curve_info_dict=lift_curve_info_dict)
        lift[j, k] = T*np.cos(a) + H*np.sin(a)
        drag[j, k] = H*np.cos(a) - T*np.sin(a) + frame_drag
        #print "drag = " + str(H*np.cos(a) - T*np.sin(a) + frame_drag)
        #print "frame drag = " + str(frame_drag)
lift *= 0.224809
drag *= 0.224809

# color_str = ['ko-', 'k*-', 'kv-', 'ks-', 'kD-']
# color_str_disc = ['ko-', 'k*-', 'kv-', 'ks-', 'kD-']
# l = []
# plt.figure(1)
# for i in xrange(len(omegas)):
#     lift[i, :][5:] *= 0.98
#     plt.plot(alphas*-1, lift[i, :]*4, color_str[i], markerfacecolor='white', markersize=8)
#     l.append('%d RPM' % omegas[i])
# for i in xrange(len(omegas)):
#     plt.plot(alphas_exp, lift_exp[i, :], color_str_disc[i], markerfacecolor='black', markersize=8)
#     l.append('%d RPM exp' % omegas[i])
# plt.xlabel('Alpha, deg')
# plt.ylabel('Lift, lb')
# plt.legend(l, loc='lower right')
# plt.ylim([0, 11])
# plt.grid()

plt.figure(1)
lift[0, :][5:] *= 1
lift[1, :][5:] *= 1
plt.plot(alphas*1, lift[0, :]*4, 'ko-', markerfacecolor='white')
plt.plot(alphas_exp*-1, lift_exp[0, :][::-1], 'ks-', markerfacecolor='white')
plt.plot(alphas*1, lift[1, :]*4, 'ko-')
plt.plot(alphas_exp*-1, lift_exp[1, :][::-1], 'ks-')
plt.xlabel(r'$\alpha,\,\mathrm{deg}$', fontsize=18)
plt.ylabel(r'$Lift,\,lb$', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tick_params(axis='both', which='minor', labelsize=14)
plt.legend(['2800 RPM BEMT', '2800 RPM exp', '4200 RPM BEMT', '4200 RPM exp'], loc='lower left')
plt.ylim([0, 11])
plt.grid()

plt.figure(2)
plt.plot(alphas*1, drag[0, :]*4, 'ko-', markerfacecolor='white')
plt.plot(alphas_exp*-1, drag_exp[0, :][::-1], 'ks-', markerfacecolor='white')
plt.plot(alphas*1, drag[1, :]*4, 'ko-')
plt.plot(alphas_exp*-1, drag_exp[1, :][::-1], 'ks-')
plt.xlabel(r'$\alpha,\,\mathrm{deg}$', fontsize=18)
plt.ylabel(r'$Drag,\,lb$', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tick_params(axis='both', which='minor', labelsize=14)
plt.legend(['2800 RPM BEMT', '2800 RPM exp', '4200 RPM BEMT', '4200 RPM exp'], loc='lower left')
plt.ylim([-5, 2])
plt.grid()

plt.show()

# l = []
# plt.figure(2)
# for i in xrange(len(omegas)):
#     plt.plot(alphas*-1, drag[i, :]*4+frame_drag, color_str[i], markerfacecolor='white')
#     l.append('%d RPM' % omegas[i])
# for i in xrange(len(omegas)):
#     plt.plot(alphas_exp, drag_exp[i, :], color_str_disc[i], markerfacecolor='white')
#     l.append('%d RPM exp' % omegas[i])
# plt.xlabel('Alpha, deg')
# plt.ylabel('Drag, lb')
# #l[1], l[2] = l[2], l[1]
# plt.legend(l, loc='lower right')
# plt.ylim([-4, 1])
# plt.grid()



# # Trim study
# trim0 = [alpha0, omega]
# alpha_trim, omega_trim, converged = trim.trim(quad, v_inf, trim0, ff_kwargs)
#
# T_ff, H_ff, P_ff = bemt.bemt_forward_flight(quad, pitch, omega_trim, alpha_trim, v_inf, n_azi_elements, alt=alt,
#                                             tip_loss=tip_loss, mach_corr=mach_corr, allowable_Re=allowable_Re,
#                                             Cl_funs=Cl_funs, Cd_funs=Cd_funs,
#                                             lift_curve_info_dict=lift_curve_info_dict)
#
# print "(a_trim_new, o_trim_new) = (%f, %f)" % (alpha_trim, omega_trim*60/2/np.pi)
# print "Hover (thrust, power) = (%f, %f)" % (sum(dT_h), P_h)
# print "FFnew (T, H, power)   = (%f, %f, %f)" % (T_ff, H_ff, P_ff)

omega = 4200 * 2*np.pi/60
a = 0.0000001
T_ff, H_ff, P_ff = bemt.bemt_forward_flight(quad, pitch, omega, a, v_inf, n_azi_elements, alt=alt,
                                            tip_loss=tip_loss, mach_corr=mach_corr, allowable_Re=allowable_Re,
                                            Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                            lift_curve_info_dict=lift_curve_info_dict)
dTh, Ph = bemt.bemt_axial(prop, pitch, omega, allowable_Re=allowable_Re, Cd_funs=Cd_funs, Cl_funs=Cl_funs)