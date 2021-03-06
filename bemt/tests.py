import numpy as np
import matplotlib.pyplot as plt
import propeller
import quadrotor
import aero_coeffs
import unit_conversion
import bemt
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
n_elements = 10
radius = unit_conversion.in2m(9.6)/2
root_cutout = 0.1 * radius
dy = float(radius-root_cutout)/n_elements
dr = float(1)/n_elements
y = root_cutout + dy*np.arange(1, n_elements+1)
r = y/radius
pitch = 0.0
airfoils = (('SDA1075_494p', 0.0, 1.0),)
#allowable_Re = []
allowable_Re = [1000000., 500000., 250000., 100000., 90000., 80000., 70000., 60000., 50000., 40000., 30000., 20000., 10000.]
vehicle_weight = 12.455
max_chord = 0.6
alt = 0
tip_loss = True
mach_corr = False

# Forward flight parameters
v_inf = 4.    # m/s
alpha0 = 0.0454  # Starting guess for trimmed alpha in radians
n_azi_elements = 5

# Mission times
time_in_hover = 300.     # Time in seconds
time_in_ff = 500.
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

###########################################
# Set design variable bounds
###########################################
# omega = 4250. * 2*np.pi/60
# chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787,
#                        0.1786, 0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154]) * radius
# chord = np.array([chord[i] for i in [0, 2, 4, 6, 8, 10, 12, 14, 15, 17]])
# twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009,
#                        19.814, 18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
# twist = np.array([twist[i] for i in [0, 2, 4, 6, 8, 10, 12, 14, 15, 17]])

chord = np.array([0.11966637, 0.21945255, 0.3191826, 0.41661201, 0.46042987, 0.46404804, 0.36930386, 0.30469107,
                  0.20583294, 0.13521735])
chord *= radius
twist = np.array([0.33294302, 0.40815866, 0.31931257, 0.34385899, 0.32057012, 0.26702793, 0.2698419, 0.23210482,
                  0.19992483, 0.04420424])
omega = 3811.29467393 * 2*np.pi/60

prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy, airfoils=airfoils, Cl_tables=Cl_tables,
                           Cd_tables=Cd_tables)

quad = quadrotor.Quadrotor(prop, vehicle_weight)

ff_kwargs = {'propeller': prop, 'pitch': pitch, 'n_azi_elements': n_azi_elements, 'allowable_Re': allowable_Re,
             'Cl_funs': Cl_funs, 'Cd_funs': Cd_funs, 'tip_loss': tip_loss, 'mach_corr': mach_corr, 'alt': alt,
             'lift_curve_info_dict': lift_curve_info_dict}


# dT_h, P_h = bemt.bemt_axial(prop, pitch, omega, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
#                             tip_loss=True, mach_corr=mach_corr, alt=alt)

trim0 = [alpha0, omega]

for i in xrange(1000):
    alpha_trim, omega_trim, converged = trim.trim(quad, v_inf, trim0, ff_kwargs)
    T_ff, H_ff, P_ff = bemt.bemt_forward_flight(quad, pitch, omega_trim, alpha_trim, v_inf, n_azi_elements, alt=alt,
                                                tip_loss=tip_loss, mach_corr=mach_corr, allowable_Re=allowable_Re,
                                                Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                                lift_curve_info_dict=lift_curve_info_dict)
    print "FFnew (T, H, power)   = (%f, %f, %f)" % (T_ff, H_ff, P_ff)
    print i

# print "(a_trim_new, o_trim_new) = (%f, %f)" % (alpha_trim, omega_trim*60/2/np.pi)
# print "Hover (thrust, power) = (%f, %f)" % (sum(dT_h), P_h)
# print "FFnew (T, H, power)   = (%f, %f, %f)" % (T_ff, H_ff, P_ff)
#
#


# ######################################### Trim analysis ########################################################
# alpha = 0.04
# trim0 = [alpha, omega]
# failed_converge = []
# for v in np.arange(0, 15., 0.1):
#     try:
#         alpha_trim, omega_trim, converged = trim.trim(quad, v, trim0, ff_kwargs)
#         T_ff, H_ff, P_ff = bemt.bemt_forward_flight(quad, pitch, omega_trim, alpha_trim, v, n_azi_elements, alt=alt,
#                                                     tip_loss=tip_loss, mach_corr=mach_corr, allowable_Re=allowable_Re,
#                                                     Cl_funs=Cl_funs, Cd_funs=Cd_funs,
#                                                     lift_curve_info_dict=lift_curve_info_dict)
#         if not converged:
#             print "trim did not converge"
#             failed_converge.append(v)
#         else:
#             print "(alpha_trim, omega_trim, V, P) = (%f, %f, %f, %f)" % (alpha_trim, omega_trim, v, P_ff)
#         trim0 = [alpha_trim, omega_trim]
#     except FloatingPointError:
#         print "Floating point error in trim"
# print "\n"
# print "Speeds that failed to converge: \n"
# print repr(failed_converge)

# alpha = 0.04
# trim0 = [alpha, omega]
# for v in [0.7, 0.9]:
#     try:
#         alpha_trim, omega_trim, converged = trim.trim(quad, v, trim0, ff_kwargs)
#         T_ff, H_ff, P_ff = bemt.bemt_forward_flight(quad, pitch, omega_trim, alpha_trim, v, n_azi_elements, alt=alt,
#                                                     tip_loss=tip_loss, mach_corr=mach_corr, allowable_Re=allowable_Re,
#                                                     Cl_funs=Cl_funs, Cd_funs=Cd_funs,
#                                                     lift_curve_info_dict=lift_curve_info_dict)
#         if not converged:
#             print "trim did not converge"
#         else:
#             print "(alpha_trim, omega_trim, V, P) = (%f, %f, %f, %f)" % (alpha_trim, omega_trim, v, P_ff)
#         trim0 = [alpha_trim, omega_trim]
#     except FloatingPointError:
#         print "Floating point error in trim"