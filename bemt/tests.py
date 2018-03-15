import numpy as np
import propeller
import quadrotor
import aero_coeffs
import unit_conversion
import bemt


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
radius = unit_conversion.in2m(9.0)/2
print "radius = " + str(radius)
root_cutout = 0.1 * radius
dy = float(radius-root_cutout)/n_elements
dr = float(1)/n_elements
y = root_cutout + dy*np.arange(1, n_elements+1)
r = y/radius
pitch = 0.0
airfoils = (('SDA1075_494p', 0.0, 1.0),)
allowable_Re = []
#allowable_Re = [1000000., 500000., 250000., 100000., 90000., 80000., 70000., 60000., 50000., 40000., 30000., 20000., 10000.]
vehicle_weight = 11.5
max_chord = 0.6
alt = 0
tip_loss = True
mach_corr = False

# Forward flight parameters
v_inf = 20.     # m/s
alpha0 = 6. * np.pi / 180  # Starting guess for trimmed alpha in radians
n_azi_elements = 10

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
if Cl_tables and allowable_Re:
    Cl_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cl_fun(Re, Cl_tables[airfoils[0][0]], Clmax[airfoils[0][0]][Re]) for Re in allowable_Re]))
    Cd_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cd_fun(Re, Cd_tables[airfoils[0][0]]) for Re in allowable_Re]))

###########################################
# Set design variable bounds
###########################################
omega = 5943.0 * 2*np.pi/60
chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787,
                       0.1786, 0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154]) * radius
chord = np.array([chord[i] for i in [0, 2, 4, 6, 8, 10, 12, 14, 15, 17]])
twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009,
                       19.814, 18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
twist = np.array([twist[i] for i in [0, 2, 4, 6, 8, 10, 12, 14, 15, 17]])

prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy, airfoils=airfoils, Cl_tables=Cl_tables,
                           Cd_tables=Cd_tables)

quad = quadrotor.Quadrotor(prop, vehicle_weight)

print "quad created"

dT_h, P_h = bemt.bemt_axial(prop, pitch, omega, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                            tip_loss=True, mach_corr=mach_corr, alt=alt)
print "axial bemt finished"
print "hover thrust = " + str(sum(dT_h))
alpha = 0.
T_trim, H_trim, P_trim, converged = bemt.bemt_forward_flight(quad, pitch, omega, alpha, v_inf, n_azi_elements, alt=alt, tip_loss=tip_loss,
                                           mach_corr=mach_corr, allowable_Re=allowable_Re, Cl_funs=Cl_funs,
                                           Cd_funs=Cd_funs)

print "hover thrust = " + str(sum(dT_h))
print "hover power = " + str(P_h)
print "ff thrust = " + str(T_trim)
print "ff power = " + str(P_trim)