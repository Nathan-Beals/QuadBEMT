from pyOpt import Optimization
from pyOpt import SLSQP
from pyOpt import NSGA2
from pyOpt import CONMIN
import numpy as np
import propeller
import bemt
import unit_conversion
import matplotlib.pyplot as plt
import aero_coeffs


def counted(f):
    def wrapped(*args, **kwargs):
        print wrapped.calls
        wrapped.calls += 1
        return f(*args, **kwargs)
    wrapped.calls = 0
    return wrapped


@counted
def objfun_optimize_chord(xn, **kwargs):

    radius = kwargs['radius']
    r = kwargs['r']
    y = kwargs['y']
    dr = kwargs['dr']
    dy = kwargs['dy']
    n_blades = kwargs['n_blades']
    airfoils = kwargs['airfoils']
    pitch = kwargs['pitch']
    target_thrust = kwargs['thrust']
    cval_max = kwargs['max_chord']
    tip_loss = kwargs['tip_loss']
    mach_corr = kwargs['mach_corr']
    allowable_Re = kwargs['allowable_Re']
    Cl_tables = kwargs['Cl_tables']
    Cd_tables = kwargs['Cd_tables']
    Cl_funs = kwargs['Cl_funs']
    Cd_funs = kwargs['Cd_funs']
    twist = kwargs['twist']

    # The algorithm works with radius-normalized chord lengths, but the BEMT takes real chord lengths, so multiply by R
    omega = xn[0]
    chord0 = xn[1] * radius
    dchord = np.array([c*radius for c in xn[2:]])

    chord = calc_chord_dist(chord0, dchord)

    f = 1000.
    fail = 0
    g = [1.0] * (2*len(r)+1)

    # Calculate geometric constraint values and return immediately if there are any failures
    g[1:] = get_geocons(chord, cval_max, radius)
    # if any(geocon > 0.0 for geocon in g[1:]):
    #     print "geocons violated"
    #     print g[1:]
    #     return f, g, fail

    prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy,
                               airfoils=airfoils, Cl_tables=Cl_tables, Cd_tables=Cd_tables)

    try:
        dT, P = bemt.bemt_axial(prop, pitch, omega, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                tip_loss=tip_loss, mach_corr=mach_corr)
    except FloatingPointError:
        fail = 1
        return f, g, fail
    except IndexError:
        return 10000, g, fail

    f = P
    print "Power = " + str(f)
    print "Thrust = %s" % str(sum(dT))

    # Evaluate thrust constraint. Target thrust must be less than computed thrust
    g[0] = target_thrust - sum(dT)

    return f, g, fail


@counted
def objfun_optimize_twist(xn, **kwargs):

    radius = kwargs['radius']
    r = kwargs['r']
    y = kwargs['y']
    dr = kwargs['dr']
    dy = kwargs['dy']
    n_blades = kwargs['n_blades']
    airfoils = kwargs['airfoils']
    pitch = kwargs['pitch']
    target_thrust = kwargs['thrust']
    tip_loss = kwargs['tip_loss']
    mach_corr = kwargs['mach_corr']
    allowable_Re = kwargs['allowable_Re']
    Cl_tables = kwargs['Cl_tables']
    Cd_tables = kwargs['Cd_tables']
    Cl_funs = kwargs['Cl_funs']
    Cd_funs = kwargs['Cd_funs']
    chord = kwargs['chord']*radius

    # The algorithm works with radius-normalized chord lengths, but the BEMT takes real chord lengths, so multiply by R
    omega = xn[0]
    twist0 = xn[1]
    dtwist = xn[2:]

    twist = calc_twist_dist(twist0, dtwist)

    f = 1000.
    fail = 0
    g = [1.0]

    prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy,
                               airfoils=airfoils, Cl_tables=Cl_tables, Cd_tables=Cd_tables)

    try:
        dT, P = bemt.bemt_axial(prop, pitch, omega, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                tip_loss=tip_loss, mach_corr=mach_corr)
    except FloatingPointError:
        fail = 1
        return f, g, fail
    except IndexError:
        return 10000, g, fail

    f = P
    print "Power = " + str(f)
    print "Thrust = %s" % str(sum(dT))

    # Evaluate thrust constraint. Target thrust must be less than computed thrust
    g[0] = target_thrust - sum(dT)

    return f, g, fail


def get_geocons(c, cval_max, R):
    geocons = [0.0] * (2*len(c))
    i = 0
    # Check if all chord values are greater than zero
    for cval in c:
        geocons[i] = -cval
        i += 1
    # Check if all c/R values are below cval_max
    for cval in c:
        geocons[i] = cval/R - cval_max
        i += 1
    return geocons


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


def optimize_chord(**k):

    omega = k['omega']
    omega_lower = k['omega_lower']
    omega_upper = k['omega_upper']
    chord0 = k['chord0']
    chord0_lower = k['chord0_lower']
    chord0_upper = k['chord0_upper']
    n_elements = k['n_elements']
    dchord = k['dchord']
    dchord_lower = k['dchord_lower']
    dchord_upper = k['dchord_upper']

    opt_prob_ft = Optimization('Rotor in Hover w/ Fixed Twist', objfun_optimize_chord)
    opt_prob_ft.addVar('omega', 'c', value=omega, lower=omega_lower, upper=omega_upper)
    opt_prob_ft.addVar('chord0', 'c', value=chord0, lower=chord0_lower, upper=chord0_upper)
    opt_prob_ft.addVarGroup('dchord', n_elements-1, 'c', value=dchord, lower=dchord_lower, upper=dchord_upper)
    opt_prob_ft.addObj('f')
    opt_prob_ft.addCon('thrust', 'i')
    opt_prob_ft.addConGroup('c_lower', n_elements, 'i')
    opt_prob_ft.addConGroup('c_upper', n_elements, 'i')
    #print opt_prob_ft

    n_blades = k['n_blades']
    root_cutout = k['root_cutout']
    radius = k['radius']
    dy = k['dy']
    dr = k['dr']
    y = k['y']
    r = k['r']
    pitch = k['pitch']
    airfoils = k['airfoils']
    thrust = k['thrust']
    max_chord = k['max_chord']
    twist = k['twist']
    allowable_Re = k['allowable_Re']
    Cl_tables = k['Cl_tables']
    Cd_tables = k['Cd_tables']
    Cl_funs = k['Cl_funs']
    Cd_funs = k['Cd_funs']

    # Routine for optimizing chord with constant twist
    slsqp1 = SLSQP()
    slsqp1.setOption('IPRINT', 1)
    slsqp1.setOption('MAXIT', 200)
    slsqp1.setOption('ACC', 1e-7)
    [fstr, xstr, inform] = slsqp1(opt_prob_ft, sens_type='FD', n_blades=n_blades, n_elements=n_elements,
                                  root_cutout=root_cutout, radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch,
                                  airfoils=airfoils, thrust=thrust, max_chord=max_chord, tip_loss=False,
                                  mach_corr=False, omega=omega, twist=twist, allowable_Re=allowable_Re,
                                  Cl_tables=Cl_tables, Cd_tables=Cd_tables, Cl_funs=Cl_funs, Cd_funs=Cd_funs)
    # #print opt_prob_ft.solution(0)

    # conmin = CONMIN()
    # conmin.setOption('IPRINT', 1)
    # [fstr, xstr, inform] = conmin(opt_prob_ft, sens_type='FD', n_blades=n_blades, n_elements=n_elements,
    #                               root_cutout=root_cutout, radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch,
    #                               airfoils=airfoils, thrust=thrust, max_chord=max_chord, tip_loss=False,
    #                               mach_corr=False, omega=omega, twist=twist, allowable_Re=allowable_Re,
    #                               Cl_tables=Cl_tables, Cd_tables=Cd_tables, Cl_funs=Cl_funs, Cd_funs=Cd_funs)
    # print opt_prob.solution(0)

    return fstr, xstr


def optimize_twist(**k):

    omega = k['omega']
    omega_lower = k['omega_lower']
    omega_upper = k['omega_upper']
    twist0 = k['twist0']
    twist0_lower = k['twist0_lower']
    twist0_upper = k['twist0_upper']
    n_elements = k['n_elements']
    dtwist = k['dtwist']
    dtwist_lower = k['dtwist_lower']
    dtwist_upper = k['dtwist_upper']

    opt_prob_fc = Optimization('Rotor in Hover w/ Fixed Chord', objfun_optimize_twist)
    opt_prob_fc.addVar('omega', 'c', value=omega, lower=omega_lower, upper=omega_upper)
    opt_prob_fc.addVar('twist0', 'c', value=twist0, lower=twist0_lower, upper=twist0_upper)
    opt_prob_fc.addVarGroup('dtwist', n_elements-1, 'c', value=dtwist, lower=dtwist_lower, upper=dtwist_upper)
    opt_prob_fc.addObj('f')
    opt_prob_fc.addCon('thrust', 'i')
    #print opt_prob_fc

    n_blades = k['n_blades']
    root_cutout = k['root_cutout']
    radius = k['radius']
    dy = k['dy']
    dr = k['dr']
    y = k['y']
    r = k['r']
    pitch = k['pitch']
    airfoils = k['airfoils']
    thrust = k['thrust']
    chord = k['chord']
    allowable_Re = k['allowable_Re']
    Cl_tables = k['Cl_tables']
    Cd_tables = k['Cd_tables']
    Cl_funs = k['Cl_funs']
    Cd_funs = k['Cd_funs']

    # Routine for optimizing twist with a constant chord
    slsqp2 = SLSQP()
    slsqp2.setOption('IPRINT', 1)
    slsqp2.setOption('MAXIT', 200)
    slsqp2.setOption('ACC', 1e-7)
    [fstr, xstr, inform] = slsqp2(opt_prob_fc, sens_type='FD', n_blades=n_blades, n_elements=n_elements,
                                  root_cutout=root_cutout, radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch,
                                  airfoils=airfoils, thrust=thrust, tip_loss=False, mach_corr=False, omega=omega,
                                  chord=chord, allowable_Re=allowable_Re, Cl_tables=Cl_tables, Cd_tables=Cd_tables,
                                  Cl_funs=Cl_funs, Cd_funs=Cd_funs)
    #print opt_prob_fc.solution(0)

    # nsga2 = NSGA2()
    # nsga2.setOption('PrintOut', 2)
    # nsga2.setOption('PopSize', 5000)
    # nsga2.setOption('maxGen', 1000)
    # nsga2.setOption('pCross_real', 0.85)
    # #nsga2.setOption('xinit', 1)
    # #nsga2.setOption('seed', 0.541)
    # nsga2(opt_prob, n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout, radius=radius, dy=dy,
    #       dr=dr, y=y, r=r, pitch=pitch, airfoils=airfoils, thrust=thrust, max_chord=max_chord, tip_loss=False,
    #       mach_corr=False)
    # print opt_prob.solution(0)

    return fstr, xstr


def main():
    n_blades = 2
    n_elements = 5
    radius = unit_conversion.in2m(9.0)/2
    root_cutout = 0.1 * radius
    dy = float(radius-root_cutout)/n_elements
    dr = float(1)/n_elements
    y = root_cutout + dy*np.arange(1, n_elements+1)
    r = y/radius
    pitch = 0.0
    airfoils = (('SDA1075_494p', 0.0, 1.0),)
    allowable_Re = []
    # allowable_Re = [1000000., 500000., 250000., 100000., 90000.]
    allowable_Re = [1000000.]
    thrust = 10.75
    max_chord = 100

    Cl_tables = {}
    Cd_tables = {}
    # Get lookup tables
    if any(airfoil[0] != 'simple' for airfoil in airfoils):
        for airfoil in airfoils:
            Cl_table, Cd_table = aero_coeffs.create_Cl_Cd_table(airfoil[0])

            Cl_tables[airfoil[0]] = Cl_table
            Cd_tables[airfoil[0]] = Cd_table

    # Create list of Cl functions. One for each Reynolds number. Cl_tables (and Cd_tables) will be empty for the
    # 'simple' case, therefore this will be skipped for the simple case. For the full table lookup case this will be
    # skipped because allowable_Re will be empty.
    Cl_funs = {}
    Cd_funs = {}
    if Cl_tables and allowable_Re:
        Cl_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cl_fun(Re, Cl_tables[airfoils[0][0]]) for Re in allowable_Re]))
        Cd_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cd_fun(Re, Cd_tables[airfoils[0][0]]) for Re in allowable_Re]))

    ###########################################
    # Set design variable bounds
    ###########################################
    # The DA4002 blade
    chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
                      0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
    chord = np.array([chord[i] for i in [0, 4, 8, 12, 16]])
    chord = chord[-1] / r
    twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
                      18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
    twist = np.array([twist[i] for i in [0, 4, 8, 12, 16]])

    chord0 = chord[0]
    twist0 = twist[0]
    dchord = np.array([chord[i+1]-chord[i] for i in xrange(len(chord)-1)])
    dtwist = np.array([twist[i+1]-twist[i] for i in xrange(len(twist)-1)])

    omega = 5943.0 * 2*np.pi/60
    omega_upper = 7000 * 2*np.pi/60
    omega_lower = 4500. * 2*np.pi/60

    # # Twist at the hub must be less than or equal to arcsin(hub_height/hub_diameter), approx 23 degrees
    twist0_lower = 0.0 * 2 * np.pi / 360
    twist0_upper = 60.0 * 2 * np.pi / 360
    # Chord values at the hub must be less than or equal to the diameter of the hub
    chord0_upper = 10
    chord0_lower = 0.05

    # dtwist_start = 0.0 * 2 * np.pi / 360
    dtwist_lower = -10.0 * 2 * np.pi / 360
    dtwist_upper = 10.0 * 2 * np.pi / 360
    dchord_lower = -0.8
    dchord_upper = 0.8

    twist_start = twist
    chord_start = chord

    iteration = 0
    fit_old = 1000.
    converged = False
    while not converged:
        iteration += 1
        fstr, xstr = optimize_twist(omega=omega, omega_lower=omega_lower, omega_upper=omega_upper, twist0=twist0,
                                    twist0_lower=twist0_lower, twist0_upper=twist0_upper, n_elements=n_elements,
                                    dtwist=dtwist, dtwist_lower=dtwist_lower, dtwist_upper=dtwist_upper, n_blades=n_blades,
                                    root_cutout=root_cutout, radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch,
                                    airfoils=airfoils, thrust=thrust, chord=chord, allowable_Re=allowable_Re,
                                    Cl_tables=Cl_tables, Cd_tables=Cd_tables, Cl_funs=Cl_funs, Cd_funs=Cd_funs)

        omega = xstr[0]
        twist = calc_twist_dist(xstr[1], xstr[2:])

        fstr, xstr = optimize_chord(omega=omega, omega_lower=omega_lower, omega_upper=omega_upper, chord0=chord0,
                                    chord0_lower=chord0_lower, chord0_upper=chord0_upper, n_elements=n_elements,
                                    dchord=dchord, dchord_lower=dchord_lower, dchord_upper=dchord_upper, n_blades=n_blades,
                                    root_cutout=root_cutout, radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch,
                                    airfoils=airfoils, thrust=thrust, max_chord=max_chord, twist=twist,
                                    allowable_Re=allowable_Re, Cl_tables=Cl_tables, Cd_tables=Cd_tables, Cl_funs=Cl_funs,
                                    Cd_funs=Cd_funs)

        omega = xstr[0]
        chord = calc_chord_dist(xstr[1], xstr[2:])

        converged = abs((fit_old - fstr) / fstr) < 0.00005
        fit_old = fstr

    print "omega = " + str(omega*60/2/np.pi)
    print "chord = " + str(chord)
    print "twist = " + str(twist)
    print "iterations = " + str(iteration)

    def get_performance(o, c, t):
        chord_meters = c * radius
        prop = propeller.Propeller(t, chord_meters, radius, n_blades, r, y, dr, dy, airfoils=airfoils,
                                   Cl_tables=Cl_tables, Cd_tables=Cd_tables)

        return bemt.bemt_axial(prop, pitch, o, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                               tip_loss=False, mach_corr=False, output='long')

    perf = get_performance(omega, chord, twist)
    print "Thrust of optimized = " + str(sum(perf[0]))
    print "Power of optimized = " + str(sum(perf[1]))

    plt.figure(1)
    plt.plot(r, chord_start, '-b')
    plt.plot(r, chord, '-r')
    plt.xlabel('radial location')
    plt.ylabel('c/R')
    plt.legend(['start', 'opt'])

    plt.figure(2)
    plt.plot(r, twist_start*180/np.pi, '-b')
    plt.plot(r, twist*180/np.pi, '-r')
    plt.xlabel('radial location')
    plt.ylabel('twist')
    plt.legend(['start', 'opt'])

    plt.show()

if __name__ == "__main__":
    main()