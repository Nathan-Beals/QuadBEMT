import unit_conversion
import propeller
import bemt
from pyOpt import Optimization
from pyOpt import NSGA2
from pyOpt import SLSQP
import numpy as np
import aero_coeffs
import matplotlib.pyplot as plt


def counted(f):
    def wrapped(*args, **kwargs):
        print wrapped.calls
        wrapped.calls += 1
        return f(*args, **kwargs)
    wrapped.calls = 0
    return wrapped


@counted
def objfun(xn, **kwargs):

    #print 'objfun called'

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
    Cl_tables = kwargs['Cl_tables']
    Cd_tables = kwargs['Cd_tables']
    Cl_funs = kwargs['Cl_funs']
    Cd_funs = kwargs['Cd_funs']
    allowable_Re = kwargs['allowable_Re']
    opt_method = kwargs['opt_method']
    omega = xn[0]
    twist0 = xn[1]
    chord0 = xn[2]*radius
    dtwist = np.array(xn[3:len(r)+2])
    dchord = np.array([x*radius for x in xn[len(r)+2:2*len(r)+2]])

    twist = calc_twist_dist(twist0, dtwist)
    chord = calc_chord_dist(chord0, dchord)

    # int_twist = np.array(xn[:len(r)])
    # int_chord = np.array([x*radius for x in xn[len(r):]])

    f = 1000.
    fail = 0
    g = [1.0] * (2*len(r)+1)

    # Calculate geometric constraint values. If a genetic algorithm is used we can fail the case immediately if there
    # are any violations. If a gradient-based algorithm is used this will cause the gradient calculation to fail so the
    # constraints must be checked normally by the optimizer.
    g[1:] = get_geocons(chord, cval_max, radius)
    if opt_method == 'nsga2':
        if any(geocon > 0.0 for geocon in g[1:]):
            print "geocons violated"
            return f, g, fail

    prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy, airfoils=airfoils, Cl_tables=Cl_tables,
                               Cd_tables=Cd_tables)

    try:
        dT, P = bemt.bemt_axial(prop, pitch, omega, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                tip_loss=tip_loss, mach_corr=mach_corr)
    except (FloatingPointError, IndexError):
        fail = 1
        return f, g, fail

    f = P
    print f
    print "Thrust = %s" % str(sum(dT))

    # Evaluate thrust constraint. Target thrust must be less than computed thrust
    g[0] = target_thrust - sum(dT)

    return f, g, fail


def get_geocons(c, cval_max, R):
    geocons = [0.0] * (2*len(c))
    i = 0
    for cval in c:
        geocons[i] = -cval
        i += 1
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


def main():
    ###########################################
    # Define some values
    ###########################################
    n_blades = 2
    n_elements = 10
    radius = unit_conversion.in2m(9.0)/2
    root_cutout = 0.1 * radius
    dy = float(radius-root_cutout)/n_elements
    dr = float(1)/n_elements
    y = root_cutout + dy*np.arange(1, n_elements+1)
    r = y/radius
    pitch = 0.0
    airfoils = (('SDA1075_494p', 0.0, 1.0),)
    #allowable_Re = [1000000., 500000., 250000., 100000., 90000., 80000.]
    allowable_Re = [1000000.]
    thrust = 6.31
    max_chord = 0.4

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
    # chord_base = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787,
    #                         0.1786, 0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
    # twist_base = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009,
    #                         19.814, 18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
    # dtwist_base = np.array([twist_base[i+1]-twist_base[i] for i in xrange(len(twist_base)-1)])
    # dchord_base = np.array([chord_base[i+1]-chord_base[i] for i in xrange(len(chord_base)-1)])
    # twist0_base = twist_base[0]
    # chord0_base = chord_base[0]

    # dtwist_start = dtwist_base
    # dchord_start = dchord_base
    # twist0_start = twist0_base
    # chord0_start = chord0_base

    dtwist_start = np.zeros(n_elements-1)
    dchord_start = np.zeros(n_elements-1)
    twist0_start = 0.0
    chord0_start = 0.0


    omega_lower = 4500.0 * 2*np.pi/60
    omega_upper = 7000.0 * 2*np.pi/60
    omega_start = 5943.0 * 2*np.pi/60

    twist0_lower = 0.0 * 2 * np.pi / 360
    twist0_upper = 60.0 * 2 * np.pi / 360

    chord0_upper = 0.1198
    chord0_lower = 0.05

    dtwist_lower = -10.0 * 2 * np.pi / 360
    dtwist_upper = 10.0 * 2 * np.pi / 360
    dchord_lower = -0.1
    dchord_upper = 0.1

    opt_prob = Optimization('Rotor in Hover', objfun)
    opt_prob.addVar('omega', 'c', value=omega_start, lower=omega_lower, upper=omega_upper)
    opt_prob.addVar('twist0', 'c', value=twist0_start, lower=twist0_lower, upper=twist0_upper)
    opt_prob.addVar('chord0', 'c', value=chord0_start, lower=chord0_lower, upper=chord0_upper)
    opt_prob.addVarGroup('dtwist', n_elements-1, 'c', value=dtwist_start, lower=dtwist_lower, upper=dtwist_upper)
    opt_prob.addVarGroup('dchord', n_elements-1, 'c', value=dchord_start, lower=dchord_lower, upper=dchord_upper)
    opt_prob.addObj('f')
    opt_prob.addCon('thrust', 'i')
    opt_prob.addConGroup('c_lower', n_elements, 'i')
    opt_prob.addConGroup('c_upper', n_elements, 'i')
    print opt_prob

    opt_method = 'nsga2'
    nsga2 = NSGA2()
    nsga2.setOption('PrintOut', 2)
    nsga2.setOption('PopSize', 5000)
    nsga2.setOption('maxGen', 500)
    nsga2.setOption('pCross_real', 0.85)
    fstr, xstr, inform = nsga2(opt_prob, n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout,
                               radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch, airfoils=airfoils, thrust=thrust,
                               max_chord=max_chord, tip_loss=False, mach_corr=False, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                               Cl_tables=Cl_tables, Cd_tables=Cd_tables, allowable_Re=allowable_Re,
                               opt_method=opt_method)
    print opt_prob.solution(0)

    # opt_method = 'slsqp'
    # slsqp = SLSQP()
    # slsqp.setOption('IPRINT', 1)
    # slsqp.setOption('MAXIT', 500)
    # slsqp.setOption('ACC', 1e-8)
    # fstr, xstr, inform = slsqp(opt_prob, sens_type='FD', n_blades=n_blades, n_elements=n_elements,
    #                            root_cutout=root_cutout, radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch,
    #                            airfoils=airfoils, thrust=thrust, max_chord=max_chord, tip_loss=True, mach_corr=True,
    #                            Cl_funs=Cl_funs, Cd_funs=Cd_funs, Cl_tables=Cl_tables, Cd_tables=Cd_tables,
    #                            allowable_Re=allowable_Re, opt_method=opt_method, c0_max=c0_max)
    #print opt_prob.solution(0)

    def get_performance(o, c, t):
        chord_meters = c * radius
        prop = propeller.Propeller(t, chord_meters, radius, n_blades, r, y, dr, dy, airfoils=airfoils,
                                   Cl_tables=Cl_tables, Cd_tables=Cd_tables)

        return bemt.bemt_axial(prop, pitch, o, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                               tip_loss=False, mach_corr=False, output='long')

    omega = xstr[0]
    twist0 = xstr[1]
    chord0 = xstr[2]
    dtwist = xstr[3:3+len(r)-1]
    dchord = xstr[3+len(r)-1:]

    twist = calc_twist_dist(twist0, dtwist)
    chord = calc_chord_dist(chord0, dchord)

    print "Opt chord = " + str(chord)
    print "Opt twist = " + str(twist)

    # twist_base = calc_twist_dist(twist0_base, dtwist_base)
    # chord_base = calc_chord_dist(chord0_base, dchord_base)

    perf_opt = get_performance(omega, chord, twist)
    #perf_base = get_performance(omega_start, chord_base, twist_base)
    print "Omega optimized = " + str(omega*60/2/np.pi)
    print "Thrust of optimized = " + str(sum(perf_opt[0]))
    print "Power of optimized = " + str(sum(perf_opt[1]))
    print "Omega base = " + str(omega_start*60/2/np.pi)
    # print "Thrust of base = " + str(sum(perf_base[0]))
    # print "Power of base = " + str(sum(perf_base[1]))

    plt.figure(1)
    #plt.plot(r, chord_base, '-b')
    plt.plot(r, chord, '-r')
    plt.xlabel('radial location')
    plt.ylabel('c/R')
    plt.legend(['base', 'opt'])

    plt.figure(2)
    #plt.plot(r, twist_base*180/np.pi, '-b')
    plt.plot(r, twist*180/np.pi, '-r')
    plt.xlabel('radial location')
    plt.ylabel('twist')
    plt.legend(['base', 'opt'])

    plt.show()

if __name__ == "__main__":
    main()