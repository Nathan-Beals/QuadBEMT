import unit_conversion
import propeller
import bemt
from pyOpt import Optimization
from pyOpt import NSGA2
from pyOpt import SLSQP
from pyOpt import ALPSO
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
    alt = kwargs['alt']
    omega = xn[0]
    twist0 = xn[1]
    chord0 = xn[2]*radius
    dtwist = np.array(xn[3:len(r)+2])
    dchord = np.array([x*radius for x in xn[len(r)+2:2*len(r)+2]])

    twist = calc_twist_dist(twist0, dtwist)
    chord = calc_chord_dist(chord0, dchord)

    f = 1000.
    fail = 0
    g = [1.0] * (2*len(r)+1)

    # Calculate geometric constraint values. If a genetic algorithm is used we can fail the case immediately if there
    # are any violations. If a gradient-based algorithm is used this will cause the gradient calculation to fail so the
    # constraints must be checked normally by the optimizer.
    g[1:] = get_geocons(chord, cval_max, radius)
    if opt_method == 'nograd':
        if any(geocon > 0.0 for geocon in g[1:]):
            print "geocons violated"
            return f, g, fail

    prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy, airfoils=airfoils, Cl_tables=Cl_tables,
                               Cd_tables=Cd_tables)

    try:
        dT, P, FM = bemt.bemt_axial(prop, pitch, omega, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                    tip_loss=tip_loss, mach_corr=mach_corr, alt=alt)
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
    #radius = 0.1930
    root_cutout = 0.1 * radius
    dy = float(radius-root_cutout)/n_elements
    dr = float(1)/n_elements
    y = root_cutout + dy*np.arange(1, n_elements+1)
    r = y/radius
    pitch = 0.0
    airfoils = (('SDA1075_494p', 0.0, 1.0),)
    allowable_Re = []
    allowable_Re = [1000000., 500000., 250000., 100000., 90000., 80000., 70000., 60000., 50000., 40000., 30000., 20000., 10000.]
    thrust = 4.61
    max_chord = 0.4
    alt = 0
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

    ###########################################
    # Set design variable bounds
    ###########################################
    omega_start = 5943.0 * 2*np.pi/60
    chord_base = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787,
                           0.1786, 0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
    chord_base = np.array([chord_base[i] for i in [0, 2, 4, 6, 8, 10, 12, 14, 15, 17]])
    twist_base = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009,
                           19.814, 18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
    twist_base = np.array([twist_base[i] for i in [0, 2, 4, 6, 8, 10, 12, 14, 15, 17]])
    dtwist_base = np.array([twist_base[i+1]-twist_base[i] for i in xrange(len(twist_base)-1)])
    dchord_base = np.array([chord_base[i+1]-chord_base[i] for i in xrange(len(chord_base)-1)])
    twist0_base = twist_base[0]
    chord0_base = chord_base[0]

    chord_start = chord_base
    twist_start = twist_base
    dtwist_start = dtwist_base
    dchord_start = dchord_base
    twist0_start = twist0_base
    chord0_start = chord0_base

    # chord = np.array([0.11958375, 0.21820237, 0.2761888, 0.29582408, 0.2670838, 0.25332996, 0.21333354, 0.17139233,
    #                   0.09985169, 0.00114413])
    # twist = np.array([0.34160335, 0.39718607, 0.37649986, 0.35195434, 0.26725317, 0.25970378, 0.24134294, 0.19232958,
    #                   0.17713771, 0.02991615])
    # omega_start = 6574.95092927 * 2*np.pi/60
    # chord_start = chord
    # twist_start = twist
    # dchord_start = np.array([chord[i+1]-chord[i] for i in xrange(len(chord)-1)])
    # dtwist_start = np.array([twist[i+1]-twist[i] for i in xrange(len(twist)-1)])
    # twist0_start = twist[0]
    # chord0_start = chord[0]

    ## Initialize everything to zeros
    # dtwist_start = np.zeros(n_elements-1)
    # dchord_start = np.zeros(n_elements-1)
    # twist0_start = 0.0
    # chord0_start = 0.0

    omega_lower = 4000.0 * 2*np.pi/60
    omega_upper = 7000.0 * 2*np.pi/60

    twist0_lower = 0.0 * 2 * np.pi / 360
    twist0_upper = 50.0 * 2 * np.pi / 360

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

    opt_method = 'nograd'
    nsga2 = NSGA2()
    nsga2.setOption('PrintOut', 2)
    nsga2.setOption('PopSize', 5000)
    nsga2.setOption('maxGen', 500)
    nsga2.setOption('pCross_real', 0.85)
    nsga2.setOption('xinit', 1)
    fstr, xstr, inform = nsga2(opt_prob, n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout,
                               radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch, airfoils=airfoils, thrust=thrust,
                               max_chord=max_chord, tip_loss=tip_loss, mach_corr=mach_corr, Cl_funs=Cl_funs,
                               Cd_funs=Cd_funs, Cl_tables=Cl_tables, Cd_tables=Cd_tables, allowable_Re=allowable_Re,
                               opt_method=opt_method, alt=alt)
    print opt_prob.solution(0)

    # opt_method = 'nograd'
    # xstart_alpso = np.concatenate((np.array([omega_start, twist0_start, chord0_start]), dtwist_start, dchord_start))
    # alpso = ALPSO()
    # alpso.setOption('xinit', 0)
    # alpso.setOption('SwarmSize', 200)
    # alpso.setOption('maxOuterIter', 100)
    # alpso.setOption('stopCriteria', 0)
    # fstr, xstr, inform = alpso(opt_prob, xstart=xstart_alpso,  n_blades=n_blades, n_elements=n_elements,
    #                            root_cutout=root_cutout, radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch,
    #                            airfoils=airfoils, thrust=thrust, max_chord=max_chord, tip_loss=tip_loss,
    #                            mach_corr=mach_corr, Cl_funs=Cl_funs, Cd_funs=Cd_funs, Cl_tables=Cl_tables,
    #                            Cd_tables=Cd_tables, allowable_Re=allowable_Re, opt_method=opt_method)
    # print opt_prob.solution(0)

    # opt_method = 'grad'
    # slsqp = SLSQP()
    # slsqp.setOption('IPRINT', 1)
    # slsqp.setOption('MAXIT', 1000)
    # slsqp.setOption('ACC', 1e-7)
    # fstr, xstr, inform = slsqp(opt_prob, sens_type='FD', n_blades=n_blades, n_elements=n_elements,
    #                            root_cutout=root_cutout, radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch,
    #                            airfoils=airfoils, thrust=thrust, max_chord=max_chord,
    #                            tip_loss=tip_loss, mach_corr=mach_corr, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
    #                            Cl_tables=Cl_tables, Cd_tables=Cd_tables, allowable_Re=allowable_Re,
    #                            opt_method=opt_method, alt=alt)
    # print opt_prob.solution(0)

    def get_performance(o, c, t):
        chord_meters = c * radius
        prop = propeller.Propeller(t, chord_meters, radius, n_blades, r, y, dr, dy, airfoils=airfoils,
                                   Cl_tables=Cl_tables, Cd_tables=Cd_tables)

        return bemt.bemt_axial(prop, pitch, o, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                               tip_loss=tip_loss, mach_corr=mach_corr, output='long', alt=alt)

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
    # print "Omega base = " + str(omega_start*60/2/np.pi)
    # print "Thrust of base = " + str(sum(perf_base[0]))
    # print "Power of base = " + str(sum(perf_base[1]))

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