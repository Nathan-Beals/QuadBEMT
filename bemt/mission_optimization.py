import unit_conversion
import propeller
import quadrotor
import bemt
import trim
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
    #print "objfun entered"
    radius = kwargs['radius']
    r = kwargs['r']
    y = kwargs['y']
    dr = kwargs['dr']
    dy = kwargs['dy']
    n_blades = kwargs['n_blades']
    airfoils = kwargs['airfoils']
    pitch = kwargs['pitch']
    vehicle_weight = kwargs['vehicle_weight']
    max_chord = kwargs['max_chord']
    tip_loss = kwargs['tip_loss']
    mach_corr = kwargs['mach_corr']
    Cl_tables = kwargs['Cl_tables']
    Cd_tables = kwargs['Cd_tables']
    Cl_funs = kwargs['Cl_funs']
    Cd_funs = kwargs['Cd_funs']
    lift_curve_info_dict = kwargs['lift_curve_info_dict']
    allowable_Re = kwargs['allowable_Re']
    opt_method = kwargs['opt_method']
    alt = kwargs['alt']
    v_inf = kwargs['v_inf']
    alpha0 = kwargs['alpha0']
    n_azi_elements = kwargs['n_azi_elements']
    mission_time = kwargs['mission_time']
    omega_h = xn[0]
    twist0 = xn[1]
    chord0 = xn[2]*radius   # Convert to meters from c/R before we use in calculations
    dtwist = np.array(xn[3:len(r)+2])
    dchord = np.array([x*radius for x in xn[len(r)+2:2*len(r)+2]])

    twist = calc_twist_dist(twist0, dtwist)
    chord = calc_chord_dist(chord0, dchord)

    f = 10000000.
    fail = 0
    g = [1.0] * (2*len(r)+1)

    # Calculate geometric constraint values. If a genetic algorithm is used we can fail the case immediately if there
    # are any violations. If a gradient-based algorithm is used this will cause the gradient calculation to fail so the
    # constraints must be checked normally by the optimizer.
    g[1:] = get_geocons(chord, max_chord, radius)
    if opt_method == 'nograd':
        if any(geocon > 0.0 for geocon in g[1:]):
            print "geocons violated"
            return f, g, fail

    prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy, airfoils=airfoils, Cl_tables=Cl_tables,
                               Cd_tables=Cd_tables)

    quad = quadrotor.Quadrotor(prop, vehicle_weight)

    try:
        dT_h, P_h = bemt.bemt_axial(prop, pitch, omega_h, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                    tip_loss=tip_loss, mach_corr=mach_corr, alt=alt)
        if P_h < 0:
            print "hover power negative"
            fail = 1
            return f, g, fail
    except FloatingPointError:
        print "Floating point error in axial BEMT"
        fail = 1
        return f, g, fail
    except IndexError:
        print "Index error in axial BEMT"
        fail = 1
        return f, g, fail

    try:
        trim0 = [alpha0, omega_h]   # Use alpha0 (supplied by user) and the hover omega as initial guesses for trim
        ff_kwargs = {'propeller': prop, 'pitch': pitch, 'n_azi_elements': n_azi_elements, 'allowable_Re': allowable_Re,
                     'Cl_funs': Cl_funs, 'Cd_funs': Cd_funs, 'tip_loss': tip_loss, 'mach_corr': mach_corr, 'alt': alt,
                     'lift_curve_info_dict': lift_curve_info_dict}

        alpha_trim, omega_trim, converged = trim.trim(quad, v_inf, trim0, ff_kwargs)
        if not converged:
            fail = 1
            return f, g, fail

        T_trim, H_trim, P_trim = bemt.bemt_forward_flight(quad, pitch, omega_trim, alpha_trim, v_inf, n_azi_elements,
                                                          alt=alt, tip_loss=tip_loss, mach_corr=mach_corr,
                                                          allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                                          lift_curve_info_dict=lift_curve_info_dict)
        if P_trim < 0:
            print "ff power negative"
            fail = 1
            return f, g, fail
    except Exception as e:
        print "{} in ff trim".format(type(e).__name__)
        fail = 1
        return f, g, fail

    # Find total energy mission_times = [time_in_hover, time_in_ff] in seconds
    energy = P_h * mission_time[0] + P_trim * mission_time[1]

    f = energy
    print "total energy = " + str(f)
    print "Thrust hover = %s" % str(sum(dT_h))
    print "(alpha, omega) = (%f, %f)" % (alpha_trim, omega_trim)

    # Evaluate performance constraints.
    g[0] = vehicle_weight/4 - sum(dT_h)

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
    v_inf = 4.     # m/s
    alpha0 = 0.0454  # Starting guess for trimmed alpha in radians
    n_azi_elements = 5

    # Mission times
    time_in_hover = 5. * 60     # Time in seconds
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
    if Cl_tables and allowable_Re:
        Cl_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cl_fun(Re, Cl_tables[airfoils[0][0]], Clmax[airfoils[0][0]][Re]) for Re in allowable_Re]))
        Cd_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cd_fun(Re, Cd_tables[airfoils[0][0]]) for Re in allowable_Re]))
        lift_curve_info_dict = aero_coeffs.create_liftCurveInfoDict(allowable_Re, Cl_tables[airfoils[0][0]])

    ###########################################
    # Set design variable bounds
    ###########################################
    omega_start = 4250. * 2*np.pi/60
    # These are c/R values for the DA4002 propeller given at the UIUC propeller database
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
    print "chord0_start = " + str(chord0_start)

    # chord = np.array([8.92386048e-02, 1.73000845e-01, 2.70523039e-01, 2.71542807e-01, 2.78749355e-01, 2.36866151e-01,
    #                   2.04103526e-01, 1.37456074e-01, 8.68094589e-02, 1.05601135e-04])
    # twist = np.array([0.00161645, 0.15105685, 0.28791442, 0.31577392, 0.28644651, 0.27418749, 0.24854514, 0.21812646,
    #                   0.19802027, 0.14972058])
    # omega_start = 3184.41320387 * 2*np.pi/60
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

    omega_lower = 2000 * 2*np.pi/60
    omega_upper = 8000.0 * 2*np.pi/60

    twist0_lower = 0.0 * 2 * np.pi / 360
    twist0_upper = 80. * 2 * np.pi / 360

    chord0_upper = 0.1198
    chord0_lower = 0.05

    dtwist_lower = -10.0 * 2 * np.pi / 360
    dtwist_upper = 10.0 * 2 * np.pi / 360
    dchord_lower = -0.1
    dchord_upper = 0.1

    opt_prob = Optimization('Mission Simulator', objfun)
    opt_prob.addVar('omega_h', 'c', value=omega_start, lower=omega_lower, upper=omega_upper)
    opt_prob.addVar('twist0', 'c', value=twist0_start, lower=twist0_lower, upper=twist0_upper)
    opt_prob.addVar('chord0', 'c', value=chord0_start, lower=chord0_lower, upper=chord0_upper)
    opt_prob.addVarGroup('dtwist', n_elements-1, 'c', value=dtwist_start, lower=dtwist_lower, upper=dtwist_upper)
    opt_prob.addVarGroup('dchord', n_elements-1, 'c', value=dchord_start, lower=dchord_lower, upper=dchord_upper)
    opt_prob.addObj('f')
    opt_prob.addCon('thrust', 'i')
    opt_prob.addConGroup('c_lower', n_elements, 'i')
    opt_prob.addConGroup('c_upper', n_elements, 'i')
    print opt_prob

    pop_size = 5000
    max_gen = 8
    opt_method = 'nograd'
    nsga2 = NSGA2()
    nsga2.setOption('PrintOut', 2)
    nsga2.setOption('PopSize', pop_size)
    nsga2.setOption('maxGen', max_gen)
    nsga2.setOption('pCross_real', 0.85)
    nsga2.setOption('xinit', 1)
    fstr, xstr, inform = nsga2(opt_prob, n_blades=n_blades, radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch,
                               airfoils=airfoils, vehicle_weight=vehicle_weight, max_chord=max_chord, tip_loss=tip_loss,
                               mach_corr=mach_corr, Cl_funs=Cl_funs, Cd_funs=Cd_funs, Cl_tables=Cl_tables,
                               Cd_tables=Cd_tables, allowable_Re=allowable_Re, opt_method=opt_method, alt=alt,
                               v_inf=v_inf, alpha0=alpha0, mission_time=mission_time, n_azi_elements=n_azi_elements,
                               pop_size=pop_size, max_gen=max_gen, lift_curve_info_dict=lift_curve_info_dict)
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
        quad = quadrotor.Quadrotor(prop, vehicle_weight)

        ff_kwargs = {'propeller': prop, 'pitch': pitch, 'n_azi_elements': n_azi_elements, 'allowable_Re': allowable_Re,
                     'Cl_funs': Cl_funs, 'Cd_funs': Cd_funs, 'tip_loss': tip_loss, 'mach_corr': mach_corr, 'alt': alt,
                     'lift_curve_info_dict': lift_curve_info_dict}
        trim0 = np.array([alpha0, o])
        alpha_trim, omega_trim, converged = trim.trim(quad, v_inf, trim0, ff_kwargs)
        T_ff, H_ff, P_ff = bemt.bemt_forward_flight(quad, pitch, omega_trim, alpha_trim, v_inf, n_azi_elements, alt=alt,
                                                    tip_loss=tip_loss, mach_corr=mach_corr, allowable_Re=allowable_Re,
                                                    Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                                    lift_curve_info_dict=lift_curve_info_dict)

        dT_h, P_h = bemt.bemt_axial(prop, pitch, o, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                    tip_loss=tip_loss, mach_corr=mach_corr, alt=alt)
        return sum(dT_h), P_h, T_ff, P_ff, alpha_trim, omega_trim

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
    print "Hover Thrust of optimized = " + str(perf_opt[0])
    print "Hover Power of optimized = " + str(perf_opt[1])
    print "FF Thrust of optimized = " + str(perf_opt[2])
    print "FF Power of optimized = " + str(perf_opt[3])
    print "Trim (alpha, omega) = (%f, %f)" % (perf_opt[4], perf_opt[5])
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