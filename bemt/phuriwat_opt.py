from pyOpt import Optimization
from pyOpt import CONMIN
from pyOpt import SLSQP
import numpy as np
import propeller
import bemt
import unit_conversion


def counted(f):
    def wrapped(*args, **kwargs):
        print wrapped.calls
        wrapped.calls += 1
        return f(*args, **kwargs)
    wrapped.calls = 0
    return wrapped


@counted
def objfun_const_twist(xn, **kwargs):

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
    #omega = kwargs['omega']
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
    if any(geocon > 0.0 for geocon in g[1:]):
        print "geocons violated"
        print g[1:]
        return f, g, fail

    prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy,
                               airfoils=airfoils)

    try:
        dT, P = bemt.bemt_axial(prop, pitch, omega, tip_loss=tip_loss, mach_corr=mach_corr)
    except FloatingPointError:
        fail = 1
        return f, g, fail

    f = P
    print f
    print "Thrust = %s" % str(sum(dT))

    # Evaluate thrust constraint. Target thrust must be less than computed thrust
    g[0] = target_thrust - sum(dT)

    return f, g, fail


def objfun_const_chord(xn, **kwargs):

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
    omega = kwargs['omega']
    chord = kwargs['chord']*radius

    # The algorithm works with radius-normalized chord lengths, but the BEMT takes real chord lengths, so multiply by R
    twist0 = xn[0]
    dtwist = xn[1:]

    twist = calc_twist_dist(twist0, dtwist)

    f = 1000.
    fail = 0
    g = [1.0]

    prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy,
                               airfoils=airfoils)

    try:
        dT, P = bemt.bemt_axial(prop, pitch, omega, tip_loss=tip_loss, mach_corr=mach_corr)
    except FloatingPointError:
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


def main():
    ###########################################
    # Define some values
    ###########################################
    n_blades = 2
    n_elements = 5
    radius = unit_conversion.in2m(9.0)/2
    root_cutout = 0.1 * radius
    dy = float(radius-root_cutout)/n_elements
    dr = float(1)/n_elements
    y = root_cutout + dy*np.arange(1, n_elements+1)
    r = y/radius
    pitch = 0.0
    #airfoils = (('SDA1075_494p', 0.0, 1.0),)
    airfoils = (('simple', 0.0, 1.0),)
    thrust = 6.14
    max_chord = 100

    ###########################################
    # Set design variable bounds
    ###########################################
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
    omega_upper = 7000 *2*np.pi/60
    omega_lower = 4500 *2*np.pi/60

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

    opt_prob_ft = Optimization('Rotor in Hover w/ Fixed Twist', objfun_const_twist)
    opt_prob_ft.addVar('omega', 'c', value=omega, lower=omega_lower, upper=omega_upper)
    opt_prob_ft.addVar('chord0', 'c', value=chord0, lower=chord0_lower, upper=chord0_upper)
    opt_prob_ft.addVarGroup('dchord', n_elements-1, 'c', value=dchord, lower=dchord_lower, upper=dchord_upper)
    opt_prob_ft.addObj('f')
    opt_prob_ft.addCon('thrust', 'i')
    opt_prob_ft.addConGroup('c_lower', n_elements, 'i')
    opt_prob_ft.addConGroup('c_upper', n_elements, 'i')
    print opt_prob_ft

    # opt_prob_fc = Optimization('Rotor in Hover w/ Fixed Chord', objfun_const_chord)
    # opt_prob_fc.addVar('twist0', 'c', value=twist0, lower=twist0_lower, upper=twist0_upper)
    # opt_prob_fc.addVarGroup('dtwist', n_elements-1, 'c', value=dtwist, lower=dtwist_lower, upper=dtwist_upper)
    # opt_prob_fc.addObj('f')
    # opt_prob_fc.addCon('thrust', 'i')
    # print opt_prob_fc

    # Routine for optimizing chord with constant twist
    slsqp1 = SLSQP()
    slsqp1.setOption('IPRINT', 1)
    slsqp1.setOption('MAXIT', 200)
    slsqp1.setOption('ACC', 1e-7)
    [fstr, xstr, inform] = slsqp1(opt_prob_ft, sens_type='FD', n_blades=n_blades, n_elements=n_elements,
                                  root_cutout=root_cutout, radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch,
                                  airfoils=airfoils, thrust=thrust, max_chord=max_chord, tip_loss=False,
                                  mach_corr=False, omega=omega, twist=twist)
    print opt_prob_ft.solution(0)
    print "fstr = " + str(fstr)
    print "xstr = " + str(xstr)
    print "inform = " + str(inform)

    ## Routine for optimizing twist with a constant chord
    # slsqp2 = SLSQP()
    # slsqp2.setOption('IPRINT', 1)
    # slsqp2.setOption('MAXIT', 200)
    # slsqp2.setOption('ACC', 1e-7)
    # slsqp2(opt_prob_fc, sens_type='FD', n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout,
    #        radius=radius, dy=dy, dr=dr, y=y, r=r, pitch=pitch, airfoils=airfoils, thrust=thrust, tip_loss=False,
    #        mach_corr=False, omega=omega, chord=chord)
    # print opt_prob_fc.solution(0)

if __name__ == "__main__":
    main()