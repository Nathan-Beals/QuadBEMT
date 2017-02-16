import unit_conversion
import propeller
import bemt
from pyOpt import Optimization
from pyOpt import NSGA2
from pyOpt import SLSQP
from pyOpt import CONMIN
import numpy as np


def counted(f):
    def wrapped(*args, **kwargs):
        print wrapped.calls
        wrapped.calls += 1
        return f(*args, **kwargs)
    wrapped.calls = 0
    return wrapped


@counted
def objfun(xn, **kwargs):

    print 'objfun called'

    int_radius = kwargs['radius']
    int_r = kwargs['r']
    int_y = kwargs['y']
    int_dr = kwargs['dr']
    int_dy = kwargs['dy']
    int_n_blades = kwargs['n_blades']
    int_airfoils = kwargs['airfoils']
    int_pitch = kwargs['pitch']
    target_thrust = kwargs['thrust']
    cval_max = kwargs['max_chord']
    omega = xn[0]
    twist0 = xn[1]
    chord0 = xn[2]*radius
    dtwist = np.array(xn[3:len(r)+2])
    dchord = np.array([x*radius for x in xn[len(r)+2:2*len(r)+2]])

    int_twist = calc_twist_dist(twist0, dtwist)
    int_chord = calc_chord_dist(chord0, dchord)

    # int_twist = np.array(xn[:len(r)])
    # int_chord = np.array([x*radius for x in xn[len(r):]])

    f = 1000
    fail = 0
    g = [1.0] * (2*len(r)+1)

    # Calculate geometric constraint values and return immediately if there are any failures
    g[1:] = get_geocons(int_chord, cval_max, int_radius)
    if any(geocon > 0.0 for geocon in g[1:]):
        print "geocons violated"
        return f, g, fail

    prop = propeller.Propeller(int_twist, int_chord, int_radius, int_n_blades, int_r, int_y, int_dr, int_dy,
                               airfoils=int_airfoils)

    try:
        dT, P = bemt.bemt_axial(prop, int_pitch, omega, tip_loss=False, mach_corr=False)
    except FloatingPointError:
        fail = 1
        return f, g, fail

    f = P
    print omega
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


###########################################
# Define some values
###########################################
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
thrust = 4.94
max_chord = 0.4

###########################################
# Set design variable bounds
###########################################
chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
                  0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
                  18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
dtwist_start = np.array([twist[i+1]-twist[i] for i in xrange(len(twist)-1)])
dchord_start = np.array([chord[i+1]-chord[i] for i in xrange(len(chord)-1)])
theta0_start = twist[0]
chord0_start = chord[0]
# chord = np.array([chord[i] for i in [0, 4, 8, 12, 17]])
# twist = np.array([twist[i] for i in [0, 4, 8, 12, 17]])
# chord = np.array([chord[i] for i in [12]])
# twist = np.array([twist[i] for i in [12]])
omega_lower = 5000.0 * 2*np.pi/60
omega_upper = 6500.0 * 2*np.pi/60
omega_start = 6100.0 * 2*np.pi/60
# twist_lower = 0.0 * 2*np.pi/360
# twist_upper = 50.0 * 2*np.pi/360
# twist_start = twist
# chord_lower = 0
# chord_upper = 0.4
# chord_start = chord

# Twist at the hub must be less than or equal to arcsin(hub_height/hub_diameter), approx 23 degrees
#theta0_start = 20.0 * 2 * np.pi / 360
theta0_lower = 0.0 * 2 * np.pi / 360
theta0_upper = 60.0 * 2 * np.pi / 360
# Chord values at the hub must be less than or equal to the diameter of the hub
#chord0_start = 0.12
chord0_upper = 10
chord0_lower = 0.05

#dtwist_start = 0.0 * 2 * np.pi / 360
dtwist_lower = -10.0 * 2 * np.pi / 360
dtwist_upper = 10.0 * 2 * np.pi / 360
#dchord_start = 0.0
dchord_lower = -0.1
dchord_upper = 0.1
# chord_lower = unit_conversion.in2m(0.25)
# chord_upper = unit_conversion.in2m(3.0)
# chord_start = unit_conversion.in2m(0.5)

opt_prob = Optimization('Rotor in Hover', objfun)
opt_prob.addVar('omega', 'c', value=omega_start, lower=omega_lower, upper=omega_upper)
opt_prob.addVar('theta0', 'c', value=theta0_start, lower=theta0_lower, upper=theta0_upper)
opt_prob.addVar('chord0', 'c', value=chord0_start, lower=chord0_lower, upper=chord0_upper)
opt_prob.addVarGroup('dtwist', n_elements-1, 'c', value=dtwist_start, lower=dtwist_lower, upper=dtwist_upper)
opt_prob.addVarGroup('dchord', n_elements-1, 'c', value=dchord_start, lower=dchord_lower, upper=dchord_upper)
# opt_prob.addVarGroup('twist', n_elements, 'c', value=twist_start, lower=twist_lower, upper=twist_upper)
# opt_prob.addVarGroup('chord', n_elements, 'c', value=chord_start, lower=chord_lower, upper=chord_upper)
opt_prob.addObj('f')
opt_prob.addCon('thrust', 'i')
opt_prob.addConGroup('c_lower', n_elements, 'i')
opt_prob.addConGroup('c_upper', n_elements, 'i')
print opt_prob

nsga2 = NSGA2()
nsga2.setOption('PrintOut', 2)
nsga2.setOption('PopSize', 372)
nsga2.setOption('maxGen', 40)
nsga2(opt_prob, n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout, radius=radius, dy=dy,
      dr=dr, y=y, r=r, pitch=pitch, airfoils=airfoils, thrust=thrust, max_chord=max_chord)
print opt_prob.solution(0)

# conmin = CONMIN()
# conmin.setOption('IPRINT', 0)
# conmin(opt_prob, sens_type='CS', n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout, radius=radius,
#        dy=dy, dr=dr, y=y, r=r, Clalpha=Clalpha, pitch=pitch, airfoils=airfoils, thrust=thrust)
# print opt_prob.solution(0)


# slsqp = SLSQP()
# slsqp.setOption('IPRINT', 1)
# slsqp(opt_prob, sens_type='FD', n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout, radius=radius, dy=dy,
#       dr=dr, y=y, r=r, pitch=pitch, airfoils=airfoils, thrust=thrust, max_chord=max_chord)
# print opt_prob.solution(0)