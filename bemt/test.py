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
    int_Clalpha = kwargs['Clalpha']
    int_pitch = kwargs['pitch']
    target_thrust = kwargs['thrust']
    omega = xn[0]
    int_twist = np.array(xn[1:len(r)+1])
    int_chord = np.array([x*radius for x in xn[len(r)+1:2*len(r)+1]])
    # int_twist = np.array(xn[:len(r)])
    # int_chord = np.array([x*radius for x in xn[len(r):]])

    f = 1000
    fail = 0
    g = [1.0] * (4*(len(r)-1)+1)

    # Calculate geometric constraint values and return immediately if there are any failures
    a1 = 2.0
    a2 = unit_conversion.deg2rad(10.0)
    g[1:] = get_geocons(int_twist, int_chord, a1, a2)
    if any(g[1:]) > 0.0:
        print "geocons violated"
        return f, g, fail

    prop = propeller.Propeller(int_twist, int_chord, int_radius, int_n_blades, int_r, int_y, int_dr, int_dy,
                               int_Clalpha, airfoils=int_airfoils)

    try:
        CT, CP, CQ, Cl, dT, pCT, pCP, Re, aoa = bemt.bemt_axial(prop, int_pitch, omega)
    except FloatingPointError:
        fail = 1
        return f, g, fail

    f = pCP
    print omega
    print "CT = %s" % str(pCT)
    print f
    print "Thrust = %s" % str(sum(dT))

    # Evaluate thrust constraint. Target thrust must be less than computed thrust
    g[0] = target_thrust - sum(dT)

    return f, g, fail


def get_geocons(t, c, a1, a2):
    """
    This function evaluates geometric constraints on the alternative and returns the list.
    :param t: twist values
    :param c: chord values
    :param a1: twist stretch factor
    :param a2: chord stretch factor
    :return: constraint values
    """
    geocons = [0.0] * (4*(len(t)-1))
    con_indx = 0
    for i in xrange(len(t)-1):
        geocons[con_indx] = t[i+1] - a2*t[i]
        con_indx += 1
    for i in xrange(len(t)-1):
        geocons[con_indx] = t[i]/a2 - t[i+1]
        con_indx += 1
    for i in xrange(len(c)-1):
        try:
            geocons[con_indx] = c[i+1] - a1*c[i]
        except IndexError:
            print con_indx
            print i
            raise
        con_indx += 1
    for i in xrange(len(c)-1):
        geocons[con_indx] = c[i]/a1 - i + 1
        con_indx += 1
    return geocons


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
Clalpha = 2 * np.pi
pitch = 0.0
airfoils = (('SDA1075', 0.0, 1.0),)
thrust = 5.22

###########################################
# Set design variable bounds
###########################################
chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
                  0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
                  18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
# chord = np.array([chord[i] for i in [0, 4, 8, 12, 17]])
# twist = np.array([twist[i] for i in [0, 4, 8, 12, 17]])
# chord = np.array([chord[i] for i in [12]])
# twist = np.array([twist[i] for i in [12]])
omega_lower = 5000.0 * 2*np.pi/60
omega_upper = 6500.0 * 2*np.pi/60
omega_start = 5900.0 * 2*np.pi/60
twist_lower = 0.0 * 2*np.pi/360
twist_upper = 50.0 * 2*np.pi/360
twist_start = twist
chord_lower = 0
chord_upper = 0.4
chord_start = chord
# chord_lower = unit_conversion.in2m(0.25)
# chord_upper = unit_conversion.in2m(3.0)
# chord_start = unit_conversion.in2m(0.5)

opt_prob = Optimization('Rotor in Hover', objfun)
opt_prob.addVar('omega', 'c', value=omega_start, lower=omega_lower, upper=omega_upper)
opt_prob.addVarGroup('twist', n_elements, 'c', value=twist_start, lower=twist_lower, upper=twist_upper)
opt_prob.addVarGroup('chord', n_elements, 'c', value=chord_start, lower=chord_lower, upper=chord_upper)
opt_prob.addObj('f')
opt_prob.addCon('thrust', 'i')
opt_prob.addConGroup('tdec', n_elements-1, 'i')
opt_prob.addConGroup('t2nd', n_elements-1, 'i')
opt_prob.addConGroup('cdec', n_elements-1, 'i')
opt_prob.addConGroup('c2nd', n_elements-1, 'i')
print opt_prob

nsga2 = NSGA2()
nsga2.setOption('PrintOut', 0)
nsga2(opt_prob, n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout, radius=radius, dy=dy,
      dr=dr, y=y, r=r, Clalpha=Clalpha, pitch=pitch, airfoils=airfoils, thrust=thrust, PrintOut=2)
print opt_prob.solution(0)