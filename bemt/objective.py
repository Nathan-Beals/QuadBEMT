import unit_conversion
import propeller
import bemt
from pyOpt import Optimization
from pyOpt import NSGA2
from pyOpt import SLSQP
from pyOpt import CONMIN
import numpy as np


def objfun(xn, **kwargs):

    print 'objfun called'
    print xn
    print len(xn)

    radius = kwargs['radius']
    r = kwargs['r']
    y = kwargs['y']
    dr = kwargs['dr']
    dy = kwargs['dy']
    n_blades = kwargs['n_blades']
    airfoils = kwargs['airfoils']
    Clalpha = kwargs['Clalpha']
    pitch = kwargs['pitch']
    target_thrust = kwargs['thrust']

    omega = xn[0]
    twist = xn[1:len(r)+1]
    chord = xn[len(r)+1:2*len(r)+1]*radius

    g = [0.0]*(4*(len(r)-1)+1)
    f = 0.0
    fail = 0

    prop = propeller.Propeller(twist, chord, radius, n_blades, r, y, dr, dy, Clalpha, airfoils=airfoils)

    try:
        CT, CP, CQ, Cl, dT, pCT, pCP, Re, aoa = bemt.bemt_axial(prop, pitch, omega)
    except FloatingPointError:
        fail = 1
        return f, g, fail

    f = CP
    print "CT = %s" % str(CT)
    print f
    print "Thrust = %s" % str(sum(dT))

    # Evaluate constraints
    a1 = 2.0
    a2 = unit_conversion.deg2rad(10.0)
    g[0] = abs(sum(dT) - target_thrust) - 0.001
    con_indx = 1
    for i in xrange(len(twist)-1):
        g[con_indx] = twist[i+1] - a2*twist[i]
        con_indx += 1
    for i in xrange(len(twist)-1):
        g[con_indx] = twist[i]/a2 - twist[i+1]
        con_indx += 1
    for i in xrange(len(chord)-1):
        try:
            g[con_indx] = chord[i+1] - a1*chord[i]
        except IndexError:
            print con_indx
            print i
            raise
        con_indx += 1
    for i in xrange(len(chord)-1):
        g[con_indx] = chord[i]/a1 - i + 1
        con_indx += 1

    #print "ObjFun exit"
    return f, g, fail

###########################################
# Define some values
###########################################
n_blades = 2
n_elements = 5
radius = unit_conversion.in2m(9.0)/2
root_cutout = 0.0 * radius
dy = float(radius-root_cutout)/n_elements
dr = float(1)/n_elements
y = root_cutout + dy*np.arange(1, n_elements+1)
r = y/radius
Clalpha = 2 * np.pi
pitch = 0.0
airfoils = (('SDA1075', 0, 1),)
thrust = 3.435

###########################################
# Set design variable bounds
###########################################
chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
                  0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
                  18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
chord = np.array([chord[i] for i in [0, 4, 8, 12, 17]])
twist = np.array([twist[i] for i in [0, 4, 8, 12, 17]])
omega_lower = 0.0
omega_upper = 8000.0 * 2*np.pi/60
omega_start = 5500.0 * 2*np.pi/60
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

# nsga2 = NSGA2()
# nsga2.setOption('PrintOut', 0)
# nsga2(opt_prob, n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout, radius=radius, dy=dy,
#       dr=dr, y=y, r=r, Clalpha=Clalpha, pitch=pitch, airfoils=airfoils, thrust=thrust)
# print opt_prob.solution(0)

# conmin = CONMIN()
# conmin.setOption('IPRINT', 0)
# conmin(opt_prob, sens_type='CS', n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout, radius=radius,
#        dy=dy, dr=dr, y=y, r=r, Clalpha=Clalpha, pitch=pitch, airfoils=airfoils, thrust=thrust)
# print opt_prob.solution(0)


# slsqp = SLSQP()
# slsqp.setOption('IPRINT', -1)
# slsqp(opt_prob, sens_type='FD', n_blades=n_blades, n_elements=n_elements, root_cutout=root_cutout, radius=radius, dy=dy,
#       dr=dr, y=y, r=r, Clalpha=Clalpha, pitch=pitch, airfoils=airfoils, thrust=thrust)
# print opt_prob.solution(0)

