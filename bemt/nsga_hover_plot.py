import numpy as np
import matplotlib.pyplot as plt
import propeller
import bemt
import unit_conversion

radius = unit_conversion.in2m(9.0)/2
n_blades = 2
n_elements = 18
root_cutout = 0.1 * radius
dy = float(radius-root_cutout)/n_elements
dr = float(1)/n_elements
y = root_cutout + dy*np.arange(1, n_elements+1)
r = y/radius
Clalpha = 2 * np.pi
pitch = 0.0
airfoils = (('SDA1075', 0.0, 1.0),)
thrust = 5.22
max_chord = 0.4

############## Test 1 ######################
# omega = 680.495097 * 60 / 2 / np.pi
# theta0 = 0.514064 * 360 / 2 / np.pi
# chord0 = 0.383916
# dtwist = np.array([-0.171182, -0.174171, -0.159025, -0.007174, -0.000560, 0.058346, 0.056624, -0.046871, -0.062376,
#                    0.005443, 0.174406, -0.059963, -0.111853, -0.014548, 0.022623, -0.023352, -0.025071]) * 360/2/np.pi
# dchord = np.array([-0.006112, 0.048116, 0.045012, 0.046198, 0.003806, 0.049170, 0.016146, 0.019119, -0.041333, 0.005851,
#                    -0.006038, 0.019241, 0.012476, 0.029806, 0.025614, -0.019556, 0.000240])
############################################

############## Test 2 ######################
# omega = 676.204313 * 60 / 2 / np.pi
# theta0 = 0.287725 * 360 / 2 / np.pi
# chord0 = 0.119485
# dtwist = np.array([-0.129910, -0.003080, 0.028230, 0.150022, -0.059270, -0.107367, -0.128715, 0.102621, 0.026984,
#                    -0.136511, 0.096974, -0.010185, -0.039903, -0.040320, -0.012895, 0.009487, -0.068005]) * 360/2/np.pi
# dchord = np.array([0.049188, 0.036173, 0.034341, 0.039901, 0.049247, 0.049173, -0.022731, 0.048395, 0.021738, 0.031147,
#                    0.045673, 0.032719, 0.044858, 0.038226, 0.006546, 0.035800, -0.013972])
############################################

############ Test 3 ########################
omega = 680.540767
theta0 = 0.320063
chord0 = 0.122870
dtwist = np.array([-0.131272, -0.048603, -0.091425, 0.171687, -0.047733, -0.088764, -0.003366, 0.150650, 0.060453,
                   -0.168429, 0.058931, 0.039991, -0.095891, -0.017007, 0.005211, -0.066564, -0.058793])
dchord = np.array([0.039651, 0.040291, -0.033465, 0.033191, 0.019851, 0.030646, -0.046596, 0.049231, 0.046226,
                   -0.039365, 0.039439, 0.024125, 0.030256, 0.025809, -0.016043, -0.012905, -0.013165])


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


twist = calc_twist_dist(theta0, dtwist)
chord = calc_chord_dist(chord0, dchord)

chord_meters = chord * radius

# Run BEMT on propeller
prop = propeller.Propeller(twist, chord_meters, radius, n_blades, r, y, dr, dy, Clalpha, airfoils=airfoils)

CT, CP, CQ, Cl, dT, pCT, pCP, Re, aoa = bemt.bemt_axial(prop, pitch, omega)

print "C_P = " + str(pCP)
print "Thrust = " + str(sum(dT))
print "Re = " + str(Re)

plt.figure(1)
plt.plot(r, chord, '-')
plt.xlabel("radial station, r")
plt.ylabel("chord length, c/R")

plt.figure(2)
plt.plot(r, twist * 360/2/np.pi, '-')
plt.xlabel("radial station, r")
plt.ylabel("twist, degrees")

plt.show()



