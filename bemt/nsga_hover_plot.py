import numpy as np
import matplotlib.pyplot as plt
import propeller
import bemt
import unit_conversion


def main():
    radius = unit_conversion.in2m(9.0)/2
    n_blades = 2
    n_elements = 18
    root_cutout = 0.1 * radius
    dy = float(radius-root_cutout)/n_elements
    dr = float(1)/n_elements
    y = root_cutout + dy*np.arange(1, n_elements+1)
    r = y/radius
    pitch = 0.0
    airfoils = (('SDA1075_494p', 0.0, 1.0),)
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
    # omega = 680.540767
    # theta0 = 0.320063
    # chord0 = 0.122870
    # dtwist = np.array([-0.131272, -0.048603, -0.091425, 0.171687, -0.047733, -0.088764, -0.003366, 0.150650, 0.060453,
    #                    -0.168429, 0.058931, 0.039991, -0.095891, -0.017007, 0.005211, -0.066564, -0.058793])
    # dchord = np.array([0.039651, 0.040291, -0.033465, 0.033191, 0.019851, 0.030646, -0.046596, 0.049231, 0.046226,
    #                    -0.039365, 0.039439, 0.024125, 0.030256, 0.025809, -0.016043, -0.012905, -0.013165])
    ############################################

    ############ Test 4 ########################
    # omega = 614.410800
    # theta0 = 0.370607
    # chord0 = 0.126276
    # dtwist = np.array([-0.158407, 0.080290, 0.088475, 0.058860, 0.136430, -0.099698, -0.095792, -0.132567, -0.144394,
    #                    0.155031, -0.112972, 0.048031, -0.049668, 0.001085, 0.012840, -0.129098, -0.049580])
    # dchord = np.array([0.012987, 0.049636, 0.042496, 0.017863, 0.040340, -0.000895, 0.016335, 0.049916, 0.020678, -0.011996,
    #                    -0.023302, -0.018517, -0.007816, 0.013603, -0.030474, 0.022230, 0.045440])
    ###########################################

    ############ Test 5 #######################
    # omega = 624.698621
    # theta0 = 0.331718
    # chord0 = 0.120594
    # dtwist = np.array([-0.037092, 0.037449, 0.087490, 0.060749, -0.033991, -0.147636, 0.010353, -0.077433, -0.156896,
    #                    0.155319, -0.071940, 0.012944, -0.043551, -0.001275, 0.027710, -0.133497, -0.045003])
    # dchord = np.array([0.035663, 0.048244, 0.026311, 0.025575, 0.032583, 0.019873, 0.011969, 0.028818, 0.031445, 0.001121,
    #                    -0.018293, -0.028388, 0.013017, 0.016265, -0.030051, 0.004988, 0.044851])
    ###########################################

    ########### No Tip loss no GeoCons ########
    # omega = 537.809942
    # theta0 = 0.538470
    # chord0 = 0.875314
    # dtwist = np.array([-0.118980, -0.021633, 0.024974, -0.174503, -0.072010, -0.173192, -0.058247, -0.102644, 0.028660,
    #                    0.161007, 0.069044, -0.062104, -0.037059, -0.081645, 0.024131, 0.108387, -0.010453])
    # dchord = np.array([0.038025, -0.015922, -0.022282, -0.002897, -0.034241, 0.030693, -0.050589, -0.007714, 0.032086,
    #                    -0.018731, 0.045765, -0.054933, 0.076813, 0.048910, -0.079941, -0.072764, 0.073052])
    ###########################################

    ########## No Tip/Mach no NegInflows ######
    # omega = 525.624058
    # theta0 = 0.355572
    # chord0 = 0.716411
    # dtwist = np.array([0.008419, 0.061331, -0.124701, -0.139102, -0.057765, 0.089032, 0.087092, -0.129683, 0.019089,
    #                    0.039940, -0.102672, -0.072554, 0.096547, -0.068370, -0.056449, -0.000731, 0.098145])
    # dchord = np.array([0.037889, -0.076113, 0.032740, -0.087678, 0.065222, 0.039185, 0.081976, -0.058099, -0.050094,
    #                    -0.046166, 0.000603, 0.036750, 0.095209, 0.068672, 0.071175, -0.093865, -0.076424])
    ###########################################


    ######### No Tip/Mach new SDA1075 Fast ####
    # omega = 550.413232
    # theta0 = 0.595659
    # chord0 = 0.347843
    # dtwist = np.array([0.083052, -0.139437, -0.131492, -0.033905, -0.014306, 0.046297, -0.091609, 0.060680, -0.055783,
    #                    -0.049641, 0.078339, -0.167067, -0.115555, 0.083237, 0.036530, 0.064380, -0.024836])
    # dchord = np.array([-0.063772, -0.006073, -0.099464, 0.026802, 0.066349, 0.017330, -0.026736, 0.058189, -0.063305,
    #                    0.070179, -0.096371, -0.052586, 0.075146, 0.027074, 0.028175, -0.067737, -0.033894])
    ###########################################

    ######## No Tip/Mach 494p, 372pop, 250 gen #
    # omega = 527.857015
    # theta0 = 0.464779
    # chord0 = 0.385671
    # dtwist = np.array([0.046148, -0.103487, -0.079213, -0.034663, -0.028742, 0.119307, -0.095230, 0.032763, -0.052037,
    #                    -0.039653, 0.044043, -0.086893, -0.038794, 0.039525, 0.026184, 0.040540, -0.039594])
    # dchord = np.array([-0.057880, 0.023178, -0.090472, 0.018907, 0.043419, 0.017060, -0.008579, 0.046066, -0.035405,
    #                    0.027569, -0.097129, -0.048276, 0.054540, 0.025390, 0.006726, -0.077676, 0.007061])
    ############################################

    ########### NoTipMach 600pop 40 gen ########
    # omega = 578.926935
    # theta0 = 0.006811
    # chord0 = 0.303445
    # dtwist = np.array([0.128549, 0.110798, -0.003510, 0.094296, 0.126790, -0.118536, -0.141363, 0.104940, 0.035874,
    #                    -0.075178, -0.082803, 0.016918, 0.022771, -0.045963, -0.086065, 0.055630, -0.094279])
    # dchord = np.array([0.046937, -0.029662, -0.059134, 0.012516, 0.012557, 0.086067, -0.046761, -0.048293, 0.002057,
    #                    0.060100, 0.017885, -0.091718, 0.010184, 0.045959, -0.032151, 0.096939, -0.062810])
    ############################################

    ########### 0.8 crossover 40 gen #############
    # omega = 527.534826
    # theta0 = 0.386225
    # chord0 = 0.376514
    # dtwist = np.array([-0.096199, -0.158534, 0.053815, 0.172098, -0.125135, 0.096198, -0.117144, 0.084459, 0.092395,
    #                    -0.098802, -0.032193, -0.020467, 0.018538, -0.040243, 0.021903, -0.051685, 0.078419])
    # dchord = np.array([-0.010890, -0.017131, -0.014639, 0.004333, -0.018299, -0.018490, -0.007401, -0.001974, -0.011334,
    #                    -0.004361, -0.010958, -0.012996, 0.001096, 0.008517, -0.018022, -0.015253, -0.014777])
    ############################################

    ########## 0.8 crossover 150 gen ############
    omega = 526.938693
    theta0 = 0.434501
    chord0 = 0.376514
    dtwist = np.array([-0.101670, -0.153288, 0.063368, 0.173516, -0.162831, 0.086124, -0.122493, 0.103849, 0.053405,
                       -0.134726, -0.000400, 0.022519, 0.006426, -0.029408, 0.016352, -0.082736, 0.061399])
    dchord = np.array([-0.010467, -0.015021, -0.019865, 0.009399, -0.017646, -0.018732, -0.014597, -0.002709, -0.010477,
                       -0.002962, -0.005953, -0.015082, -0.004879, 0.010107, -0.015014, -0.015507, -0.013003])


    ########### DA4002 ########################
    # omega = 5943 * 2*np.pi/60
    # twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
    #                   18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
    # chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
    #                   0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
    # dtwist = np.array([twist[i+1]-twist[i] for i in xrange(len(twist)-1)])
    # dchord = np.array([chord[i+1]-chord[i] for i in xrange(len(chord)-1)])
    # theta0 = twist[0]
    # chord0 = chord[0]
    ###########################################

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
    prop = propeller.Propeller(twist, chord_meters, radius, n_blades, r, y, dr, dy, airfoils=airfoils)

    dT, dP, P, Cd, Cl, ures, chord_meters, dL, inflow, inflow_ang, eff_aoa, dFx, dFz, Re = bemt.bemt_axial(prop, pitch, omega, tip_loss=False, mach_corr=False)
    #dT, P = bemt.bemt_axial(prop, pitch, omega, tip_loss=False, mach_corr=False)

    print "omega = " + str(omega)
    print "P = " + str(P)
    print "dP = " + str(dP)
    print "chord = " + str(chord_meters)
    print "Cd = " + str(Cd)
    print "Cl = " + str(Cl)
    print "Cl/Cd = " + str(Cl/Cd)
    print "dL = " + str(dL)
    print "inflow = " + str(inflow)
    print "phi = " + str(inflow_ang)
    print "alpha = " + str(eff_aoa*360/2/np.pi)
    print "ures = " + str(ures)
    print "twist = " + str(twist)
    print "dFx = " + str(dFx)
    print "dFz = " + str(dFz)
    print "Thrust = " + str(sum(dT))
    print "Power = " + str(sum(dP))
    print "Re = " + str(Re)
    #
    plt.figure(1)
    plt.plot(r, chord, '-')
    plt.xlabel("radial station, r")
    plt.ylabel("chord")
    #
    plt.figure(2)
    plt.plot(r, twist * 360/2/np.pi, '-')
    plt.xlabel("radial station, r")
    plt.ylabel("twist, degrees")
    #
    plt.show()

if __name__ == '__main__':
    main()



