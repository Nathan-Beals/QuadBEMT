import numpy as np
#import matplotlib.pyplot as plt
import propeller
import bemt
import unit_conversion
import lookup_table
import aero_coeffs


def main():
    radius = unit_conversion.in2m(9.0)/2
    n_blades = 2
    n_elements = 5
    root_cutout = 0.1 * radius
    dy = float(radius-root_cutout)/n_elements
    dr = float(1)/n_elements
    y = root_cutout + dy*np.arange(1, n_elements+1)
    r = y/radius
    pitch = 0.0
    allowable_Re = [100000., 90000., 80000., 70000., 60000.]
    airfoils = (('SDA1075_494p', 0.0, 1.0),)
    #airfoils = (('simple', 0.0, 1.0),)

    Cl_tables = {}
    Cd_tables = {}
    # Get lookup tables
    if any(airfoil[0] != 'simple' for airfoil in airfoils):
        for airfoil in airfoils:
            Cl_table, Cd_table = aero_coeffs.create_Cl_Cd_table(airfoil[0])

            Cl_tables[airfoil[0]] = Cl_table
            Cd_tables[airfoil[0]] = Cd_table

    # Create list of Cl functions. One for each Reynolds number
    if Cl_tables:
        Cl_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cl_fun(Re, Cl_tables[airfoils[0][0]]) for Re in allowable_Re]))
        Cd_funs = dict(zip(allowable_Re, [aero_coeffs.get_Cd_fun(Re, Cd_tables[airfoils[0][0]]) for Re in allowable_Re]))

    ########### DA4002 ########################
    # omega = 5943 * 2*np.pi/60
    # twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
    #                   18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
    # chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
    #                   0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
    # dtwist_base = np.array([twist[i+1]-twist[i] for i in xrange(len(twist)-1)])
    # dchord_base = np.array([chord[i+1]-chord[i] for i in xrange(len(chord)-1)])
    # twist0_base = twist[0]
    # chord0_base = chord[0]
    ###########################################

    ######### Low Element DA4002 ##############
    omega = 5943.0 * 2*np.pi/60
    chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
                      0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
    chord = np.array([chord[i] for i in [0, 4, 8, 12, 16]])
    twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
                      18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
    twist = np.array([twist[i] for i in [0, 4, 8, 12, 16]])

    chord0 = chord[0]
    twist0 = twist[0]
    dchord = np.array([chord[i+1]-chord[i] for i in xrange(len(chord)-1)])
    dtwist = np.array([twist[i+1]-twist[i] for i in xrange(len(twist)-1)])

    chord0_base = np.array(chord0)
    twist0_base = np.array(twist0)
    dchord_base = np.array(dchord)
    dtwist_base = np.array(dtwist)

    ####### Low Element SLSQP Fixed Chord Correct Solidity From DA4002 #####
    # twist0 = 0.501503
    # dtwist = np.array([-0.058172, -0.055531, -0.102258, -0.061148])
    ############################################################

    ###### Low Element SLSQP Fixed Twist Correct Solidity #######
    # chord0 = 0.123033
    # dchord = np.array([-0.050000, 0.050000, 0.050000, 0.050000])

    ########### Hyperbolic 5 element chord ######################
    chord0 = 0.604286
    dchord = np.array([-0.236460, -0.103451, -0.058034, -0.037141])
    #############################################################

    ######## From Hyper 5 with optimized omega ###################
    # omega_opt = 616.666829
    # chord0 = 0.160608
    # dchord = np.array([-0.089870, -0.056311, 0.225470, -0.030641])
    ##############################################################

    ######## From Hyper 5 w omega and simple airfoil #############
    # airfoils = (('simple', 0.0, 1.0),)
    # omega_opt = 471.238898
    # chord0 = 1.027156
    # dchord = np.array([-0.744996, -0.023607, -0.032289, -0.008823])
    ##############################################################

    ####### NSGA Simple Hover 5000pop 1000gen ####################
    #airfoils = (('simple', 0.0, 1.0),)
    #allowable_Re = np.array([100000.])
    # omega = 471.413361
    # twist0 = 0.899588
    # chord0 = 0.467506
    # dtwist = np.array([-0.127577, -0.171353, -0.160939, -0.052663, 0.026401, -0.058003, -0.011649, 0.009380, -0.003034])
    # dchord = np.array([-0.031194, -0.088289, 0.070828, -0.072690, -0.099504, -0.017530, -0.055469, -0.039089, -0.011930])

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

    twist = calc_twist_dist(twist0, dtwist)
    chord = calc_chord_dist(chord0, dchord)
    chord_meters = chord * radius

    # twist_base = calc_twist_dist(twist0_base, dtwist_base)
    # chord_base = calc_chord_dist(chord0_base, dchord_base)
    # chord_base_meters = chord_base * radius
    ##Base prop
    #base_prop = propeller.Propeller(twist_base, chord_base_meters, radius, n_blades, r, y, dr, dy, airfoils=(('SDA1075_494p', 0.0, 1.0),))

    # Optimized prop
    prop = propeller.Propeller(twist, chord_meters, radius, n_blades, r, y, dr, dy, airfoils=airfoils,
                               Cl_tables=Cl_tables, Cd_tables=Cd_tables)

    dT, dP, P, Cd, Cl, ures, chord_meters, dL, inflow, inflow_ang, eff_aoa, dFx, dFz, Re, Re_act, Re_app, Re_app_act = \
        bemt.bemt_axial(prop, pitch, omega, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs, tip_loss=False,
                        mach_corr=False, output='long')

    # base_vals = bemt.bemt_axial(base_prop, pitch, omega, tip_loss=False, mach_corr=False, output='long')
    # base_vals_newbemt = bemt.bemt_axial_alt(base_prop, pitch, omega, tip_loss=False, mach_corr=False)

    # opt_vals_newbemt = bemt.bemt_axial_alt(prop, pitch, omega, tip_loss=False, mach_corr=False)

    # print "Thrust_base = " + str(sum(base_vals[0]))
    # print "Power_base = " + str(sum(base_vals[1]))
    # print "Thrust_opt = " + str(sum(dT))
    # print "Power_opt = " + str(sum(dP))
    #
    # plt.figure(5)
    # plt.plot(r, chord, '-b')
    # plt.plot(r, chord_base, '-r')
    # plt.xlabel("radial station, r")
    # plt.ylabel("chord")
    # plt.legend(['optimized', 'base'])
    #
    # plt.figure(4)
    # plt.plot(r, twist * 360/2/np.pi, '-b')
    # plt.plot(r, twist_base * 360/2/np.pi, '-r')
    # plt.xlabel("radial station, r")
    # plt.ylabel("twist, degrees")
    # plt.legend(['optimized', 'base'], loc=2)
    #
    # plt.figure(1)
    # plt.plot(r, base_vals_newbemt[2], '-b')
    # plt.plot(r, opt_vals_newbemt[2], '-r')
    # plt.xlabel('radial location, r')
    # plt.ylabel('induced power')
    # plt.legend(['base', 'opt'])
    #
    # plt.figure(2)
    # plt.plot(r, base_vals_newbemt[3], '-b')
    # plt.plot(r, opt_vals_newbemt[3], '-r')
    # plt.xlabel('radial location, r')
    # plt.ylabel('profile power')
    # plt.legend(['base', 'opt'])
    #
    # plt.figure(3)
    # plt.plot(r, base_vals_newbemt[1], '-b')
    # plt.plot(r, opt_vals_newbemt[1], '-r')
    # plt.xlabel('radial location, r')
    # plt.ylabel('power')
    # plt.legend(['base', 'opt'])
    #
    # plt.show()

    print "omega = " + str(omega)
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
    print "Re_act = " + str(Re_act)
    print "Re_app = " + str(Re_app)
    print "Re_app_act = " + str(Re_app_act)

    # plt.figure(1)
    # plt.plot(r, chord, '-b')
    # #plt.plot(r, chord_base, '-r')
    # plt.xlabel("radial station, r")
    # plt.ylabel("chord")
    # #plt.legend(['optimized', 'base'])
    #
    # plt.figure(2)
    # plt.plot(r, twist * 360/2/np.pi, '-b')
    # #plt.plot(r, twist_base * 360/2/np.pi, '-r')
    # plt.xlabel("radial station, r")
    # plt.ylabel("twist, degrees")
    # #plt.legend(['optimized', 'base'], loc=2)

    #
    # plt.figure(3)
    # plt.plot(r, inflow, '-b')
    # plt.plot(r, base_vals[8], '-r')
    # plt.xlabel("radial station, r")
    # plt.ylabel("inflow")
    # plt.legend(['optimized', 'base'], loc=2)
    #
    # plt.figure(4)
    # plt.plot(r, Cl/Cd, '-b')
    # plt.plot(r, base_vals[4]/base_vals[3], '-r')
    # plt.xlabel("radial station, r")
    # plt.ylabel("Cl/Cd")
    # plt.legend(['optimized', 'base'], loc=2)
    #
    # plt.figure(5)
    # plt.plot(r, Re, '-b')
    # plt.plot(r, base_vals[13], '-r')
    # plt.xlabel("radial station, r")
    # plt.ylabel("Re")
    # plt.legend(['optimized', 'base'], loc=2)
    #
    # plt.figure(6)
    # plt.plot(r, ures, '-b')
    # plt.plot(r, base_vals[5], '-r')
    # plt.xlabel("radial station, r")
    # plt.ylabel("ures")
    # plt.legend(['optimized', 'base'], loc=2)
    #
    # plt.figure(7)
    # plt.plot(r, dP, '-b')
    # plt.plot(r, base_vals[1], '-r')
    # plt.xlabel("radial station, r")
    # plt.ylabel("Power")
    # plt.legend(['optimized', 'base'], loc=2)
    #
    #plt.show()

if __name__ == '__main__':
    main()



