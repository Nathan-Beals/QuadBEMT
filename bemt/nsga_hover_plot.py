import numpy as np
import matplotlib.pyplot as plt
import propeller
import bemt
import unit_conversion
import lookup_table
import aero_coeffs


def main():

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

    radius = unit_conversion.in2m(9.0)/2
    n_blades = 2
    n_elements = 10
    root_cutout = 0.1 * radius
    dy = float(radius-root_cutout)/n_elements
    dr = float(1)/n_elements
    y = root_cutout + dy*np.arange(1, n_elements+1)
    r = y/radius
    pitch = 0.0
    allowable_Re = [1000000., 500000., 250000., 100000., 90000., 80000., 70000., 60000., 50000., 40000., 30000., 20000., 10000.]
    airfoils = (('SDA1075_494p', 0.0, 1.0),)
    tip_loss = True

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

    ######### DA4002 ##############
    omega = 5943.0 * 2*np.pi/60
    chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
                      0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
    # chord = np.array([chord[i] for i in [0, 4, 8, 12, 16]])
    # chord = np.array([chord[i] for i in [0, 3, 6, 9, 11, 13, 15, 17]])
    chord = np.array([chord[i] for i in [0, 2, 4, 6, 8, 10, 12, 14, 15, 17]])
    twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
                      18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
    # twist = np.array([twist[i] for i in [0, 4, 8, 12, 16]])
    # twist = np.array([twist[i] for i in [0, 3, 6, 9, 11, 13, 15, 17]])
    twist = np.array([twist[i] for i in [0, 2, 4, 6, 8, 10, 12, 14, 15, 17]])
    base = (omega, chord, twist, 'DA4002')

    #################### (OPT # 1) NSGA 10 pt 5000p 500g Reallowed=All T=4.61 xinit=base tiploss ############
    chord = np.array([0.11470868, 0.21117422, 0.30252118, 0.39582241, 0.42498534, 0.41241736, 0.32928823, 0.23259937,
                      0.17993009, 0.11100332])
    twist = np.array([0.55653743, 0.45683294, 0.38657942, 0.39545834, 0.3270524, 0.3053103, 0.27289751, 0.23934658,
                      0.20193655, 0.034115])
    omega = 5286.48230542 * 2*np.pi/60
    opt1 = (omega, chord, twist, 'NSGA10e_5000p_500g_Reallowed=all_T=4.61_tiploss')

    ############### (OPT # 2) SLSQP Opt of DA4002 10 element. Reallow=All, T=4.61, Tiploss ##################
    chord = np.array([0.11994618, 0.19944238, 0.2781612, 0.33993244, 0.33206394, 0.28471642, 0.22285041, 0.18446874,
                      0.12841128, 0.01537926])
    twist = np.array([0.51724241, 0.4268415, 0.39644603, 0.33147393, 0.31163335, 0.28957352, 0.25814615, 0.24327558,
                      0.20237492, 0.00653631])
    omega = 5943.55347904 * 2*np.pi/60
    opt2 = (omega, chord, twist, 'SLSQP of DA4002_10e_Reallowed=all_T=4.61_tiploss')

    ############### (OPT # 3) SLSQP Opt of OPT # 1, Reallow=All, T=4.61, Tiploss #######################################
    chord = np.array([0.12000001, 0.22000002, 0.30478926, 0.4047892, 0.41469428, 0.35391398, 0.31830687, 0.21970192,
                      0.16508119, 0.11102533])
    twist = np.array([0.59057564, 0.4386199, 0.42691716, 0.35621971, 0.33793345, 0.30917999, 0.27913501, 0.25924289,
                      0.21518868, 0.04164781])
    omega = 5286.15517496 * 2*np.pi/60
    opt3 = (omega, chord, twist, 'SLSQP of '+opt1[3])

    cases2run = (opt1, opt3)

    def get_performance(o, c, t, tip_loss, mach_corr=False):
        chord_meters = c * radius
        prop = propeller.Propeller(t, chord_meters, radius, n_blades, r, y, dr, dy, airfoils=airfoils,
                                   Cl_tables=Cl_tables, Cd_tables=Cd_tables)

        return bemt.bemt_axial(prop, pitch, o, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                               tip_loss=tip_loss, mach_corr=mach_corr, output='long')

    def print_performance(c, p):
        print "Case: " + c[3]
        print "Thrust = " + str(sum(p[0]))
        print "Power = " + str(sum(p[1]))
        print "dP = " + str(p[1])
        print "Cd = " + str(p[2])
        print "Cl = " + str(p[3])
        print "Cl/Cd = " + str(p[3]/p[2])
        print "Re = " + str(p[-1])
        print "\n"

    perf = [get_performance(case[0], case[1], case[2], tip_loss) for case in cases2run]
    for i in xrange(len(cases2run)):
        print_performance(cases2run[i], perf[i])

    plt.figure(1)
    l = []
    for i, case in enumerate(cases2run):
        plt.plot(r, case[1])
        l.append("Case " + str(i+1))
    plt.xlabel("radial station, r")
    plt.ylabel("chord")
    plt.legend(l)

    plt.figure(2)
    for i, case in enumerate(cases2run):
        plt.plot(r, case[2]*360/2/np.pi)
        l.append("Case " + str(i+1))
    plt.xlabel("radial station, r")
    plt.ylabel("twist, degrees")
    plt.legend(l)

    plt.show()
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



