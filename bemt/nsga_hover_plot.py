import numpy as np
import matplotlib.pyplot as plt
import propeller
import quadrotor
import bemt
import unit_conversion
import lookup_table
import aero_coeffs
import trim


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

    # DA4002
    omega = 5943.0 * 2*np.pi/60
    chord = np.array([0.1198, 0.1128, 0.1436, 0.1689, 0.1775, 0.1782, 0.1773, 0.1782, 0.1790, 0.1787, 0.1787, 0.1786,
                      0.1785, 0.1790, 0.1792, 0.1792, 0.1692, 0.0154])
    chord = np.array([chord[i] for i in [0, 2, 4, 6, 8, 10, 12, 14, 15, 17]])
    twist = np.array([42.481, 44.647, 41.154, 37.475, 34.027, 30.549, 27.875, 25.831, 23.996, 22.396, 21.009, 19.814,
                      18.786, 17.957, 17.245, 16.657, 13.973, 2.117]) * 2 * np.pi / 360
    twist = np.array([twist[i] for i in [0, 2, 4, 6, 8, 10, 12, 14, 15, 17]])
    base = (omega, chord, twist, 'DA4002')

    # 5k pop, 8 gen, 5 azi_elem, 9.6 in prop, 12.455 N weight
    chord = np.array([0.09024478, 0.18743201, 0.2589227, 0.2692717, 0.36417748, 0.39572612, 0.46280335, 0.48632645,
                      0.449613, 0.50439701])
    twist = np.array([6.37001388e-05, 1.07349719e-01, 2.64969331e-01, 3.46340529e-01, 3.53158041e-01, 2.99419272e-01,
                      1.36439696e-01, 1.45026520e-01, 2.09532374e-01, 1.83957955e-01])
    omega = 4294.52427027 * 2*np.pi/60
    case1 = (omega, chord, twist, '8gen_5azi')

    cases2run = (case1,)

    def get_performance(o, c, t):
        """
        :param o: hover rotational speed in rad/s
        :param c: chord distribution normalized by the rotor radius
        :param t: twist distribution in radians
        :return:
        """
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
        return sum(dT_h), P_h, T_ff, H_ff, P_ff, alpha_trim, omega_trim

    def print_performance(c, p):
        print "Case: " + c[3]
        print "Hover thrust = " + str(p[0])
        print "Hover power  = " + str(p[1])
        print "FF thrust    = " + str(p[2])
        print "FF drag      = " + str(p[3])
        print "FF power     = " + str(p[4])
        print "alpha trim   = " + str(p[5])
        print "omega trim   = " + str(p[6])

    perf_list = [get_performance(case[0], case[1], case[2]) for case in cases2run]
    for i, perf in enumerate(perf_list):
        print_performance(cases2run[i], perf)

    # l_case1 = [r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.6}$',
    #            r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.4}$',
    #            r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.3}$',
    #            r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.2}$',
    #            r'$\mathrm{DA4002}$']
    #
    # l_case2 = [r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.6}$',
    #            r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.4}$',
    #            r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.3}$',
    #            r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.2}$',
    #            r'$\mathrm{Carroll}$']
    #
    # l_case2_nocarroll = [r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.6}$',
    #                      r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.4}$',
    #                      r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.3}$',
    #                      r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.2}$']
    #
    # color_str = ['k*-', 'kv-', 'ks-', 'ko-', 'kD-']

    # plt.figure(1)
    # l = []
    # for i, case in enumerate(cases2run):
    #     plt.plot(r, perf[i][-1]/1000., color_str[i], markerfacecolor='white')
    #     l.append(l_case2_nocarroll[i])
    # plt.xlabel("radial station, r", fontsize=18)
    # plt.ylabel(r'$\mathrm{Re}\,\mathrm{(10}^\mathrm{3}\mathrm{)}$', fontsize=18)
    # plt.legend(l, loc='upper left')
    # plt.tick_params(axis='both', which='major', labelsize=14)
    # plt.tick_params(axis='both', which='minor', labelsize=14)
    # plt.xlim([0.0, 1.0])
    # #
    # # plt.figure(2)
    # # for i, case in enumerate(cases2run):
    # #     plt.plot(r, case[2]*360/2/np.pi)
    # #     l.append("Case " + str(i+1))
    # # plt.xlabel("radial station, r", fontsize=18)
    # # plt.ylabel(r'$\beta,\,\mathrm{degrees}$', fontsize=18)
    # # plt.legend(l)
    # # plt.tick_params(axis='both', which='major', labelsize=14)
    # # plt.tick_params(axis='both', which='minor', labelsize=14)
    #
    # plt.figure(3)
    # plt.plot(r, opt3[1], 'k*-', markerfacecolor='white')
    # plt.plot(r, opt18[1], 'kv-', markerfacecolor='white')
    # plt.plot(r, opt13[1], 'ks-', markerfacecolor='white')
    # plt.plot(r, opt19[1], 'ko-', markerfacecolor='white')
    # plt.plot(r, base[1], 'kD-', markerfacecolor='white')
    # plt.xlabel("radial station, r", fontsize=18)
    # plt.ylabel(r'$\mathrm{c/R}$', fontsize=18)
    # plt.legend(l_case1, loc='upper left')
    # plt.tick_params(axis='both', which='major', labelsize=14)
    # plt.tick_params(axis='both', which='minor', labelsize=14)
    # plt.ylim([0.0, 0.5])
    # plt.xlim([0.0, 1.0])
    # plt.grid()
    #
    # plt.figure(4)
    # plt.plot(r, opt3[2]*360/2/np.pi, 'k*-', markerfacecolor='white')
    # plt.plot(r, opt18[2]*360/2/np.pi, 'kv-', markerfacecolor='white')
    # plt.plot(r, opt13[2]*360/2/np.pi, 'ks-', markerfacecolor='white')
    # plt.plot(r, opt19[2]*360/2/np.pi, 'ko-', markerfacecolor='white')
    # plt.plot(r, base[2]*360/2/np.pi, 'kD-', markerfacecolor='white')
    # plt.xlabel("radial station, r", fontsize=18)
    # plt.ylabel("twist, deg", fontsize=18)
    # plt.legend(l_case1)
    # plt.tick_params(axis='both', which='major', labelsize=14)
    # plt.tick_params(axis='both', which='minor', labelsize=14)
    # plt.xlim([0.0, 1.0])
    # plt.grid()
    #
    # carroll_chord = [0.07, 0.12, 0.13, 0.14, 0.2, 0.185, 0.16, 0.141, 0.099, 0.06]
    # carroll_twist = [0.0, 17., 22., 16.8, 11.5, 10., 9.5, 9., 8, 1.5]
    # carroll_r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    #
    # plt.figure(5)
    # plt.plot(r, opt6[1], 'k*-', markerfacecolor='white')
    # plt.plot(r, opt7[1], 'kv-', markerfacecolor='white')
    # plt.plot(r, opt16[1], 'ks-', markerfacecolor='white')
    # plt.plot(r, opt21[1], 'ko-', markerfacecolor='white')
    # plt.plot(carroll_r, carroll_chord, 'k^-', markerfacecolor='white')
    # plt.xlabel("radial station, r", fontsize=18)
    # plt.ylabel(r'$\mathrm{c/R}$', fontsize=18)
    # plt.legend(l_case2, loc='upper left')
    # plt.tick_params(axis='both', which='major', labelsize=14)
    # plt.tick_params(axis='both', which='minor', labelsize=14)
    # plt.ylim([0.0, 0.5])
    # plt.xlim([0.0, 1.0])
    # plt.grid()
    #
    # plt.figure(6)
    # plt.plot(r, opt6[2]*360/2/np.pi, 'k*-', markerfacecolor='white')
    # plt.plot(r, opt7[2]*360/2/np.pi, 'kv-', markerfacecolor='white')
    # plt.plot(r, opt16[2]*360/2/np.pi, 'ks-', markerfacecolor='white')
    # plt.plot(r, opt21[2]*360/2/np.pi, 'ko-', markerfacecolor='white')
    # plt.plot(carroll_r, carroll_twist, 'k^-', markerfacecolor='white')
    # plt.xlabel("radial station, r", fontsize=18)
    # plt.ylabel("twist, deg", fontsize=18)
    # plt.legend(l_case2, loc='upper right')
    # plt.tick_params(axis='both', which='major', labelsize=14)
    # plt.tick_params(axis='both', which='minor', labelsize=14)
    # plt.xlim([0.0, 1.0])
    # plt.grid()
    #
    # plt.show()

if __name__ == '__main__':
    main()



