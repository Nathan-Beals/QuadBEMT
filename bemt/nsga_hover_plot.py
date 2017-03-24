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
    #radius = 0.193
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
    mach_corr = False
    alt = 0

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

    # DA4002
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

    # (OPT # 1) NSGA 10 pt 5000p 500g Reallowed=All T=4.61 xinit=base tiploss c/Rmax=0.6
    chord = np.array([0.11470868, 0.21117422, 0.30252118, 0.39582241, 0.42498534, 0.41241736, 0.32928823, 0.23259937,
                      0.17993009, 0.11100332])
    twist = np.array([0.55653743, 0.45683294, 0.38657942, 0.39545834, 0.3270524, 0.3053103, 0.27289751, 0.23934658,
                      0.20193655, 0.034115])
    omega = 5286.48230542 * 2*np.pi/60
    opt1 = (omega, chord, twist, 'NSGA10e_5000p_500g_Reallowed=all_T=4.61_tiploss')

    # (OPT # 2) SLSQP Opt of DA4002 10 element. Reallow=All, T=4.61, Tiploss
    chord = np.array([0.11994618, 0.19944238, 0.2781612, 0.33993244, 0.33206394, 0.28471642, 0.22285041, 0.18446874,
                      0.12841128, 0.01537926])
    twist = np.array([0.51724241, 0.4268415, 0.39644603, 0.33147393, 0.31163335, 0.28957352, 0.25814615, 0.24327558,
                      0.20237492, 0.00653631])
    omega = 5943.55347904 * 2*np.pi/60
    opt2 = (omega, chord, twist, 'SLSQP of DA4002_10e_Reallowed=all_T=4.61_tiploss')

    # (OPT # 3) SLSQP Opt of OPT # 1, Reallow=All, T=4.61, Tiploss
    chord = np.array([0.12000001, 0.22000002, 0.30478926, 0.4047892, 0.41469428, 0.35391398, 0.31830687, 0.21970192,
                      0.16508119, 0.11102533])
    twist = np.array([0.59057564, 0.4386199, 0.42691716, 0.35621971, 0.33793345, 0.30917999, 0.27913501, 0.25924289,
                      0.21518868, 0.04164781])
    omega = 5286.15517496 * 2*np.pi/60
    opt3 = (omega, chord, twist, 'SLSQP of '+opt1[3])

    # (OPT # 4 NSGA 10 pt 5000p 500g Reallowed=All T=5.5 tiploss alt=4000m c/R_max = 0.6
    chord = np.array([0.08653958, 0.18652537, 0.28627588, 0.38037405, 0.47607346, 0.46306288, 0.36738897, 0.31966668,
                      0.29608263, 0.26964711])
    twist = np.array([0.00173169, 0.17526932, 0.34294626, 0.368123, 0.35016916, 0.2610559, 0.24640237, 0.23694115,
                      0.18387004, 0.08772996])
    omega = 2496.5002026 * 2*np.pi/60
    opt4 = (omega, chord, twist, 'Carrol w/ c/Rmax=0.6')

    # (OPT # 5 NSGA 10pt 5000p 500g Reallowed=All T=5.5 tiploss alt=4000m c/R_max = 0.4
    chord = np.array([0.08527494, 0.18305166, 0.28003061, 0.37723735, 0.39791253, 0.38551352, 0.30090341, 0.20180039,
                      0.14683824, 0.06464163])
    twist = np.array([0.00114182, 0.16927357, 0.33062163, 0.36775474, 0.32781283, 0.30751832, 0.27097053, 0.24367986,
                      0.21398367, 0.10101229])
    omega = 2600.92469603 * 2*np.pi/60
    opt5 = (omega, chord, twist, 'Carrol w c/Rmax=0.4')

    # (OPT # 6 SLSQP optimization of OPT # 4
    chord = np.array([0.08658001, 0.18656597, 0.28631939, 0.3795963, 0.47300912, 0.45898628, 0.36328787, 0.31595086,
                      0.2926356, 0.26650533])
    twist = np.array([0.00173185, 0.17528111, 0.34303816, 0.3674234, 0.35081355, 0.26280001, 0.24807192, 0.23788198,
                      0.18459802, 0.08882935])
    omega = 2496.49300497 * 2*np.pi/60
    opt6 = (omega, chord, twist, 'SLSQP of '+opt4[3])

    # (OPT # 7 SLSQP optimization of OPT # 5
    chord = np.array([0.08527585, 0.18305299, 0.28003252, 0.37723067, 0.39788281, 0.38547729, 0.30086853, 0.20177347,
                      0.14681765, 0.06462935])
    twist = np.array([0.00114194, 0.16927492, 0.33062551, 0.36773377, 0.32782983, 0.30753222, 0.27097839, 0.24371668,
                      0.21400531, 0.10102211])
    omega = 2600.92449624 * 2*np.pi/60
    opt7 = (omega, chord, twist, 'SLSQP of '+opt5[3])

    # OPT # 8 NSGA 10 pt 5000p 250g Reallowed=All T=4.61 xinit=base tiploss c/Rmax=0.4
    chord = np.array([0.11963147, 0.21948599, 0.31325439, 0.39324083, 0.38238664, 0.31217141, 0.27986064, 0.22972166,
                      0.18952114, 0.10681847])
    twist = np.array([0.58226601, 0.41130811, 0.42698568, 0.36964679, 0.3344689, 0.28206486, 0.26203572, 0.26472563,
                      0.17337351, 0.04165566])
    omega = 574.732670
    opt8 = (omega, chord, twist, 'NSGA_10pt_5000p_250g_T=4.61_c/Rmax=0.4')

    # OPT # 9 NSGA 10 pt 5000p 250g Reallowed=All T=4.61 xinit=base tiploss c/Rmax=0.3
    chord = np.array([0.10277016, 0.19955381, 0.25497272, 0.25707051, 0.21844591, 0.26933938, 0.21127499, 0.18302549,
                      0.10701967, 0.0085749])
    twist = np.array([0.33535394, 0.41390913, 0.35894807, 0.32476601, 0.2650822, 0.2682921, 0.25730397, 0.21412539,
                      0.19555273, 0.07476783])
    omega = 681.173949
    opt9 = (omega, chord, twist, 'NSGA_10pt_5000p_250g_T=4.61_c/Rmax=0.3')

    # OPT # 10 SLSQP optimization of OPT # 8
    chord = np.array([0.1198, 0.21958079, 0.31958069, 0.38922467, 0.36971416, 0.29535462, 0.25967067, 0.21561951,
                      0.18162215, 0.10694134])
    twist = np.array([0.60635079, 0.43899644, 0.40277763, 0.34947613, 0.32675009, 0.29279386, 0.27404534, 0.25758374,
                      0.21384189, 0.03945741])
    omega = 5488.38927605 * 2*np.pi/60
    opt10 = (omega, chord, twist, 'SLSQP of '+opt8[3])

    # OPT # 11 SLSQP optimization of OPT # 9
    chord = np.array([0.10639693, 0.20386555, 0.26456196, 0.26621299, 0.23369806, 0.26408675, 0.21013131, 0.1687175, 0.10552476, 0.00674878])
    twist = np.array([0.32810083, 0.3527571, 0.33500691, 0.31100313, 0.25881002, 0.2835209, 0.25685164, 0.23318679, 0.18140129, 0.04924958])
    omega = 6504.83447976 * 2*np.pi/60
    opt11 = (omega, chord, twist, 'SLSQP of '+opt9[3])

    # OPT # 12 NSGA 10pt 5000p 500g Reallowed=All T=4.61 xinit=base tiploss c/Rmax=0.3
    chord = np.array([0.11958375, 0.21820237, 0.2761888, 0.29582408, 0.2670838, 0.25332996, 0.21333354, 0.17139233,
                      0.09985169, 0.00114413])
    twist = np.array([0.34160335, 0.39718607, 0.37649986, 0.35195434, 0.26725317, 0.25970378, 0.24134294, 0.19232958,
                      0.17713771, 0.02991615])
    omega = 6574.95092927 * 2*np.pi/60
    opt12 = (omega, chord, twist, 'NSGA_10pt_5000p_500g_T=4.61_c/Rmax=0.3')

    # OPT # 13 SLSQP optimization of OPT # 12
    chord = np.array([0.11956503, 0.21708697, 0.274589, 0.29301011, 0.26383857, 0.24929219, 0.20911063, 0.16797264,
                      0.09876164, 0.0011191])
    twist = np.array([0.34100869, 0.39172876, 0.37041507, 0.34754131, 0.26809427, 0.26250597, 0.2461286, 0.19832232,
                      0.17977324, 0.03239786])
    omega = 6574.96318233 * 2*np.pi/60
    opt13 = (omega, chord, twist, 'SLSQP of '+opt12[3])

    # OPT # 14 NSGA 10 pt 5000p 500g Reallowed=All T=4.61 xinit=base tiploss c/Rmax=0.2
    chord = np.array([0.11979988, 0.16749867, 0.19254645, 0.19463329, 0.19635013, 0.19565444, 0.19250374, 0.16106175,
                      0.09234817, 0.00055174])
    twist = np.array([0.64485005, 0.48662595, 0.37141195, 0.30465231, 0.26518135, 0.26743179, 0.24975756, 0.2205021,
                      0.17419048, 0.02354747])
    omega = 6960.37743102 * 2*np.pi/60
    opt14 = (omega, chord, twist, 'NSGA_10pt_5000p_500g_T=4.61_c/Rmax=0.2')

    # OPT # 15 NSGA 10 pt 5000p 500g Reallowed=All T=4.61 xinit=base tiploss c/Rmax=0.4
    chord = np.array([0.11382682, 0.21274313, 0.30536217, 0.3878935, 0.37140305, 0.30969498, 0.25016208, 0.19799559,
                      0.16869266, 0.10448205])
    twist = np.array([0.60318707, 0.42869181, 0.42709282, 0.35208938, 0.32339698, 0.29094724, 0.2505386, 0.24947303,
                      0.20516782, 0.03193417])
    omega = 5615.24908926 * 2*np.pi/60
    opt15 = (omega, chord, twist, 'NSGA_10pt_5000p_500g_T=4.61_c/Rmax=0.4')

    # OPT # 16 NSGA 10 pt 5000p 500g Reallowed=All T=5.5 alt=4000 tiploss c/Rmax=0.3
    chord = np.array([8.92386048e-02, 1.73000845e-01, 2.70523039e-01, 2.71542807e-01, 2.78749355e-01, 2.36866151e-01,
                      2.04103526e-01, 1.37456074e-01, 8.68094589e-02, 1.05601135e-04])
    twist = np.array([0.00161645, 0.15105685, 0.28791442, 0.31577392, 0.28644651, 0.27418749, 0.24854514, 0.21812646,
                      0.19802027, 0.14972058])
    omega = 3184.41320387 * 2*np.pi/60
    opt16 = (omega, chord, twist, 'Carrol w c/Rmax=0.3')

    # OPT # 17 NSGA 10 pt 5000p 500g Reallowed=All T=5.5 alt=4000 tiploss c/Rmax=0.2
    chord = np.array([0.0518402, 0.100106, 0.19271316, 0.19690797, 0.18563632, 0.17909439, 0.14834613, 0.13261967,
                      0.10188759, 0.00190945])
    twist = np.array([0.00116852, 0.16823478, 0.31311048, 0.27888512, 0.25888298, 0.24650332, 0.22996178, 0.22292014,
                      0.17389608, 0.02976441])
    omega = 3542.74531607 * 2*np.pi/60
    opt17 = (omega, chord, twist, 'Carroll w c/Rmax=0.2')

    # OPT # 18 SLSQP optimization of OPT # 15
    chord = np.array([0.1198, 0.2198, 0.3198, 0.39998273, 0.36354119, 0.31198366, 0.25089796, 0.19514172, 0.13489857,
                      0.10450839])
    twist = np.array([0.60193722, 0.43850153, 0.41493225, 0.35232941, 0.32482648, 0.29761368, 0.27135973, 0.24905772,
                      0.2040977, 0.02962213])
    omega = 5615.17802792 * 2*np.pi/60
    opt18 = (omega, chord, twist, 'SLSQP of '+opt15[3])

    # OPT # 19 SLSQP optimization of OPT # 14
    chord = np.array([0.11979995, 0.18609334, 0.19681079, 0.19770371, 0.19843832, 0.19814065, 0.19409736, 0.15772375,
                      0.09672598, 0.00023612])
    twist = np.array([0.55991373, 0.41820809, 0.36433488, 0.3015608, 0.26968276, 0.2302029, 0.25095793, 0.23438399,
                      0.17403859, 0.01622641])
    omega = 6960.70746602 * 2*np.pi/60
    opt19 = (omega, chord, twist, 'SLSQP of '+opt14[3])

    # OPT # 20 SLSQP optimization of OPT # 16
    chord = np.array([9.00000025e-02, 1.89999642e-01, 2.89999645e-01, 2.99999093e-01, 2.44846538e-01, 2.00366933e-01,
                      1.64549740e-01, 1.36163105e-01, 9.28467808e-02, 2.29006079e-08])
    twist = np.array([0.5745666, 0.40003502, 0.34104694, 0.3222615, 0.28755732, 0.26060618, 0.23983699, 0.22522367,
                      0.19440473, 0.13845183])
    omega = 3183.71953258 * 2*np.pi/60
    opt20 = (omega, chord, twist, 'SLSQP of '+opt16[3])

    # OPT # 21 SLSQP optimization of OPT # 17
    chord = np.array([0.05189348, 0.10199043, 0.19297108, 0.1938842, 0.17973795, 0.1731011, 0.14317486, 0.12792554,
                      0.09828124, 0.00184187])
    twist = np.array([0.03350517, 0.19980878, 0.33337885, 0.29418723, 0.26880152, 0.25118232, 0.230463, 0.2182604,
                      0.1766044, 0.03184142])
    omega = 3542.73176324 * 2*np.pi/60
    opt21 = (omega, chord, twist, 'SLSQP of '+opt17[3])

    cases2run = (opt3, opt18, opt13, opt19, base)

    def get_performance(o, c, t):
        chord_meters = c * radius
        prop = propeller.Propeller(t, chord_meters, radius, n_blades, r, y, dr, dy, airfoils=airfoils,
                                   Cl_tables=Cl_tables, Cd_tables=Cd_tables)

        return bemt.bemt_axial(prop, pitch, o, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                               tip_loss=tip_loss, mach_corr=mach_corr, output='long', alt=alt)

    def print_performance(c, p):
        print "Case: " + c[3]
        print "Thrust = " + str(sum(p[0]))
        print "Power = " + str(sum(p[1]))
        print "dP = " + str(p[1])
        print "Cd = " + str(p[2])
        print "Cl = " + str(p[3])
        print "Cl/Cd = " + str(p[3]/p[2])
        print "Re = " + str(p[-3])
        print "CT = " + str(p[-2])
        print "Pprofile = " + str(p[-2])
        print "Pinduced = " + str(p[-1])
        print r
        print "\n"

    perf = [get_performance(case[0], case[1], case[2]) for case in cases2run]
    for i in xrange(len(cases2run)):
        print_performance(cases2run[i], perf[i])

    l_case1 = [r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.6}$',
               r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.4}$',
               r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.3}$',
               r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.2}$',
               r'$\mathrm{DA4002}$']

    l_case2 = [r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.6}$',
               r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.4}$',
               r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.3}$',
               r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.2}$']

    # plt.figure(1)
    # l = []
    # for i, case in enumerate(cases2run):
    #     plt.plot(r, case[1])
    #     l.append("Case " + str(i+1))
    # plt.xlabel("radial station, r", fontsize=18)
    # plt.ylabel(r'$\mathrm{c/R}$', fontsize=18)
    # plt.legend(l)
    # plt.tick_params(axis='both', which='major', labelsize=14)
    # plt.tick_params(axis='both', which='minor', labelsize=14)
    #
    # plt.figure(2)
    # for i, case in enumerate(cases2run):
    #     plt.plot(r, case[2]*360/2/np.pi)
    #     l.append("Case " + str(i+1))
    # plt.xlabel("radial station, r", fontsize=18)
    # plt.ylabel(r'$\beta,\,\mathrm{degrees}$', fontsize=18)
    # plt.legend(l)
    # plt.tick_params(axis='both', which='major', labelsize=14)
    # plt.tick_params(axis='both', which='minor', labelsize=14)

    plt.figure(3)
    plt.plot(r, opt3[1], 'k*-', markerfacecolor='white')
    plt.plot(r, opt18[1], 'kv-', markerfacecolor='white')
    plt.plot(r, opt13[1], 'ks-', markerfacecolor='white')
    plt.plot(r, opt19[1], 'ko-', markerfacecolor='white')
    plt.plot(r, base[1], 'kD-', markerfacecolor='white')
    plt.xlabel("radial station, r", fontsize=18)
    plt.ylabel(r'$\mathrm{c/R}$', fontsize=18)
    plt.legend(l_case1)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.ylim([0.0, 0.5])
    plt.grid()

    plt.figure(4)
    plt.plot(r, opt3[2]*360/2/np.pi, 'k*-', markerfacecolor='white')
    plt.plot(r, opt18[2]*360/2/np.pi, 'kv-', markerfacecolor='white')
    plt.plot(r, opt13[2]*360/2/np.pi, 'ks-', markerfacecolor='white')
    plt.plot(r, opt19[2]*360/2/np.pi, 'ko-', markerfacecolor='white')
    plt.plot(r, base[2]*360/2/np.pi, 'kD-', markerfacecolor='white')
    plt.xlabel("radial station, r", fontsize=18)
    plt.ylabel(r'$\beta,\,\mathrm{deg}$', fontsize=18)
    plt.legend(l_case1)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.grid()

    plt.figure(5)
    plt.plot(r, opt6[1], 'k*-', markerfacecolor='white')
    plt.plot(r, opt7[1], 'kv-', markerfacecolor='white')
    plt.plot(r, opt16[1], 'ks-', markerfacecolor='white')
    plt.plot(r, opt21[1], 'ko-', markerfacecolor='white')
    plt.xlabel("radial station, r", fontsize=18)
    plt.ylabel(r'$\mathrm{c/R}$', fontsize=18)
    plt.legend(l_case2)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.ylim([0.0, 0.5])
    plt.grid()

    plt.figure(6)
    plt.plot(r, opt6[2]*360/2/np.pi, 'k*-', markerfacecolor='white')
    plt.plot(r, opt7[2]*360/2/np.pi, 'kv-', markerfacecolor='white')
    plt.plot(r, opt16[2]*360/2/np.pi, 'ks-', markerfacecolor='white')
    plt.plot(r, opt21[2]*360/2/np.pi, 'ko-', markerfacecolor='white')
    plt.xlabel("radial station, r", fontsize=18)
    plt.ylabel(r'$\beta,\,\mathrm{deg}$', fontsize=18)
    plt.legend(l_case2)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.grid()

    plt.show()

if __name__ == '__main__':
    main()



