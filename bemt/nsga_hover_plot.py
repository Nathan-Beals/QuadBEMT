import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import propeller
import quadrotor
import bemt
import old_ff_bemt
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

    ####################################### Hover opts using obj geogonrecheck #######################################
    # Hover opt 100 gen, 1000 pop, 12.455 N weight, 9.6 in prop
    chord = np.array([0.11019156, 0.20021497, 0.30021309, 0.38116296, 0.41117733, 0.34220286, 0.31976733, 0.3051378,
                      0.2461775, 0.16171859])
    twist = np.array([0.4114017, 0.27568074, 0.42871935, 0.37570118, 0.32963633, 0.31803996, 0.20201996, 0.2792073,
                      0.23461276, 0.11703326])
    omega = 3765.2288005 * 2*np.pi/60
    hover1 = (omega, chord, twist, 'hover100gen')

    # Hover opt 500 gen, 1000 pop, 12.455 N weight, 9.6 in prop
    chord = np.array([0.11923604, 0.2168746, 0.31540216, 0.39822882, 0.42919, 0.35039799, 0.3457828, 0.28567224, 0.23418368, 0.13502483])
    twist = np.array([0.45316866, 0.38457724, 0.38225075, 0.34671967, 0.33151445, 0.28719111, 0.25679667, 0.25099005, 0.19400679, 0.10926302])
    omega = 3811.03596674 * 2*np.pi/60
    hover2 = (omega, chord, twist, 'hover500gen')

    ####################################### v_inf = 4 cases ###########################################################

    # 300 pop, 500 gen, 9.6, 12.455 c_max = 0.6, dFz[-1] = 0.
    chord = np.array([0.11968948,  0.21937987,  0.3182458,  0.41816814,  0.46347988, 0.37271461,  0.33244453,
                       0.28706656,  0.22163902,  0.12675553])
    twist = np.array([0.34646128,  0.459587,  0.35604962,  0.34112132,  0.27325105, 0.26739119,  0.25000727,
                       0.19915033,  0.15261586, -0.01996389])
    omega = 4131.4987072* 2*np.pi/60
    case13 = (omega, chord, twist, 'case13')

    # 300 pop 700 gen, 9.6, 12.455, c_max = 0.6, dFz[np.argwhere(r>0.97)] = 0.
    chord = np.array([0.11967194, 0.21647242, 0.31550402, 0.4138184, 0.50218225, 0.48664647, 0.40026445, 0.30073221, 0.23394046, 0.13399367])
    twist = np.array([4.64143128e-01, 3.07029419e-01, 3.67075807e-01, 3.39429422e-01, 3.03971603e-01, 2.73597532e-01,
                      2.52153328e-01, 2.36219090e-01, 1.74550201e-01, 5.51910537e-05])
    omega = 3860.12156849* 2*np.pi/60
    case14 = (omega, chord, twist, 'case14')

    # 300 pop 700 gen 9.6 12.455 c_max = 5., dchord = 0.2 (aka totally unrestricted)
    chord = np.array([0.11911694, 0.31906466, 0.51906199, 0.71471174, 0.91466236, 1.02746187, 1.04465427, 0.95496686,
                      0.84514798, 0.66650262])
    twist = np.array([0.56684028, 0.5399532, 0.37894101, 0.39985625, 0.33103878, 0.26743957, 0.22936254, 0.18269451,
                      0.12498615, -0.04437022])
    omega = 3290.65852292* 2*np.pi/60
    case16 = (omega, chord, twist, 'case16')

    # 300 pop 900 gen c_max = 0.6 no eff_aoa restriction (same as case 13/14 but with 900 gen)
    chord = np.array([ 0.11976461,  0.21674507,  0.31632813,  0.41492634,  0.50915059, 0.49242833,  0.40604632,  0.30740114,  0.23627649,  0.13816009])
    twist = np.array([ 0.4572405 ,  0.30366688,  0.36459117,  0.33370541,  0.30116492, 0.27257351,  0.25236179,  0.23919269,  0.16748809, -0.00700692])
    omega = 3860.17077361* 2*np.pi/60
    case17 = (omega, chord, twist, 'case17')

    # 300 pop 500 gen c_max = 0.6 (should be the same as case13, doublechecking)
    chord = np.array([0.11979654,  0.21971113,  0.31768566,  0.41654234,  0.51464254, 0.55749166,  0.50071197,  0.59703481,  0.54194578,  0.53298076])
    twist = np.array([0.40945918,  0.31129442,  0.39272177,  0.3771759 ,  0.26534598, 0.30003162,  0.21020565,  0.19785936,  0.1209569 , -0.05170223])
    omega = 3797.95349418* 2*np.pi/60
    case18 = (omega, chord, twist, 'case18')

    # 300 pop 500 gen c_max = 5.0, 9.6, 12.455
    chord = np.array([ 0.1197982 ,  0.21965552,  0.31948792,  0.417718  ,  0.51206442, 0.55625209,  0.47125151,  0.55641091,  0.59843211,  0.53537045])
    twist = np.array([ 0.17666704,  0.26000819,  0.32499536,  0.3212985 ,  0.32544299, 0.27487363,  0.25279071,  0.19595582,  0.12082239, -0.0520852 ])
    omega = 3790.0222464* 2*np.pi/60
    case19 = (omega, chord, twist, 'case19')

    # 300 pop 600 gen c_max = 0.6 9.6in 12.455
    chord = np.array([0.11852411, 0.21848123, 0.31645576, 0.41543821, 0.51131543, 0.55250242, 0.49552009, 0.59441046,
                      0.5379174, 0.54141083])
    twist = np.array([0.4099726, 0.30146709, 0.3815585, 0.36601236, 0.28159023, 0.29640811, 0.20990898, 0.19841352,
                      0.12222948, -0.05159499])
    omega = 3802.42024921* 2*np.pi/60
    case20 = (omega, chord, twist, 'case20')

    # 300 pop 500 gen 9.6 12.455 c_max = 0.6 dL[-1] = 0.
    chord = np.array([0.1184386 ,  0.21841442,  0.31825553,  0.41738728,  0.46364181, 0.48349788,  0.38852228,
                      0.30379532,  0.23479559,  0.16161961])
    twist = np.array([0.39676715,  0.43825037,  0.34666288,  0.32926906,  0.33052567, 0.27683669,  0.25849766,
                      0.23486052,  0.18739048,  0.02909696])
    omega = 3816.66551278* 2*np.pi/60
    case21 = (omega, chord, twist, 'case21')

    # 300 pop 600 gen 9.6 12.455 c_max = 0.6 dL[-1] = 0
    chord = np.array([0.11922888,  0.21920113,  0.31900133,  0.4169311, 0.46194409, 0.48770574, 0.39343732, 0.3026776,
                      0.23230887, 0.14034282])
    twist = np.array([0.39725249,  0.44465958,  0.35515042,  0.327588, 0.32424015, 0.27099176, 0.25265596, 0.23940818,
                      0.19681692, 0.02419742])
    omega = 3816.60942736* 2*np.pi/60
    case22 = (omega, chord, twist, 'case22')

    # 300 pop 700 gen 9.6 12.455 c_max = 0.6 dL[-1] = 0.
    chord = np.array([0.1184555, 0.21844679, 0.31842026, 0.41231438, 0.46518974, 0.48510572, 0.39062965, 0.29569951,
                      0.22929197, 0.1387922])
    twist = np.array([0.39653315, 0.43499077, 0.34511864, 0.32267895, 0.32432533, 0.27068676, 0.25234788, 0.2396517,
                      0.20244573, 0.029898])
    omega = 3825.66098485* 2*np.pi/60
    case23 = (omega, chord, twist, 'case23')

    # 300 pop 500 gen 9.6 12.455 c_max = 5. dL[-1] = 0.
    chord = np.array([0.11961287, 0.21874277, 0.31825403, 0.41824134, 0.51428782, 0.6040633, 0.67264788, 0.69568852,
                      0.6078391 , 0.53489543])
    twist = np.array([0.50077726, 0.37533403, 0.3591938, 0.31893474, 0.3021634, 0.28436863, 0.20427255, 0.1664158,
                      0.12501316, 0.00465653])
    omega = 3824.3909652* 2*np.pi/60
    case25 = (omega, chord, twist, 'case25')

    # 300 pop 900 gen 9.6 12.455 c_max = 0.6 dL[-1] = 0
    chord = np.array([0.11979653, 0.21978687, 0.31978536, 0.41962515, 0.47246546, 0.49262045, 0.39833725, 0.30804084,
                      0.23791201, 0.13931635])
    twist = np.array([0.39651539, 0.44928119, 0.35940906, 0.32635117, 0.32290322, 0.26977521, 0.25144634, 0.23789696,
                      0.20069098, 0.02809155])
    omega = 3800.79067169* 2*np.pi/60
    case26 = (omega, chord, twist, 'case26')

    # 1000 pop 500 gen 9.6 12.455 c_max = 0.6 dL[-1] = 0.
    chord = np.array([0.1196944, 0.21770649, 0.31717874, 0.41113149, 0.49500415, 0.46606351, 0.36794244, 0.30387524, 0.23692036, 0.13915433])
    twist = np.array([0.42859819, 0.43809706, 0.31993687, 0.32664475, 0.31451747, 0.26411071, 0.25858693, 0.2296344, 0.21382795, 0.08701286])
    omega = 3841.54315307* 2*np.pi/60
    case27 = (omega, chord, twist, 'case27')

    # 300 p 500 g 9.6 12.455 c_max = 0.4 dL[-1] = 0.
    chord = np.array([0.1178891, 0.21722331, 0.31551886, 0.35617918, 0.39895906, 0.39711467, 0.36235461, 0.26470141,
                      0.20096506, 0.14229118])
    twist = np.array([0.34490299, 0.28131579, 0.35149267, 0.33722407, 0.26644888, 0.26262655, 0.22241742, 0.23250699,
                      0.1860713, 0.01269026])
    omega = 4157.06784876* 2*np.pi/60
    case28 = (omega, chord, twist, 'case28')

    # 300p 500g 9.6 12.455 c_max = 0.3 dL[-1] = 0.
    chord = np.array([0.09235703, 0.19232431, 0.28593088, 0.29529476, 0.29956715, 0.29818168, 0.24476486, 0.23165514,
                      0.14004843, 0.11367408])
    twist = np.array([0.19878182, 0.22847754, 0.33446732, 0.28813435, 0.27464274, 0.23748252, 0.22892405, 0.1737965,
                      0.16414805, 0.06473553])
    omega = 4777.12987159* 2*np.pi/60
    case29 = (omega, chord, twist, 'case29')

    # 300p 500g 9.6 12.455 c_max = 0.2 dL[-1] = 0.
    chord = np.array([0.11974937, 0.1988364, 0.19845687, 0.18358676, 0.19946962, 0.18837444, 0.15840381, 0.149858,
                      0.10195713, 0.00286152])
    twist = np.array([0.36476035, 0.1984514, 0.31483518, 0.23364511, 0.25401822, 0.18120287, 0.212223, 0.19510032,
                      0.15543919, 0.07810287])
    omega = 5704.80803182* 2*np.pi/60
    case30 = (omega, chord, twist, 'case30')

    # 300 p 500g 9.6 12.455 c_max 0.6 UNIFORM INFLOW
    chord = np.array([0.10391889, 0.18799523, 0.27614229, 0.37562084, 0.45449965, 0.48123414, 0.58121353, 0.58874668,
                      0.56506879, 0.47377971])
    twist = np.array([0.04186956, 0.19010911, 0.21211266, 0.36218881, 0.36069853, 0.2886057, 0.23230848, 0.26068052,
                      0.23741129, 0.13761006])
    omega = 3421.57237585* 2*np.pi/60



    ######################################## "Forward flight optimization case" ################################
    # 300 pop 500 gen cmax=0.6 9.6in 12.455 times = 500 s ff, 1 s hover
    chord = np.array([0.11979663, 0.21965327, 0.31961028, 0.41960177, 0.51513626, 0.59998041, 0.59952532, 0.59950103,
                      0.57028389, 0.47230215])
    twist = np.array([0.47928485, 0.30532313, 0.31923982, 0.29516638, 0.28386677, 0.23312418, 0.17156788, 0.12472249,
                      0.06041558, -0.03990621])
    omega = 4545.36006558 * 2*np.pi/60
    ff1 = (omega, chord, twist, 'ff1')

    ####################################### "Hover optimization case (1 second in ff) ###########################
    chord = np.array([0.11302328, 0.21298767, 0.31283867, 0.36367498, 0.28204006, 0.3687244, 0.3455101, 0.25136865,
                      0.23114779, 0.13553282])
    twist = np.array([0.20810447, 0.35332739, 0.42475926, 0.33883661, 0.27949311, 0.31701892, 0.28471683, 0.25705176,
                      0.18672455, 0.0562048 ])
    omega = 3903.31609145* 2*np.pi/60
    hover3 = (omega, chord, twist, 'hover3')


    cases2run = (case26, case27, case23)

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
        print "Hover omega  = " + str(c[0]*60/2/np.pi)
        print "Hover thrust = " + str(p[0])
        print "Hover power  = " + str(p[1])
        print "FF thrust    = " + str(p[2])
        print "FF drag      = " + str(p[3])
        print "FF power     = " + str(p[4])
        print "alpha trim   = " + str(p[5])
        print "omega trim   = " + str(p[6]*60/2/np.pi)
        energy = mission_time[0]*p[1] + mission_time[1]*p[4]
        print "mission energy = " + str(energy)

    perf_list = [get_performance(case[0], case[1], case[2]) for case in cases2run]
    for i, perf in enumerate(perf_list):
        print_performance(cases2run[i], perf)

    color_str = ['k*-', 'kv-', 'ks-', 'ko-', 'kD-']

    plt.figure(1)
    l = []
    for i, case in enumerate(cases2run):
        plt.plot(r, case[1], color_str[i], markerfacecolor='white')
        l.append(case[-1])
    plt.xlabel("radial station, r", fontsize=18)
    plt.ylabel(r'$\mathrm{c/R}$', fontsize=18)
    plt.legend(l, loc='upper left')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 0.5])
    plt.grid()

    plt.figure(2)
    for i, case in enumerate(cases2run):
        plt.plot(r, case[2]*360/2/np.pi, color_str[i], markerfacecolor='white')
        l.append(case[-1])
    plt.xlabel("radial station, r", fontsize=18)
    plt.ylabel(r'$\beta,\,\mathrm{degrees}$', fontsize=18)
    plt.legend(l)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.grid()

    # data = perf_list[0][-1]
    # psi = np.linspace(0, 2*np.pi, n_azi_elements)
    # fig = plt.figure()
    # ax = Axes3D(fig)
    # r, th = np.meshgrid(r, psi)
    # plt.subplot(projection="polar")
    # plt.pcolormesh(th, r, data)
    # plt.colorbar()
    # #plt.pcolormesh(th, z, r)
    # plt.plot(psi, r, color='k', ls='none')
    # plt.grid()

    plt.show()

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
    #

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



##################################### Cases from old code versions that weren't working properly #################
# # 5k pop, 10 gen, 9.6in, 12.445 N
# chord = np.array([0.0728792, 0.15130816, 0.24288949, 0.34110123, 0.41877626, 0.50875582, 0.50219216, 0.48461647,
#                   0.52715237, 0.4960476])
# twist = np.array([0.38480448, 0.22860905, 0.35666174, 0.2094294, 0.33141239, 0.34117174, 0.26380779, 0.20220869,
#                   0.05000656, 0.0440292])
# omega = 3998.0227186 * 2*np.pi / 60
# case1 = (omega, chord, twist, 'case1')
#
# # 1000 pop, 50 gen, 9.6 in 12.455 N
# chord = np.array([0.11661063, 0.20755812, 0.28378134, 0.35876439, 0.43670325, 0.53006448, 0.5603994, 0.56488541,
#                   0.58958217, 0.572237])
# twist = np.array([4.85285561e-01, 3.12838874e-01, 3.19185491e-01, 3.81078639e-01, 2.50041106e-01, 2.51847501e-01,
#                   1.67247053e-01, 1.69140731e-01, 1.59020072e-01, -3.04859048e-04])
# omega = 4057.36371477 * 2*np.pi/60
# case2 = (omega, chord, twist, 'case2')
#
# # 300 pop, 15 gen, 9.6 in, 12.455 N
# chord = np.array([0.11515094,  0.21356443,  0.2860178 ,  0.30196478,  0.37619129, 0.39739965,  0.30603692,
#                   0.30535553,  0.22078035,  0.12639688])
# twist = np.array([0.43410038,  0.3728036 ,  0.29237855,  0.34650477,  0.38169184, 0.31704846,  0.2876572,
#                   0.21076985,  0.20945838,  0.07399843])
# omega = 3993.95692615 * 2*np.pi/60
# case3 = (omega, chord, twist, 'case3')
#
# # 300 pop, 500 gen, 9.6 in, 12.455 N, c_max = 0.6, no c_tip_max
# chord = np.array([0.11333318,  0.20714241,  0.30688952,  0.40460222,  0.50182865, 0.55843652,  0.58074275,
#                   0.59463337,  0.56960695,  0.57561376])
# twist = np.array([0.43713452,  0.37610116,  0.28035933,  0.34455469,  0.29613552, 0.23117021,  0.22684089,
#                   0.16601488,  0.11331936,  0.01326311])
# omega = 3995.62869717 * 2*np.pi/60
# case4 = (omega, chord, twist, 'case4')
#
# # 300 pop, 500 gen, 9.6 in, 12.455 N, c_max = 0.6, c_tip_max = 0.25
# chord = np.array([0.11939826, 0.21928043, 0.31654785, 0.41624897, 0.46007747, 0.44671773, 0.3593443, 0.31909281,
#                   0.22553974, 0.1303472])
# twist = np.array([0.29267164, 0.3625469, 0.28505877, 0.3554414, 0.30034785, 0.26578228, 0.25055408, 0.20833788,
#                   0.1816911, 0.03696819])
# omega = 3952.30288602* 2*np.pi/60
# case5 = (omega, chord, twist, 'case5')
#
# # 300 pop, 500 gen, 9.6 in, 12.455 N, c_max = 5. c_tip_max = 5.
# chord = np.array([0.11974794,  0.21912069,  0.31812644,  0.41705443,  0.51620713, 0.61611694,  0.7018769,
#                   0.72303464,  0.74187665,  0.82406679])
# twist = np.array([0.27600223,  0.34765329,  0.4189613,  0.31230476,  0.31763595, 0.27352963,  0.18479028,
#                   0.18164489,  0.1301707,  0.01280449])
# omega = 3769.03133552* 2*np.pi/60
# case6 = (omega, chord, twist, 'case6')
#
# # 300 pop, 500 gen, 9.6 in, 12.455 N, c_max = 0.6, c_tip_max = 0.4
# chord = np.array([0.1197414, 0.21951981, 0.31937234, 0.41567804, 0.51567661, 0.45103824, 0.43736838, 0.34144522,
#                   0.27313875, 0.25051002])
# twist = np.array([0.37309394, 0.48562675, 0.3347374, 0.31337371, 0.33948754, 0.2676794, 0.26660352, 0.23390463,
#                   0.18071602, 0.0097657])
# omega = 3757.58133102* 2*np.pi/60
# case7 = (omega, chord, twist, 'case7')
#
# # 300 pop, 500 gen, 9.6 in, 12.455 N, c_max = 0.4, c_tip_max = 0.25
# chord = np.array([0.11917634, 0.21906258, 0.31715898, 0.39373037, 0.38919806, 0.35545338, 0.31881517, 0.2784707,
#                   0.19671025, 0.14658203])
# twist = np.array([0.43826917, 0.44700471, 0.37168474, 0.28746841, 0.27842648, 0.27734615, 0.23716972, 0.20165707,
#                   0.1808409, 0.00838674])
# omega = 4193.86615577 * 2*np.pi/60
# case8 = (omega, chord, twist, 'case8')
#
# # 300 pop, 700 gen, 9.6 in, 12.455 N, c_max = 0.6, c_tip_max = 0.25
# chord = np.array([0.11697466, 0.21676395, 0.31676083, 0.41367311, 0.45804153, 0.46873766, 0.36884461, 0.32855824,
#                   0.22906588, 0.13536223])
# twist = np.array([0.30572188, 0.37334749, 0.29231291, 0.37587646, 0.30210497, 0.282821, 0.25690732, 0.22940615,
#                   0.19667221, 0.0233916 ])
# omega = 3808.17064497 * 2*np.pi/60
# case9 = (omega, chord, twist, 'case9')
#
# # 300 pop, 500 gen, 9.6 in, 12.455 N, c_max = 0.6, c_tip_max = 5., F[np.isnan(F)] = 0.99999999
# chord = np.array([0.11810595, 0.21795989, 0.31654866, 0.41634815, 0.50951731, 0.56384868, 0.59367274, 0.56612156,
#                   0.59401312, 0.59259751])
# twist = np.array([0.46647849, 0.44827483, 0.36121518, 0.25119389, 0.26948641, 0.29670099, 0.19325489, 0.20519149,
#                   0.11709039, 0.0138717])
# omega = 3895.22520691 * 2*np.pi/60
# case10 = (omega, chord, twist, 'case10')
#
# # 300 pop, 500 gen, 9.6 in, 12.455 N, c_max = 100., dchord_upper/lower = 0.2/-0.2
# chord = np.array([0.118906,  0.31877433,  0.51877154,  0.71871369,  0.91774941, 1.11350574,  1.3134808, 1.5095779,
#                   1.70781305, 1.90730912])
# twist = np.array([0.40483507,  0.4336412,  0.3700821,  0.36742537,  0.32636117, 0.31836313,  0.21515102, 0.15104657,
#                   0.12512185, 0.04945353])
# omega = 3262.61032833 * 2*np.pi/60
# case11 = (omega, chord, twist, 'case11')
#
# # 300 pop, 700 gen, 9.6, 12.455, c_max = 0.3, c_tip = 0.25
# chord = np.array([0.10561703, 0.1985356, 0.29743209, 0.26521539, 0.2981148, 0.29414362, 0.27450984, 0.22486947,
#                   0.19163144, 0.12626986])
# twist = np.array([0.45382792, 0.38828041, 0.32190524, 0.27353347, 0.24903454, 0.24899439, 0.21955375, 0.20586851,
#                   0.15679694, 0.0129808])
# omega = 4647.55026389 * 2*np.pi/60
# case12 = (omega, chord, twist, 'case12')