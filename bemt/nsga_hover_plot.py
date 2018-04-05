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
    time_in_hover = 0     # Time in seconds
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
    hover1 = (omega, chord, twist, v_inf, 'hover100gen')

    # Hover opt 500 gen, 1000 pop, 12.455 N weight, 9.6 in prop
    chord = np.array([0.11923604, 0.2168746, 0.31540216, 0.39822882, 0.42919, 0.35039799, 0.3457828, 0.28567224, 0.23418368, 0.13502483])
    twist = np.array([0.45316866, 0.38457724, 0.38225075, 0.34671967, 0.33151445, 0.28719111, 0.25679667, 0.25099005, 0.19400679, 0.10926302])
    omega = 3811.03596674 * 2*np.pi/60
    hover2 = (omega, chord, twist, v_inf, 'hover500gen')

    # 300p 1100 g 9.6 12.455
    chord = np.array([0.11843164, 0.21662496, 0.31015779, 0.37157618, 0.41826784, 0.40700329, 0.32408767, 0.23006844,
                      0.2031113, 0.11133428])
    twist = np.array([0.59411782, 0.5060331, 0.41068162, 0.38048604, 0.33646003, 0.32812993, 0.27305065, 0.22576999,
                      0.19610176, 0.05793477])
    omega = 3787.65021182 * 2*np.pi/60
    hover3 = (omega, chord, twist, v_inf, 'hover 300p 1100g')

    ####################################### v_inf = 4 cases ###########################################################
    ##################### c_max = 0.6 cases ####################################################
    # 300/500 9.6 12.455 cmax=0.6 dL[-1] = 0 in ff only pMut = 0.4
    chord = np.array([0.11847653, 0.21846786, 0.31817481, 0.41473315, 0.51108386, 0.44612053, 0.38662758, 0.29248156,
                      0.20472331, 0.13940618])
    twist = np.array([0.4724454, 0.4628394, 0.33895807, 0.2780318, 0.2888433, 0.28405327, 0.22388795, 0.23338565,
                      0.18997092, 0.06096851])
    omega = 3984.70883238 * 2*np.pi/60
    v_inf = 4.0
    g500_cmax06_pmut04 = (omega, chord, twist, v_inf, '500g cmax0.6 pmut0.4')

    # 300/500 9.6 12.455 cmax=0.6 dL[-1] = 0 in ff only pMut = 0.3
    chord = np.array([0.11490002, 0.2147637, 0.31385954, 0.41377984, 0.50030368, 0.51857104, 0.59275099, 0.59282674,
                      0.58309676, 0.53218839])
    twist = np.array([0.53360365, 0.44923495, 0.28449606, 0.27628026, 0.3258491, 0.28968536, 0.2093388, 0.20263226,
                      0.1248437, 0.02248061])
    omega = 3823.06586764 * 2*np.pi/60
    v_inf = 4.0
    g500_cmax06_pmut03 = (omega, chord, twist, v_inf, '500g cmax0.6 pmut0.3')

    # 300/500 9.6 12.455 cmax=0.6 dL[-1] = 0 in ff only pMut = 0.25
    chord = np.array([0.11913403, 0.21902208, 0.31710577, 0.41478176, 0.44725369, 0.41141313, 0.34371539, 0.28588431,
                      0.20752417, 0.13072573])
    twist = np.array([0.55005268, 0.47791234, 0.36330333, 0.26885214, 0.28296035, 0.29166506, 0.25659206, 0.22793896,
                      0.19599009, 0.02891601])
    omega = 3942.18295126 * 2*np.pi/60
    v_inf = 4.0
    g500_cmax06_pmut025 = (omega, chord, twist, v_inf, '500g cmax0.6 pmut0.25')

    # 300/900 9.6 12.455 cmax=0.6 dL[-1] = 0 in ff only pMut = 0.3
    chord = np.array([0.11938031, 0.2192617, 0.31926072, 0.41698046, 0.50348454, 0.51664837, 0.59075465, 0.59079784,
                      0.58601493, 0.53233059])
    twist = np.array([0.51577126, 0.44931604, 0.2817431, 0.3231065, 0.31115893, 0.28789016, 0.20750974, 0.19441336,
                      0.11639786, 0.02451927])
    omega = 3855.99949332 * 2*np.pi/60
    v_inf = 4.0
    g900_cmax06_pMut03 = (omega, chord, twist, v_inf, '900g cmax0.6 pMut0.3')

    # 300/500 9.6 12.455 cmax=0.6 dL[-1] = 0 in ff only
    chord = np.array([0.11966637, 0.21900429, 0.31680873, 0.41423814, 0.46098928, 0.46460738, 0.36895834, 0.30434101, 0.20548289, 0.13486721])
    twist = np.array([0.33294302, 0.40815838, 0.31966698, 0.3442134, 0.32092452, 0.26738233, 0.27019631, 0.23245922, 0.20029737, 0.09543751])
    omega = 3811.29319425 * 2*np.pi/60
    v_inf = 4.0
    g500_cmax06 = (omega, chord, twist, v_inf, '500g cmax0.6')

    # 300/700 9.6 12.455 cmax = 0.6 dL[-1] = 0 in ff only
    chord = np.array([0.11966429, 0.21945035, 0.31908323, 0.41651265, 0.46036133, 0.46428148, 0.3695573, 0.30493997,
                      0.20608201, 0.13484094])
    twist = np.array([0.33294302, 0.40738511, 0.31856103, 0.34310745, 0.31982095, 0.26627876, 0.26909274, 0.23335901,
                      0.20090297, 0.09586083])
    omega = 3811.32410953 * 2*np.pi/60
    v_inf = 4.0
    g700_cmax06 = (omega, chord, twist, v_inf, '700g cmax0.6')

    # 300/900 9.6 12.455 cmax=0.6 dL[-1] = in in ff only
    chord = np.array([0.11966637, 0.21902615, 0.31677977, 0.41441104, 0.47634846, 0.46720008, 0.3715869, 0.30466464,
                      0.20706076, 0.13528156])
    twist = np.array([0.33294302, 0.40823459, 0.31943387, 0.34398026, 0.32068984, 0.26726229, 0.2700911, 0.23229424,
                      0.19400504, 0.04906076])
    omega = 3811.43278646 * 2*np.pi/60
    v_inf = 4.0
    g900_cmax06 = (omega, chord, twist, v_inf, '900g cmax0.6')

    # 300/1100 9.6 12.455 cmax = 0.6 dL[-1] = 0 in ff only
    chord = np.array([0.11966637, 0.21945255, 0.3191826, 0.41661201, 0.46042987, 0.46404804, 0.36930386, 0.30469107,
                      0.20583294, 0.13521735])
    twist = np.array([0.33294302, 0.40815866, 0.31931257, 0.34385899, 0.32057012, 0.26702793, 0.2698419, 0.23210482,
                      0.19992483, 0.04420424])
    twist = np.array([0.49, 0.40815866, 0.31931257, 0.34385899, 0.32057012, 0.26702793, 0.2698419, 0.23210482,
                      0.19992483, 0.04420424])
    omega = 3811.29467393 * 2*np.pi/60
    v_inf = 4.0
    g1100_cmax06 = (omega, chord, twist, v_inf, '1100g cmax0.6')

    # 1000/500 9.6 12.455 cmax = 0.6 dL[-1] in ff only
    chord = np.array([0.11620593, 0.21605971, 0.31560612, 0.41512654, 0.50377189, 0.50106609, 0.40115501, 0.30785641, 0.23370032, 0.1352783])
    twist = np.array([0.44441691, 0.40087482, 0.37615744, 0.30005601, 0.3003249, 0.26681709, 0.24090271, 0.24450708, 0.19593291, 0.02180752])
    omega = 3812.76567949 * 2*np.pi/60
    v_inf = 4.0
    g500_p1000_cmax06 = (omega, chord, twist, v_inf, '500g 1000p cmax0.6')

    ################################ c_max = 0.4 cases #################################33
    # 300/500 9.6 12.455 cmax = 0.4 dL[-1] = 0 in ff only pmut = 0.3
    chord = np.array([0.11461409, 0.21153649, 0.30793342, 0.39988646, 0.38276784, 0.34831652, 0.30505384, 0.26654542,
                      0.1890812, 0.12839928])
    twist = np.array([0.40672333, 0.34149858, 0.37968077, 0.28476059, 0.29889374, 0.27699235, 0.26043234, 0.24346932,
                      0.2077853, 0.06095271])
    omega = 4016.96958559 * 2*np.pi/60
    v_inf = 4.0
    g500_cmax04_pmut03 = (omega, chord, twist, v_inf, '500g cmax0.6 pmut0.3')

    # 300/500 9.6 12.455 cmax = 0.4 dL[-1] = 0 in ff only
    chord = np.array([0.11915047, 0.21849478, 0.31788447, 0.38013234, 0.33717869, 0.37913264, 0.29024744, 0.25839581,
                      0.18715189, 0.12003223])
    twist = np.array([0.35875486, 0.30029531, 0.36665571, 0.29751113, 0.27435418, 0.26480507, 0.20252785, 0.21349503,
                      0.18433429, 0.04248182])
    omega = 4292.27791073 * 2*np.pi/60
    v_inf = 4.0
    g500_cmax04 = (omega, chord, twist, v_inf, '500g cmax0.4')

    # 300/700 9.6 12.455 cmax = 0.4 dL[-1] = 0 in ff only
    chord = np.array([0.11928499, 0.21862929, 0.31801899, 0.38026685, 0.33731321, 0.37926716, 0.29038195, 0.25853033,
                      0.1872864, 0.12016675])
    twist = np.array([0.35875486, 0.30029531, 0.36665571, 0.29769775, 0.27449779, 0.26494719, 0.20266997, 0.21363714,
                      0.1844764, 0.01673915])
    omega = 4292.27791073 * 2*np.pi/60
    v_inf = 4.0
    g700_cmax04 = (omega, chord, twist, v_inf, '700g cmax0.4')

    # 300/900 9.6 12.455 cmax = 0.4 dL[-1] = 0 in ff only
    chord = np.array([0.11929206, 0.21863637, 0.31802607, 0.38027393, 0.33732028, 0.37927423, 0.29038903, 0.25853741,
                      0.18729348, 0.12017382])
    twist = np.array([0.35874391, 0.30028436, 0.3714652, 0.30370865, 0.2805517, 0.27100259, 0.20872537, 0.20964449,
                      0.16755424, 0.00394849])
    omega = 4292.27791073 * 2*np.pi/60
    v_inf = 4.0
    g900_cmax04 = (omega, chord, twist, v_inf, '900g cmax0.4')

    # 300/1100 9.6 12.455 cmax = 0.4 dL[-1] = 0 in ff only
    chord = np.array([0.11928697, 0.21863128, 0.31802097, 0.38026884, 0.33731519, 0.37926914, 0.29038393, 0.25853231,
                      0.18728838, 0.12016873])
    twist = np.array([0.35875486, 0.30029531, 0.3666576, 0.30198579, 0.27877385, 0.26922488, 0.20694766, 0.20176056,
                      0.18324988, 0.0087248])
    omega = 4292.27791073 * 2*np.pi/60
    v_inf = 4.0
    g1100_cmax04 = (omega, chord, twist, v_inf, '1100g cmax0.4')

    # 300/1100 9.6 12.455 cmax = 0.4 dL[-1] = 0 in ff only pmut = 0.3
    chord = np.array([0.11790448, 0.20916995, 0.30206876, 0.39116066, 0.38369303, 0.34657802, 0.33523101, 0.29802089,
                      0.22755861, 0.12910791])
    twist = np.array([0.41070948, 0.33380376, 0.37425208, 0.30071594, 0.31253589, 0.26232999, 0.25380047, 0.24380831,
                      0.17474276, 0.05712327])
    omega = 3984.81388187 * 2*np.pi/60
    v_inf = 4.0
    g1100_cmax04_pmut03 = (omega, chord, twist, v_inf, '1100g cmax0.4 pmut0.3')


    ############################## c_max = 0.3 cases ############################
    # 300/500 9.6 12.455 cmax = 0.3 dL[-1] = 0 in ff only pMut = 0.3
    chord = np.array([0.11584583, 0.20992076, 0.27050657, 0.29213525, 0.27971303, 0.27731739, 0.23130046, 0.2228163,
                      0.13167057, 0.03200796])
    twist = np.array([0.42504762, 0.36559227, 0.31472603, 0.25987081, 0.25742304, 0.24685365, 0.19613074, 0.19764744,
                      0.18090024, 0.16342444])
    omega = 4848.35704233 * 2*np.pi/60
    v_inf = 4.0
    g500_cmax03_pmut03 = (omega, chord, twist, v_inf, '300/500 cmax0.3 pMut0.3')

    # 300/500 9.6 12.455 cmax = 0.3 dL[-1] = 0 in ff only pMut = 0.25
    chord = np.array([0.09803243, 0.19788593, 0.29439608, 0.25822472, 0.29112695, 0.28688834, 0.23469024, 0.2356359,
                      0.16944293, 0.11012919])
    twist = np.array([0.40496939, 0.31186813, 0.34677394, 0.28614832, 0.27585495, 0.2528887, 0.22575301, 0.17916481,
                      0.15323425, 0.04918692])
    omega = 4673.82174614 * 2*np.pi/60
    v_inf = 4.0
    g500_cmax03_pmut025 = (omega, chord, twist, v_inf, '300/500 cmax0.3 pMut0.25')

    # 300/500 9.6 12.455 cmax = 0.3 dL[-1] = 0 in ff only
    chord = np.array([0.11546019, 0.21433848, 0.2852159, 0.29467352, 0.23804887, 0.2915607, 0.25084762, 0.2104048,
                      0.1894686, 0.1248765])
    twist = np.array([0.59284111, 0.42972097, 0.28802944, 0.30903828, 0.24766486, 0.25454117, 0.20488933, 0.18536425,
                      0.13422406, -0.00468971])
    omega = 4814.18056096 * 2*np.pi/60
    v_inf = 4.0
    g500_cmax03 = (omega, chord, twist, v_inf, '500g cmax0.3')

    # 300/700 9.6 12.455 cmax = 0.3 dL[-1] = 0 in ff only
    chord = np.array([0.11757832, 0.21497178, 0.28584921, 0.29908738, 0.24563781, 0.29811567, 0.25321042, 0.2127676,
                      0.17992469, 0.12679849])
    twist = np.array([0.59262868, 0.41820139, 0.29148644, 0.30144142, 0.2555974, 0.2588507, 0.22751418, 0.20648936,
                      0.15780436, -0.01070732])
    omega = 4638.71361799 * 2*np.pi/60
    v_inf = 4.0
    g700_cmax03 = (omega, chord, twist, v_inf, '700g cmax0.3')

    # 300/1100 9.6 12.455 cmax = 0.3 dL[-1] = 0 in ff only
    chord = np.array([0.11668114, 0.21436524, 0.28524266, 0.2978438, 0.24439383, 0.296601, 0.25625421, 0.21000138,
                      0.17114882, 0.11125641])
    chord = np.array([0.11668114, 0.21436524, 0.28524266, 0.2978438, 0.2981, 0.296601, 0.25625421, 0.21000138,
                      0.17114882, 0.11125641])
    twist = np.array([0.59260295, 0.41823816, 0.2911823, 0.29818251, 0.25210241, 0.2553557, 0.22401918, 0.20289526,
                      0.154201, 0.03183863])
    omega = 4633.37690235 * 2*np.pi/60
    v_inf = 4.0
    g1100_cmax03 = (omega, chord, twist, v_inf, '1100g cmax0.3')

    ############################### c_max = 0.2 cases #############################################
    # 300/900 9.6 12.455 cmax = 0.2 dL[-1] = 0 in ff only
    chord = np.array([0.1197121, 0.17193583, 0.1987021, 0.19814803, 0.19581124, 0.17572306, 0.18155826, 0.15528053,
                      0.10534952, 0.00548909])
    twist = np.array([0.37537592, 0.27609915, 0.34542111, 0.23700681, 0.2126302, 0.20417805, 0.18649412, 0.18049172,
                      0.15703695, 0.1039155])
    omega = 5709.50549273 * 2*np.pi/60
    v_inf = 4.0
    g900_cmax02 = (omega, chord, twist, v_inf, '900g cmax0.2')

    # 300/1100 9.6 12.455 cmax = 0.2 dL[-1] = 0 in ff only
    chord = np.array([0.11771472, 0.16970633, 0.19623047, 0.19561193, 0.19388567, 0.1759398, 0.17140728, 0.16272408,
                      0.10491537, 0.00500641])
    twist = np.array([0.36833199, 0.27075819, 0.2496, 0.23122227, 0.21652866, 0.20627188, 0.20134915, 0.17757133,
                      0.15418487, 0.04510659])
    omega = 5773.07817791 * 2*np.pi/60
    v_inf = 4.0
    g1100_cmax02 = (omega, chord, twist, v_inf, '1100g cmax0.2')

    # 300/500 9.6 12.455 cmax = 0.2 dL[-1] = 0 in ff only pMut = 0.3
    chord = np.array([0.11463782, 0.1448912, 0.19949987, 0.1968132, 0.18420485, 0.18413558, 0.15749732, 0.14232675,
                      0.10599133, 0.00628147])
    twist = np.array([0.39664003, 0.36818507, 0.19837917, 0.24127658, 0.20634473, 0.17108318, 0.18002921, 0.16616339,
                      0.13558457, 0.28454228])
    omega = 6133.32992167 * 2*np.pi/60
    v_inf = 4.0
    g500_cmax02_pmut03 = (omega, chord, twist, v_inf, '500g cmax0.2 pmut0.3')

    ######################################## "Forward flight optimization case" ################################
    # 300 pop 500 gen cmax=0.6 9.6in 12.455 times = 500 s ff, 1 s hover
    chord = np.array([0.11979663, 0.21965327, 0.31961028, 0.41960177, 0.51513626, 0.59998041, 0.59952532, 0.59950103,
                      0.57028389, 0.47230215])
    twist = np.array([0.47928485, 0.30532313, 0.31923982, 0.29516638, 0.28386677, 0.23312418, 0.17156788, 0.12472249,
                      0.06041558, -0.03990621])
    omega = 4545.36006558 * 2*np.pi/60
    v_inf = 4.0
    ff1 = (omega, chord, twist, v_inf, 'ff1')

    # 300 pop 500 gen cmax = 0.6 9.6 12.455 times = 500 s ff, 1 s hover
    chord = np.array([0.11978794, 0.21978419, 0.31954872, 0.41954125, 0.49704328, 0.50900264, 0.59325853, 0.59518854,
                      0.57548322, 0.47580079])
    twist = np.array([0.44303895, 0.28359202, 0.28218263, 0.29483998, 0.27398385, 0.24548995, 0.18876882, 0.13530204,
                      0.06926693, -0.04513989])
    omega = 4399.65264214 * 2*np.pi/60
    v_inf = 4.0
    ff2 = (omega, chord, twist, v_inf, 'ff2')

    # 300 p 900g cmax = 0.6 12.455 times = 500 s ff, 0.0001 s hover
    chord = np.array([0.11979085, 0.21978386, 0.31978378, 0.41978215, 0.51411119, 0.59926748, 0.59286911, 0.58544553, 0.5728953 , 0.47323209])
    twist = np.array([0.46984061, 0.29660078, 0.32316654, 0.30887717, 0.26967884, 0.22564501, 0.17985311, 0.13988012, 0.06513897, -0.03602759])
    omega = 5944.36943583 * 2*np.pi/60
    v_inf = 4.0
    ff3 = (omega, chord, twist, v_inf, 'ff3')

    # 300p 900g cmax = 0.6 9.6 12.455 times = 500s ff 0.001 s hover dL[-1] = 0 in ff only
    chord = np.array([0.11979398, 0.21979299, 0.31967194, 0.41948423, 0.51016813, 0.56521479, 0.59022916, 0.59994763,
                      0.59873803, 0.50478805])
    twist = np.array([0.30485815, 0.40221196, 0.31091777, 0.31019199, 0.29634665, 0.27401099, 0.20369896, 0.16051177,
                      0.08139644, -0.02436476])
    omega = 5477.241562 * 2*np.pi/60
    v_inf = 4.0
    ff4 = (omega, chord, twist, v_inf, 'ff limit')

    # 300/500 cmax =  0.6 9.6 12.455 vinf = 6. times = 500 s ff 0.001 s hover dL[-1] = 0 in ff only
    chord = np.array([0.11979624, 0.21976408, 0.31974935, 0.41967903, 0.51788119, 0.59944258, 0.59977493, 0.59995535,
                      0.58049264, 0.48084003])
    twist = np.array([0.36930057, 0.27143611, 0.2400825, 0.22244242, 0.18631815, 0.13399833, 0.11680813, 0.09436245,
                      0.06988677, -0.0189504])
    omega = 6609.02070934 * 2*np.pi/60
    v_inf = 6.0
    ff5 = (omega, chord, twist, v_inf, 'vinf = 6.')

    # 300/500 cmax = 0.6 9.6 12.455 v_inf = 2. times = 500 s ff 0.001 s hover dL[-1] in ff only
    chord = np.array([0.11979981, 0.21978519, 0.31978496, 0.41975751, 0.51958182, 0.47188887, 0.37544677, 0.29815084,
                      0.24787631, 0.15479815])
    twist = np.array([0.35107145, 0.37770426, 0.37681237, 0.37769542, 0.34357769, 0.31141169, 0.28524753, 0.26463772,
                      0.22771365, 0.06708662])
    omega = 4086.9824735 * 2*np.pi/60
    v_inf = 2.
    ff6 = (omega, chord, twist, v_inf, 'vinf = 2.')

    # 300/500 cmax = 0.6 9.6 12.455 v_inf = 1. times = 500 s ff 0.001 s hover dL[-1] in ff only
    chord = np.array([0.11963162, 0.21962728, 0.31956784, 0.41920849, 0.51582344, 0.42956998, 0.35785945, 0.27662849,
                      0.1923358, 0.14756468])
    twist = np.array([0.47608788, 0.53789671, 0.40021503, 0.38981403, 0.35069393, 0.32087076, 0.2943965, 0.2645475,
                      0.23755495, 0.06437609])
    omega = 6706.28019172 * 2*np.pi/60
    v_inf = 1.
    ff7 = (omega, chord, twist, v_inf, 'vinf = 1.')

    # 300/500 cmax = 0.6 9.6 12.455 v_inf = 0.1 times = 500 s ff 0.001 s hover dL[-1] in ff only
    chord = np.array([0.11944163, 0.21936326, 0.31933019, 0.40598688, 0.43794792, 0.37885534, 0.3230275, 0.26089917,
                      0.17860018, 0.11382533])
    twist = np.array([0.58327805, 0.44591819, 0.39134574, 0.34996863, 0.34628957, 0.31437616, 0.28297183, 0.25683314,
                      0.23517746, 0.06281419])
    omega = 5725.74774288 * 2*np.pi/60
    v_inf = 0.1
    ff8 = (omega, chord, twist, v_inf, 'vinf = 0.1')

    # 300/500 cmax = 0.6 9.6 12.455 v_inf = 0.01 times = 500 s ff 0.001 s hover dL[-1] in ff only
    chord = np.array([0.11975646, 0.21975493, 0.31788282, 0.38737281, 0.4448589, 0.36635049, 0.32984978, 0.26529463,
                      0.18422089, 0.14068652])
    twist = np.array([0.27465113, 0.43094561, 0.41543167, 0.37682445, 0.34445291, 0.31150413, 0.26999638, 0.26712771,
                      0.24018623, 0.08428127])
    omega = 7542.28176122 * 2*np.pi/60
    v_inf = 0.01
    ff9 = (omega, chord, twist, v_inf, 'vinf = 0.01')

    # 300/1100 cmax = 0.6 9.6 12.455 v_inf = 3. times = 500 s ff 0 s hover
    chord = np.array([0.11979983, 0.21979982, 0.31979978, 0.41978462, 0.5196591, 0.49930335, 0.39930724, 0.30866294,
                      0.2400338, 0.15735026])
    twist = np.array([0.45538258, 0.35159309, 0.32967499, 0.34870401, 0.31720839, 0.29386115, 0.27453983, 0.25594989,
                      0.2203728, 0.04584362])
    omega = 6802.30952139 * 2*np.pi/60
    v_inf = 3.
    g1100_ff_v3 = (omega, chord, twist, v_inf, 'ff v_inf=3')

    # 300/1100 cmax = 0.6 9.6 12.455 v_inf = 4. times = 500 s ff 0 s hover
    chord = np.array([0.11977906, 0.21910944, 0.3190903, 0.41883535, 0.51881464, 0.59997696, 0.59014961, 0.5819978,
                      0.57340307, 0.47359526])
    twist = np.array([0.44163572, 0.28092681, 0.31437201, 0.29890268, 0.26736614, 0.22448697, 0.18050267, 0.13971361,
                      0.06405428, -0.02630616])
    omega = 5866.74389564 * 2*np.pi/60
    v_inf = 4.
    g1100_ff_v4 = (omega, chord, twist, v_inf, 'ff v_inf=4')

    # 300/1100 cmax = 0.6 9.6 12.455 v_inf = 0.0001 times = 500 s ff 0 s hover "ff hover case"
    chord = np.array([0.11975438, 0.21975433, 0.31975248, 0.41970334, 0.51408979, 0.4442331, 0.34730048, 0.28507365,
                      0.1924555, 0.14695505])
    twist = np.array([0.68071731, 0.5726948, 0.4300643, 0.38847116, 0.36154393, 0.33356663, 0.30090916, 0.28245851,
                      0.2451, 0.07067683])
    omega = 6164.52809266 * 2*np.pi/60
    v_inf = 0.0001
    g1100_ffhover_v0001 = (omega, chord, twist, v_inf, 'ff hover v=.0001')

    # 300/1100 cmax = 0.6 9.6 12.455 v_inf = 0 times = 500 s ff 0 s hover "ff hover case"
    chord = np.array([0.11979986, 0.2196452, 0.31964355, 0.38886306, 0.42386373, 0.37300352, 0.32076483, 0.2561999,
                      0.17709234, 0.13731552])
    twist = np.array([0.43188647, 0.36893034, 0.42330566, 0.34895208, 0.33926483, 0.31310162, 0.28969622, 0.26322425,
                      0.23332532, 0.05879448])
    twist = np.array([0.456578, 0.439125, 0.42330566, 0.34895208, 0.33926483, 0.31310162, 0.28969622, 0.26322425,
                      0.23332532, 0.05879448])
    omega = 6939.10085731 * 2*np.pi/60
    v_inf = 0.
    g1100_ffhover_v0 = (omega, chord, twist, v_inf, 'ff hover v=0')

    ####################################### "Hover optimization case (0.001 second in ff) ###########################
    # 300/900 cmax = 0.6 9.6 12.455 times = 500 s hover, 0.001 s ff
    chord = np.array([0.11748975, 0.21223798, 0.3106284, 0.37953187, 0.4282579, 0.42372692, 0.34257662, 0.28768868,
                      0.22834561, 0.14139695])
    chord = np.array([0.11748975, 0.21223798, 0.3106284, 0.37953187, 0.4282579, 0.3828, 0.34257662, 0.28768868,
                      0.22834561, 0.14139695])
    twist = np.array([0.60065045, 0.54381023, 0.43543153, 0.35946649, 0.31969083, 0.32926997, 0.29780046, 0.26173905,
                      0.1959914, 0.02531718])
    omega = 3703.29999273 * 2*np.pi/60
    v_inf = 4.0
    hover7 = (omega, chord, twist, v_inf, 'hover limit')

    # 300/1100 cmax = 0.6 9.6 12.455 times = 500 s hover, 0 s ff, vinf = 0.0001
    chord = np.array([0.1191768, 0.21907918, 0.31029663, 0.32772206, 0.32907991, 0.40881134, 0.37411193, 0.28505367,
                      0.2155725, 0.13574664])
    twist = np.array([0.36923813, 0.34959923, 0.42455455, 0.39449948, 0.32269602, 0.30763444, 0.27348481, 0.22840145,
                      0.18712643, 0.02251118])
    omega = 3799.51206753 * 2*np.pi/60
    v_inf = 0.0001
    g1100_hoveropt_v0001 = (omega, chord, twist, v_inf, 'hover opt v=.0001')

    # 300/1100 cmax = 0.6 9.6 12.455 times = 500 s hover, 0 s ff, vinf = 0
    chord = np.array([0.10824171, 0.20823286, 0.30368035, 0.40258367, 0.4282665, 0.45646601, 0.38513469, 0.28630271,
                      0.21888292, 0.1448447])
    twist = np.array([0.1407348, 0.28393958, 0.38122109, 0.46555475, 0.32793555, 0.32259156, 0.29322366, 0.27949152,
                      0.20683867, 0.03624003])
    omega = 3560.16113105 * 2*np.pi/60
    v_inf = 0.
    g1100_hoveropt_v0 = (omega, chord, twist, v_inf, 'hover opt v=0')

    cases2run = (g1100_ffhover_v0, ff7, ff6, g1100_ff_v3, g1100_ff_v4)

    def get_performance(this_case):
        """
        :param o: hover rotational speed in rad/s
        :param c: chord distribution normalized by the rotor radius
        :param t: twist distribution in radians
        :return:
        """
        o, c, t, v, _ = this_case
        chord_meters = c * radius
        prop = propeller.Propeller(t, chord_meters, radius, n_blades, r, y, dr, dy, airfoils=airfoils,
                                   Cl_tables=Cl_tables, Cd_tables=Cd_tables)
        quad = quadrotor.Quadrotor(prop, vehicle_weight)
        ff_kwargs = {'propeller': prop, 'pitch': pitch, 'n_azi_elements': n_azi_elements, 'allowable_Re': allowable_Re,
                     'Cl_funs': Cl_funs, 'Cd_funs': Cd_funs, 'tip_loss': tip_loss, 'mach_corr': mach_corr, 'alt': alt,
                     'lift_curve_info_dict': lift_curve_info_dict}
        trim0 = np.array([alpha0, o])
        alpha_trim, omega_trim, converged = trim.trim(quad, v, trim0, ff_kwargs)
        T_ff, H_ff, P_ff = bemt.bemt_forward_flight(quad, pitch, omega_trim, alpha_trim, v, n_azi_elements, alt=alt,
                                                    tip_loss=tip_loss, mach_corr=mach_corr, allowable_Re=allowable_Re,
                                                    Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                                    lift_curve_info_dict=lift_curve_info_dict)

        dT_h, P_h = bemt.bemt_axial(prop, pitch, o, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                    tip_loss=tip_loss, mach_corr=mach_corr, alt=alt)
        return sum(dT_h), P_h, T_ff, H_ff, P_ff, alpha_trim, omega_trim

    def print_performance(c, p):
        print "Case: " + c[-1]
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

    perf_list = [get_performance(case) for case in cases2run]
    for i, perf in enumerate(perf_list):
        print_performance(cases2run[i], perf)

    color_str = ['k*-', 'kv-', 'ks-', 'ko-', 'kD-']*2

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
    plt.ylim([0.0, 0.7])
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

# 300 pop, 500 gen, 9.6, 12.455 c_max = 0.6, dFz[-1] = 0.
# chord = np.array([0.11968948,  0.21937987,  0.3182458,  0.41816814,  0.46347988, 0.37271461,  0.33244453,
#                    0.28706656,  0.22163902,  0.12675553])
# twist = np.array([0.34646128,  0.459587,  0.35604962,  0.34112132,  0.27325105, 0.26739119,  0.25000727,
#                    0.19915033,  0.15261586, -0.01996389])
# omega = 4131.4987072 * 2*np.pi/60
# case13 = (omega, chord, twist, 'case13')
#
# # 300 pop 700 gen, 9.6, 12.455, c_max = 0.6, dFz[np.argwhere(r>0.97)] = 0.
# chord = np.array([0.11967194, 0.21647242, 0.31550402, 0.4138184, 0.50218225, 0.48664647, 0.40026445, 0.30073221, 0.23394046, 0.13399367])
# twist = np.array([4.64143128e-01, 3.07029419e-01, 3.67075807e-01, 3.39429422e-01, 3.03971603e-01, 2.73597532e-01,
#                   2.52153328e-01, 2.36219090e-01, 1.74550201e-01, 5.51910537e-05])
# omega = 3860.12156849 * 2*np.pi/60
# case14 = (omega, chord, twist, 'case14')
#
# # 300 pop 700 gen 9.6 12.455 c_max = 5., dchord = 0.2 (aka totally unrestricted)
# chord = np.array([0.11911694, 0.31906466, 0.51906199, 0.71471174, 0.91466236, 1.02746187, 1.04465427, 0.95496686,
#                   0.84514798, 0.66650262])
# twist = np.array([0.56684028, 0.5399532, 0.37894101, 0.39985625, 0.33103878, 0.26743957, 0.22936254, 0.18269451,
#                   0.12498615, -0.04437022])
# omega = 3290.65852292 * 2*np.pi/60
# case16 = (omega, chord, twist, 'case16')
#
# # 300 pop 900 gen c_max = 0.6 no eff_aoa restriction (same as case 13/14 but with 900 gen)
# chord = np.array([ 0.11976461,  0.21674507,  0.31632813,  0.41492634,  0.50915059, 0.49242833,  0.40604632,  0.30740114,  0.23627649,  0.13816009])
# twist = np.array([ 0.4572405 ,  0.30366688,  0.36459117,  0.33370541,  0.30116492, 0.27257351,  0.25236179,  0.23919269,  0.16748809, -0.00700692])
# omega = 3860.17077361 * 2*np.pi/60
# case17 = (omega, chord, twist, 'case17')
#
# # 300 pop 500 gen c_max = 0.6 (should be the same as case13, doublechecking)
# chord = np.array([0.11979654,  0.21971113,  0.31768566,  0.41654234,  0.51464254, 0.55749166,  0.50071197,  0.59703481,  0.54194578,  0.53298076])
# twist = np.array([0.40945918,  0.31129442,  0.39272177,  0.3771759 ,  0.26534598, 0.30003162,  0.21020565,  0.19785936,  0.1209569 , -0.05170223])
# omega = 3797.95349418 * 2*np.pi/60
# case18 = (omega, chord, twist, 'case18')
#
# # 300 pop 500 gen c_max = 5.0, 9.6, 12.455
# chord = np.array([ 0.1197982 ,  0.21965552,  0.31948792,  0.417718  ,  0.51206442, 0.55625209,  0.47125151,  0.55641091,  0.59843211,  0.53537045])
# twist = np.array([ 0.17666704,  0.26000819,  0.32499536,  0.3212985 ,  0.32544299, 0.27487363,  0.25279071,  0.19595582,  0.12082239, -0.0520852 ])
# omega = 3790.0222464 * 2*np.pi/60
# case19 = (omega, chord, twist, 'case19')
#
# # 300 pop 600 gen c_max = 0.6 9.6in 12.455
# chord = np.array([0.11852411, 0.21848123, 0.31645576, 0.41543821, 0.51131543, 0.55250242, 0.49552009, 0.59441046,
#                   0.5379174, 0.54141083])
# twist = np.array([0.4099726, 0.30146709, 0.3815585, 0.36601236, 0.28159023, 0.29640811, 0.20990898, 0.19841352,
#                   0.12222948, -0.05159499])
# omega = 3802.42024921 * 2*np.pi/60
# case20 = (omega, chord, twist, 'case20')
#
# # 300 pop 500 gen 9.6 12.455 c_max = 0.6 dL[-1] = 0.
# chord = np.array([0.1184386 ,  0.21841442,  0.31825553,  0.41738728,  0.46364181, 0.48349788,  0.38852228,
#                   0.30379532,  0.23479559,  0.16161961])
# twist = np.array([0.39676715,  0.43825037,  0.34666288,  0.32926906,  0.33052567, 0.27683669,  0.25849766,
#                   0.23486052,  0.18739048,  0.02909696])
# omega = 3816.66551278 * 2*np.pi/60
# case21 = (omega, chord, twist, 'case21')
#
# # 300 pop 600 gen 9.6 12.455 c_max = 0.6 dL[-1] = 0
# chord = np.array([0.11922888,  0.21920113,  0.31900133,  0.4169311, 0.46194409, 0.48770574, 0.39343732, 0.3026776,
#                   0.23230887, 0.14034282])
# twist = np.array([0.39725249,  0.44465958,  0.35515042,  0.327588, 0.32424015, 0.27099176, 0.25265596, 0.23940818,
#                   0.19681692, 0.02419742])
# omega = 3816.60942736 * 2*np.pi/60
# case22 = (omega, chord, twist, 'case22')
#
# # 300 pop 700 gen 9.6 12.455 c_max = 0.6 dL[-1] = 0.
# chord = np.array([0.1184555, 0.21844679, 0.31842026, 0.41231438, 0.46518974, 0.48510572, 0.39062965, 0.29569951,
#                   0.22929197, 0.1387922])
# twist = np.array([0.39653315, 0.43499077, 0.34511864, 0.32267895, 0.32432533, 0.27068676, 0.25234788, 0.2396517,
#                   0.20244573, 0.029898])
# omega = 3825.66098485 * 2*np.pi/60
# case23 = (omega, chord, twist, 'case23')
#
# # 300 pop 500 gen 9.6 12.455 c_max = 5. dL[-1] = 0.
# chord = np.array([0.11961287, 0.21874277, 0.31825403, 0.41824134, 0.51428782, 0.6040633, 0.67264788, 0.69568852,
#                   0.6078391 , 0.53489543])
# twist = np.array([0.50077726, 0.37533403, 0.3591938, 0.31893474, 0.3021634, 0.28436863, 0.20427255, 0.1664158,
#                   0.12501316, 0.00465653])
# omega = 3824.3909652 * 2*np.pi/60
# case25 = (omega, chord, twist, 'case25')
#
# # 300 pop 900 gen 9.6 12.455 c_max = 0.6 dL[-1] = 0
# chord = np.array([0.11979653, 0.21978687, 0.31978536, 0.41962515, 0.47246546, 0.49262045, 0.39833725, 0.30804084,
#                   0.23791201, 0.13931635])
# twist = np.array([0.39651539, 0.44928119, 0.35940906, 0.32635117, 0.32290322, 0.26977521, 0.25144634, 0.23789696,
#                   0.20069098, 0.02809155])
# omega = 3800.79067169 * 2*np.pi/60
# case26 = (omega, chord, twist, 'case26')
#
# # 1000 pop 500 gen 9.6 12.455 c_max = 0.6 dL[-1] = 0.
# chord = np.array([0.1196944, 0.21770649, 0.31717874, 0.41113149, 0.49500415, 0.46606351, 0.36794244, 0.30387524, 0.23692036, 0.13915433])
# twist = np.array([0.42859819, 0.43809706, 0.31993687, 0.32664475, 0.31451747, 0.26411071, 0.25858693, 0.2296344, 0.21382795, 0.08701286])
# omega = 3841.54315307 * 2*np.pi/60
# case27 = (omega, chord, twist, 'case27')
#
# # 300 p 500 g 9.6 12.455 c_max = 0.4 dL[-1] = 0.
# chord = np.array([0.1178891, 0.21722331, 0.31551886, 0.35617918, 0.39895906, 0.39711467, 0.36235461, 0.26470141,
#                   0.20096506, 0.14229118])
# twist = np.array([0.34490299, 0.28131579, 0.35149267, 0.33722407, 0.26644888, 0.26262655, 0.22241742, 0.23250699,
#                   0.1860713, 0.01269026])
# omega = 4157.06784876 * 2*np.pi/60
# case28 = (omega, chord, twist, 'case28')
#
# # 300p 500g 9.6 12.455 c_max = 0.3 dL[-1] = 0.
# chord = np.array([0.09235703, 0.19232431, 0.28593088, 0.29529476, 0.29956715, 0.29818168, 0.24476486, 0.23165514,
#                   0.14004843, 0.11367408])
# twist = np.array([0.19878182, 0.22847754, 0.33446732, 0.28813435, 0.27464274, 0.23748252, 0.22892405, 0.1737965,
#                   0.16414805, 0.06473553])
# omega = 4777.12987159 * 2*np.pi/60
# case29 = (omega, chord, twist, 'case29')
#
# # 300p 500g 9.6 12.455 c_max = 0.2 dL[-1] = 0.
# chord = np.array([0.11974937, 0.1988364, 0.19845687, 0.18358676, 0.19946962, 0.18837444, 0.15840381, 0.149858,
#                   0.10195713, 0.00286152])
# twist = np.array([0.36476035, 0.1984514, 0.31483518, 0.23364511, 0.25401822, 0.18120287, 0.212223, 0.19510032,
#                   0.15543919, 0.07810287])
# omega = 5704.80803182 * 2*np.pi/60
# case30 = (omega, chord, twist, 'case30')
#
# # 300 p 500g 9.6 12.455 c_max 0.6 UNIFORM INFLOW
# chord = np.array([0.10391889, 0.18799523, 0.27614229, 0.37562084, 0.45449965, 0.48123414, 0.58121353, 0.58874668,
#                   0.56506879, 0.47377971])
# twist = np.array([0.04186956, 0.19010911, 0.21211266, 0.36218881, 0.36069853, 0.2886057, 0.23230848, 0.26068052,
#                   0.23741129, 0.13761006])
# omega = 3421.57237585 * 2*np.pi/60
# case31 = (omega, chord, twist, 'case31')
#
# # 300p 900g 9.6 12.455 cmax = 0.6 dL[-1] = 0 in both ff and hover
# chord = np.array([0.11970866, 0.21964796, 0.31853108, 0.41463554, 0.48939962, 0.40686682, 0.40969556, 0.31039211,
#                   0.23005764, 0.14189208])
# twist = np.array([0.14411419, 0.31856365, 0.3965076, 0.33163697, 0.30031886, 0.28954363, 0.25050936, 0.24496484,
#                   0.20212016, 0.04741441])
# omega = 3808.61192073 * 2*np.pi/60
# case32 = (omega, chord, twist, 'case32')
#
# # 300p 900g 9.6 12.455 cmax = 0.4 dL[-1] = 0 in both ff and hover
# chord = np.array([0.11898545, 0.21666921, 0.31555315, 0.39153882, 0.39930649, 0.39440524, 0.34675179, 0.31040638,
#                   0.22027865, 0.13270832])
# twist = np.array([0.3537641, 0.44071777, 0.31258455, 0.34557684, 0.29172076, 0.26002447, 0.25374399, 0.2129489,
#                   0.19001869, 0.01672236])
# omega = 4033.79848792 * 2*np.pi/60
# case33 = (omega, chord, twist, 'case33')
#
# # 300p 900g 9.6 12.455 cmax = 0.2 dL[-1] = 0 in both ff and hover
# chord = np.array([0.09063594, 0.14664063, 0.19840765, 0.19864742, 0.19942876, 0.19605742, 0.15698326, 0.14076759,
#                   0.10133328, 0.00135014])
# twist = np.array([0.44599459, 0.33372559, 0.27601083, 0.26067019, 0.22883813, 0.19942305, 0.21523722, 0.18224397,
#                   0.15517816, 0.10341055])
# omega = 5672.00840106 * 2*np.pi/60
# case34 = (omega, chord, twist, 'case34')
#
# # case 32 rerun
# chord = np.array([0.11970866, 0.21964796, 0.31853108, 0.41463554, 0.48939962, 0.40686682, 0.40969556, 0.31039211,
#                   0.23005764, 0.14189208])
# twist = np.array([0.14411419, 0.31856365, 0.3965076, 0.33163697, 0.30031886, 0.28954363, 0.25050936, 0.24496484,
#                   0.20212016, 0.04741441])
# omega = 3808.61192073 * 2*np.pi/60
# case32b = (omega, chord, twist, 'case32b')
#
# # 300p 1100 g 9.6 12.455 cmax = 0.6 dL[-1] = 0 in ff and hover
# chord = np.array([0.11973573, 0.21967441, 0.31809143, 0.41691464, 0.49708984, 0.41769755, 0.40226111, 0.30265944,
#                   0.23442085, 0.14098098])
# twist = np.array([0.14411027, 0.3185936, 0.39672945, 0.33185516, 0.29927734, 0.28975136, 0.2545207, 0.242937,
#                   0.20009231, 0.0445423])
# omega = 3808.60966422 * 2*np.pi/60
# case35 = (omega, chord, twist, 'case35')

# chord = np.array([0.11302328, 0.21298767, 0.31283867, 0.36367498, 0.28204006, 0.3687244, 0.3455101, 0.25136865,
#                   0.23114779, 0.13553282])
# twist = np.array([0.20810447, 0.35332739, 0.42475926, 0.33883661, 0.27949311, 0.31701892, 0.28471683, 0.25705176,
#                   0.18672455, 0.0562048 ])
# omega = 3903.31609145 * 2*np.pi/60
# hover3 = (omega, chord, twist, 'hover3')
#
# # 300 pop 500 gen c_max = 0.6 9.6 12.455 times = 500 s hover, 1 ff
# chord = np.array([0.11590125, 0.21392332, 0.31352763, 0.41232689,  0.44686989, 0.5297269, 0.53100342, 0.46727961,
#                   0.5509099, 0.49032908])
# twist = np.array([0.22640486, 0.35004697, 0.40536963, 0.39863093,  0.34831122, 0.33632961, 0.31779411, 0.2370445,
#                   0.19725237, 0.19587309])
# omega = 3318.1009846 * 2*np.pi/60
# hover4 = (omega, chord, twist, 'hover4')
#
# # 300 pop 900 gen cmax = 0.6 9.6 12.455 times = 500 s hover, 1 s ff dL[-1] = 0 in ff and hover
# chord = np.array([0.11381568, 0.21302089, 0.30221703, 0.39052546, 0.42872102, 0.37317609, 0.31394942, 0.24854337,
#                   0.19700262, 0.12940342])
# twist = np.array([0.22708904, 0.37496226, 0.40598074, 0.37456734, 0.31760063, 0.3044958, 0.25391602, 0.24288921,
#                   0.19809018, 0.24893911])
# omega = 3905.72433364 * 2*np.pi/60
# hover5 = (omega, chord, twist, 'hover5')
#
# # 300 p 900 g cmax = 0.6 9.6 12.455 times = 500 s hover 0.001 s ff dL[-1] = 0 in ff and hover
# chord = np.array([ 0.11608058,  0.21604469,  0.29597747,  0.38671991,  0.42487698, 0.36535316,  0.31896227,  0.25117344,  0.17412246,  0.1326255 ])
# twist = np.array([ 0.23066542,  0.33838942,  0.40105013,  0.37432279,  0.33179812, 0.30570999,  0.27443754,  0.25515191,  0.22798693,  0.2482371 ])
# omega = 3828.25590818 * 2*np.pi/60
# hover6 = (omega, chord, twist, 'hover6')