import lookup_table
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d


def create_Cl_Cd_table(airfoil):
    [alphas, Re, CL, CD] = lookup_table.create_table(airfoil)
    Cl_table = ((alphas, Re), CL)
    Cd_table = ((alphas, Re), CD)
    return Cl_table, Cd_table


def get_Clalpha_Cl0(Re, Cl_table):
    """
    Get lift curve slope and lift coefficient at zero angle of attack for a single Reynolds number.
    :param Re: Reynolds number at which to get Clalpha and Cl0
    :param Cl_table: Table in the form of ((alphas, Res), CLs)
    :return: Lift curve slope and lift coefficient at zero angle of attack
    """
    Cl_dict = dict(zip(zip(Cl_table[0][0], Cl_table[0][1]), Cl_table[1]))

    # If Re is not a numpy array, convert it to one
    if type(Re) is not np.ndarray:
        try:
            Re = np.array([float(Re)])
        except TypeError:
            Re = np.array(Re)

    # Find the nearest table Reynolds number(s) to the input Reynolds number(s)
    nearest_Re = closest_Re(Re, set(Cl_table[0][1]))

    points = np.concatenate((zip((0.0,)*len(Re), nearest_Re), zip((5.0,)*len(Re), nearest_Re)))
    vals = np.array([Cl_dict[tuple(point)] for point in points])
    Cl0 = vals[:len(points)/2]
    Cl5deg = vals[len(points)/2:]
    Clalpha = (Cl5deg - Cl0) / (5*2*np.pi/360)
    return Clalpha, Cl0


def get_liftCurveInfo(Re, table):
    """
    Calculate the lift curve slope, lift coefficient at zero angle of attack, and zero lift angle of attack along the
    span of the blade using the table provided (corresponding to a specific airfoil). Finds the point in the table
    nearest to each Reynolds number in the Re input array and returns the values corresponding to that point in the
    table (i.e., there is no interpolation).
    :param Re: Reynolds numbers along the span of the blade.
    :param table: Cl table in the form ((angle of attack, Reynolds number), Cl) corresponding to a specific airfoil.
    If a 'simple' airfoil is used this will be an empty tuple and trigger the first 'if' statement.
    :return: Clalpha, Cl0, and alpha0 along the span of the blade.
    """

    if not table:
        return np.ones(len(Re))*2*np.pi, np.zeros(len(Re)), np.zeros(len(Re))

    Cl_dict = dict(zip(zip(table[0][0], table[0][1]), table[1]))
    nearest_Re = closest_Re(Re, set(table[0][1]))
    alfas = [0.0]*len(nearest_Re) + [5.0]*len(nearest_Re)
    points = zip(alfas, np.concatenate((nearest_Re, nearest_Re)))
    vals = np.array([Cl_dict[tuple(point)] for point in points])
    Cl0 = vals[:len(points)/2]
    Cl5deg = vals[len(points)/2:]
    Clalpha = (Cl5deg - Cl0)/(5*2*np.pi/360)
    alpha0 = -Cl0 / Clalpha
    return Clalpha, Cl0, alpha0


def get_Cl_fun(Re, Cl_table):
    Clalpha, Cl0, alpha0 = get_liftCurveInfo(Re, Cl_table)

    def Cl_fun(alpha):
        """
        :param alpha: Angle of attack in radians at which to obtain the lift coefficient at a given Reynolds number
        :return: Lift coefficient
        """
        Cl = Cl0 + Clalpha*alpha
        return Cl

    return Cl_fun


def get_Cd_fun(Re, Cd_table):
    alphas = []
    Cds = []
    max_aoa = max(Cd_table[0][0])
    print max(Cd_table[1])
    for i in xrange(len(Cd_table[0][1])):
        if isclose(Cd_table[0][1][i], Re) and -10.0 <= Cd_table[0][0][i] and not np.isnan(Cd_table[1][i]) and -100 < Cd_table[1][i] < 100:
            alphas.append(Cd_table[0][0][i])
            Cds.append(Cd_table[1][i])
    alphas = np.array(alphas)
    Cds = np.array(Cds)
    min_Cd = min(Cds)

    z = np.polyfit(alphas*2*np.pi/360, Cds, 2)
    Cd_poly = np.poly1d(z)

    def Cd_fun(alpha):
        Cd = Cd_poly(alpha)
        try:
            if Cd < min_Cd:
                Cd = min_Cd
        except ValueError:
            low_indices = Cd < min_Cd
            Cd[low_indices] = min_Cd
        return Cd

    return Cd_fun


def interp2d_fun(table):
    alpha = table[0][0]
    Re = table[0][1]
    coeff = table[1]
    print "len(alpha) = " + str(len(alpha))
    print "len(Re) = " + str(len(Re))
    print "len"
    return interp2d(Re, alpha, coeff)


def closest_Re(Re, Res_in_table):
    return np.array([min(Res_in_table, key=lambda x: abs(x-Re))])


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def main():
    airfoil = 'SDA1075_494p'
    Cl_table, Cd_table = create_Cl_Cd_table(airfoil)
    alphas = np.linspace(-10., 40., 60)
    alphas_rad = np.linspace(-10., 40., 60) * np.pi / 180
    #target_Re = [1000000., 500000., 250000., 100000., 90000., 80000., 70000., 60000., 50000., 40000., 30000., 20000., 10000.]
    target_Re = 100000.
    # plt.figure(1)
    # for Re in target_Re:
    #     Cl_fun = get_Cl_fun(Re, Cl_table)
    #     Cd_fun = get_Cd_fun(Re, Cd_table)
    #     Cls = Cl_fun(alphas_rad)
    #     Cds = Cd_fun(alphas_rad)
    #     plt.plot(alphas, Cds)
    #     for i in xrange(len(alphas)):
    #         if Cds[i] < 0.0:
    #             print "Cd = %s at Re = %s and alpha = %s" % (str(Cds[i]), str(Re), alphas[i])
    # plt.show()

    #Get values from table
    indices = [i for i in xrange(len(Cl_table[0][1])) if isclose(Cl_table[0][1][i], target_Re)]
    alphas_table = np.array([Cl_table[0][0][i] for i in indices])
    Cls_table = np.array([Cl_table[1][i] for i in indices])
    Cds_table = np.array([Cd_table[1][i] for i in indices])

    Cl_fun = get_Cl_fun(target_Re, Cl_table)
    #Cd_fun = get_Cd_fun(target_Re, Cd_table)
    Cd_fun = interp2d_fun(Cd_table)
    Cls = Cl_fun(alphas_rad)
    #Cds = Cd_fun(np.ones(len(alphas_rad))*target_Re, alphas_rad)

    # plt.figure(1)
    # plt.plot(alphas, Cls, alphas_table, Cls_table)
    # plt.xlabel('angle of attack')
    # plt.ylabel('Cl')
    #plt.legend(['Function', 'Table'])

    plt.figure(2)
    plt.plot(alphas, Cds),
    #plt.plot(alphas_table, Cds_table)
    plt.xlabel('angle of attack')
    plt.ylabel('Cd')
    #plt.legend(['Function', 'Table'])

    plt.show()


if __name__ == "__main__":
    main()
