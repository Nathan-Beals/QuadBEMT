import lookup_table
import numpy as np
import matplotlib.pyplot as plt


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


def get_Cl_fun(Re, Cl_table):
    Clalpha, Cl0 = get_Clalpha_Cl0(Re, Cl_table)
    print Clalpha
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
    for i in xrange(len(Cd_table[0][1])):
        if isclose(Cd_table[0][1][i], Re) and not np.isnan(Cd_table[1][i]):
            alphas.append(Cd_table[0][0][i])
            Cds.append(Cd_table[1][i])

    i_upper = [i for i in xrange(len(alphas)) if 5.0 <= alphas[i]]
    alphas_upper = np.array([alphas[i] for i in i_upper])
    Cds_upper = np.array([Cds[i] for i in i_upper])
    z_upper = np.polyfit(alphas_upper*2*np.pi/360, Cds_upper, 6)

    i_lower = [i for i in xrange(len(alphas)) if alphas[i] < 5.0]
    alphas_lower = np.array([alphas[i] for i in i_lower])
    Cds_lower = np.array([Cds[i] for i in i_lower])
    z_lower = np.polyfit(alphas_lower*2*np.pi/360, Cds_lower, 4)

    Cd_fun_upper = np.poly1d(z_upper)
    Cd_fun_lower = np.poly1d(z_lower)

    def Cd_fun(alpha):
        alpha_deg = alpha * 360 / 2 / np.pi
        Cd = []
        for i, a in enumerate(alpha_deg):
            if a < 5.0:
                Cd.append(Cd_fun_lower(alpha[i]))
            else:
                Cd.append(Cd_fun_upper(alpha[i]))
        Cd = np.array(Cd)
        return Cd
    return Cd_fun


def closest_Re(Re, Res_in_table):
    return np.array([min(Res_in_table, key=lambda x: abs(x-Re))])


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def main():
    airfoil = 'SDA1075_494p'
    Cl_table, Cd_table = create_Cl_Cd_table(airfoil)
    target_Re = 100000.
    Cl_fun = get_Cl_fun(target_Re, Cl_table)
    Cd_fun = get_Cd_fun(target_Re, Cd_table)

    alphas = np.linspace(-10., 20., 40)
    alphas_rad = alphas * np.pi / 180

    # Compute simplified function values
    Cls = Cl_fun(alphas_rad)
    Cds = Cd_fun(alphas_rad)

    # Get values from table
    indices = [i for i in xrange(len(Cl_table[0][1])) if isclose(Cl_table[0][1][i], target_Re)]
    alphas_table = np.array([Cl_table[0][0][i] for i in indices])
    Cls_table = np.array([Cl_table[1][i] for i in indices])
    Cds_table = np.array([Cd_table[1][i] for i in indices])

    plt.figure(1)
    plt.plot(alphas, Cls, alphas_table, Cls_table)
    plt.xlabel('angle of attack')
    plt.ylabel('Cl')
    #plt.legend(['Function', 'Table'])

    plt.figure(2)
    plt.plot(alphas, Cds, alphas_table, Cds_table)
    plt.xlabel('angle of attack')
    plt.ylabel('Cd')
    #plt.legend(['Function', 'Table'])

    plt.show()


if __name__ == "__main__":
    main()
