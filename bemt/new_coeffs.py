import lookup_table
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import itertools


def create_Cl_Cd_table(airfoil):
    alphas, Re, CL, CD, Clmax = lookup_table.create_table(airfoil)
    Cl_table = ((alphas, Re), CL)
    Cd_table = ((alphas, Re), CD)
    return Cl_table, Cd_table, Clmax


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


def get_Cl_fun(Re, Cl_table, Clmax):
    Clalpha, Cl0, alpha0 = get_liftCurveInfo(Re, Cl_table)

    def Cl_fun(alpha):
        """
        :param alpha: Angle of attack in radians at which to obtain the lift coefficient at a given Reynolds number
        :return: Lift coefficient
        """
        Cl = Cl0 + Clalpha*alpha
        Cl[Cl > Clmax] = Clmax
        return Cl

    return Cl_fun


def get_Cd_fun(Re, Cd_table):
    alphas = []
    Cds = []
    print max(Cd_table[1])
    for i in xrange(len(Cd_table[0][1])):
        if isclose(Cd_table[0][1][i], Re) and not np.isnan(Cd_table[1][i]):
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


def closest_Re(Re, Res_in_table):
    if type(Re) is not np.ndarray:
        try:
            Re = np.array([float(Re)])
        except TypeError:
            Re = np.array(Re)
    return np.array([min(Res_in_table, key=lambda x: abs(x-R)) for R in Re])


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def main():
    airfoil = 'SDA1075_494p'
    Cl_table, Cd_table, Clmax = create_Cl_Cd_table(airfoil)
    alphas = np.linspace(-5., 25., 100)
    alphas_rad = np.linspace(-5., 25., 100) * np.pi / 180
    target_Re = [250000., 100000., 60000.,  40000., 20000.]
    markers = ['kD-', 'ko-', 'ks-', 'kv-', 'k*-']

    results = {}
    for Re in target_Re:
        Cl_fun = get_Cl_fun(Re, Cl_table, Clmax[Re])
        Cd_fun = get_Cd_fun(Re, Cd_table)
        Cls = Cl_fun(alphas_rad)
        Cds = Cd_fun(alphas_rad)
        results[Re] = (Cls, Cds)

    # #Get values from table
    # Cl_dict = dict(zip(zip(Cl_table[0][0], Cl_table[0][1]), Cl_table[1]))
    # Cd_dict = dict(zip(zip(Cd_table[0][0], Cd_table[0][1]), Cd_table[1]))
    # indices = [i for i in xrange(len(Cl_table[0][1])) if isclose(Cl_table[0][1][i], target_Re)]
    # alphas_table = np.array([Cl_table[0][0][i] for i in indices])
    # Cls_table = np.array([Cl_table[1][i] for i in indices])
    # Cds_table = np.array([Cd_table[1][i] for i in indices])

    plt.figure(1)
    for Re, marker in zip(reversed(target_Re), reversed(markers)):
        plt.plot(alphas, results[Re][0], marker, markevery=6, markerfacecolor='white')
    plt.xlim([-5., 25.])
    plt.ylim([-1.0, 2.0])
    plt.xlabel(r'$\alpha,\,\mathrm{deg}$', fontsize=18)
    plt.ylabel(r'$\mathrm{C}_\mathrm{l}$', fontsize=18)
    plt.legend(['Re=20k', 'Re=40k', 'Re=60k', 'Re=100k', 'Re=250k'], loc='upper left')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.grid()

    plt.figure(2)
    for Re, marker in zip(reversed(target_Re), reversed(markers)):
        plt.plot(alphas, results[Re][1], marker, markevery=6, markerfacecolor='white')
    plt.xlim([-5., 25.])
    plt.ylim([0.0, 0.35])
    plt.xlabel(r'$\alpha,\,\mathrm{deg}$', fontsize=18)
    plt.ylabel(r'$\mathrm{C}_\mathrm{d}$', fontsize=18)
    plt.legend(['Re=20k', 'Re=40k', 'Re=60k', 'Re=100k', 'Re=250k'], loc='upper left')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.grid()

    plt.show()


if __name__ == "__main__":
    main()
