import numpy as np
import lookup_table
from scipy.interpolate import griddata
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import itertools


def create_Cl_Cd_table(airfoil):
    """
    Create "tables" (really these are tuples of the form ((angle of attack, Reynolds number), Cl)) of lift and drag
    coefficient as a function of angle of attack and Reynolds number.
    :param airfoil: The name of the airfoil. This name should have a corresponding folder in BEMT/bemt/airfoil_data with
    output from XFOIL.
    :return: The Cl and Cd "tables" as described above.
    """
    alphas, Re, CL, CD, Clmax = lookup_table.create_table(airfoil)
    Cl_table = ((alphas, Re), CL)
    Cd_table = ((alphas, Re), CD)
    return Cl_table, Cd_table, Clmax


def get_Cl(aoa, Re, table):
    """
    :param aoa: Angle of attack in radians, numpy array along the span
    :param Re: Reynolds number, numpy array along the span
    :param table: Cl_table for a specific airfoil in the form ((angle of attack, Reynolds number), Cl). If the table is
    empty (i.e., 'simple' airfoil) use Cl = 2*pi*alpha
    :return: Cl values along the span
    """
    aoa_deg = aoa * 360 / 2 / np.pi
    if not table:
        return 2*np.pi*aoa
    else:
        print "griddata about to be called"
        raise Exception
        Cl = griddata(table[0], table[1], zip(aoa_deg, Re), method='nearest')
        return Cl


def get_Cd(aoa, Re, table):
    """
    :param aoa: Angle of attack in radians, numpy array along the span
    :param Re: Reynolds number along the span
    :param table: Cd_table for a specific airfoil in the form ((angle of attack, Reynolds number), Cd). If the table is
    empty (i.e., 'simple' airfoil) use quadratic Cd.
    :return: Cd along the span
    """
    aoa_deg = aoa * 360 / 2 / np.pi
    if not table:
        return 0.02 - 0.0216*aoa + 0.400*aoa**2
    else:
        "griddata about to be called"
        Cd = griddata(table[0], table[1], zip(aoa_deg, Re), method='nearest')
        return Cd


def get_liftCurveInfo(Re, table):
    """
    Calculate the lift curve slope, lift coefficient at zero angle of attack, and zero lift angle of attack along the
    span of the blade using the table provided (corresponding to a specific airfoil). Finds the point in the table
    nearest to each Reynolds number in the Re input array and returns the values corresponding to that point in the
    table (i.e., there is no interpolation).
    :param Re: Reynolds numbers along the span of the blade. Can  be a numpy array (entire span) or float.
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
    """

    :param Re: Reynolds number. Should be a float.
    :param Cl_table: Table of the form ((alphas, Reynolds numbers), CLs)
    :param Clmax: Dictionary of the form {Reynolds number: Clmax} for all Reynolds numbers in Cl_table
    :return:
    """
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