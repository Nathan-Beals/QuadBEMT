import numpy as np
from lookup_table import create_table, interpolate
from scipy.interpolate import griddata
import unit_conversion
import resource


class Propeller(object):
    def __init__(self, twist, chord, radius, n_blades, r, y, dr, dy, solidity=None, airfoils=None, Cl_tables={},
                 Cd_tables={}):
        """
        :param twist: Twist distribution, numpy array
        :param chord: Chord distribution, numpy array
        :param radius: blade radius, float
        :param n_blades: number of blades, integer
        :param r: non-dimensional radial element locations
        :param y: dimensional radial element locations
        :param dr: non-dimensional distance between elements
        :param dy: dimensional distance between elements
        :param solidity: optional argument for use in test cases where solidity is constant along the blade
        :param airfoils: in the form (('airfoil1', r_start, r_end), ('airfoil2', r_start, r_end), etc)
        :return:
        """
        if solidity is None:
            self.chord = chord
            self.solidity = n_blades*chord/np.pi/radius
        else:
            # Handle case for when we just want to pass a solidity
            self.solidity = solidity
            self.chord = solidity*np.pi*radius/n_blades

        if airfoils is None:
            self.airfoils = (('simple', 0., 1.),)
        elif sorted(airfoils) == sorted((('simple', 0.0, 1.0),)):
            self.airfoils = airfoils
        else:
            self.airfoils = airfoils
            self.Cl_tables = Cl_tables
            self.Cd_tables = Cd_tables
            if not Cl_tables:
                for airfoil in self.airfoils:
                    # Create the lookup table. The alpha values returned from create_table are in degrees, therefore the
                    # table will need to be queried in degrees
                    alpha, Re, CL, CD = create_table(airfoil[0])

                    self.Cl_tables[airfoil[0]] = ((alpha, Re), CL)
                    self.Cd_tables[airfoil[0]] = ((alpha, Re), CD)

        self.twist = twist
        self.radius = radius
        self.n_blades = n_blades
        self.r = r
        self.y = y
        self.dr = dr
        self.dy = dy

        # airfoil_fun_r is a list of airfoil names where airfoil_fun_r[0] is the airfoil section at r = 0 and
        # airfoil_fun_r[-1] is the airfoil section name at r = 1
        self.airfoil_fun_r = []
        for r_loc in self.r:
            self.airfoil_fun_r.append(self.get_airfoil(r_loc))

    def get_airfoil(self, r_loc):
        for airfoil in self.airfoils:
            if airfoil[1] <= r_loc <= airfoil[2]:
                return airfoil[0]
        return 'SDA1075_494p'

    def get_Clalpha_alpha0(self, Re):
        """
        Calculate the values of lift curve slope for the airfoil sections along the span of the blade
        according to a spanwise distribution of Reynolds numbers. This works only for constant airfoil section blades.
        :param Re: Reynolds number along the span of the blade.
        :return: Lift curve slopes for the airfoil sections along the span of the blade.
        """
        airfoil = self.airfoils[0][0]
        if airfoil == 'simple':
            return np.ones(len(Re))*2*np.pi, np.zeros(len(Re))
        table = self.Cl_tables[airfoil]
        alfas = [0.0]*len(Re) + [5.0]*len(Re)
        points2interp = zip(alfas, np.concatenate((Re, Re)))
        points = griddata(table[0], table[1], points2interp, method='nearest')
        Cl0 = points[:len(points)/2]
        Cl5deg = points[len(points)/2:]
        Clalpha = (Cl5deg - Cl0)/(5*2*np.pi/360)
        alpha0 = -Cl0 / Clalpha
        return Clalpha, alpha0

    def get_Cl(self, aoa, Re, method='regular'):
        """
        :param aoa: Angle of attack in radians, numpy array along the span
        :param Re: Reynolds number, numpy array along the span
        :return: Cl values along the span
        """
        aoa_deg = aoa * 360 / 2 / np.pi
        if self.airfoils[0][0] == 'simple':
            return 2*np.pi*aoa
        else:
            airfoil = self.airfoils[0][0]
            table = self.Cl_tables[airfoil]
            Cl = griddata(table[0], table[1], zip(aoa_deg, Re), method='nearest')
            return Cl

    def get_Cd(self, aoa, Re):
        """
        :param aoa: Angle of attack in radians, numpy array along the span
        :param Re: Reynolds number along the span
        :return: Cd along the span
        """
        aoa_deg = aoa * 360 / 2 / np.pi
        if self.airfoils[0][0] == 'simple':
            return 0.02 - 0.0216*aoa + 0.400*aoa**2
        else:
            airfoil = self.airfoils[0][0]
            table = self.Cd_tables[airfoil]
            Cd = griddata(table[0], table[1], zip(aoa_deg, Re), method='nearest')
            return Cd

    ############################################################################################################
    # Methods for limited Reynolds number availability.
    ############################################################################################################

    def get_Clalpha_Cl0(self, Re):
        """
        Get lift curve slope and lift coefficient at zero angle of attack for a single Reynolds number.
        :param Re: Reynolds number at which to get Clalpha and Cl0
        :return: Lift curve slope and lift coefficient at zero angle of attack
        """
        Cl_table = self.Cl_tables[self.airfoils[0][0]]
        Cl_dict = dict(zip(zip(Cl_table[0][0], Cl_table[0][1]), Cl_table[1]))

        # If Re is not a numpy array, convert it to one
        if type(Re) is not np.ndarray:
            try:
                Re = np.array([float(Re)])
            except TypeError:
                Re = np.array(Re)

        # Find the nearest table Reynolds number(s) to the input Reynolds number(s)
        nearest_Re = self.closest_Re(Re, set(Cl_table[0][1]))

        points = np.concatenate((zip((0.0,)*len(Re), nearest_Re), zip((5.0,)*len(Re), nearest_Re)))
        vals = np.array([Cl_dict[tuple(point)] for point in points])
        Cl0 = vals[:len(points)/2]
        Cl5deg = vals[len(points)/2:]
        Clalpha = (Cl5deg - Cl0) / (5*2*np.pi/360)
        alpha0 = -Cl0 / Clalpha
        return Clalpha, Cl0, alpha0

    def get_Cl_fun(self, Re):
        Clalpha, Cl0, alpha0 = self.get_Clalpha_Cl0(Re)

        def Cl_fun(alpha):
            """
            :param alpha: Angle of attack in radians at which to obtain the lift coefficient at a given Reynolds number
            :return: Lift coefficient
            """
            Cl = Cl0 + Clalpha*alpha
            return Cl
        return Cl_fun

    # def get_Cd_fun(self, Re):
    #     Cd_table = self.Cd_tables[self.airfoils[0][0]]
    #     alphas = []
    #     Cds = []
    #     for i in xrange(len(Cd_table[0][1])):
    #         if self.isclose(Cd_table[0][1][i], Re) and not np.isnan(Cd_table[1][i]):
    #             alphas.append(Cd_table[0][0][i])
    #             Cds.append(Cd_table[1][i])
    #     alphas = np.array(alphas)
    #     Cds = np.array(Cds)
    #     z = np.polyfit(alphas*2*np.pi/360, Cds, 6)
    #     Cd_fun = np.poly1d(z)
    #     return Cd_fun

    def get_Cd_fun(self, Re):
        Cd_table = self.Cd_tables[self.airfoils[0][0]]
        alphas = []
        Cds = []
        for i in xrange(len(Cd_table[0][1])):
            if self.isclose(Cd_table[0][1][i], Re) and not np.isnan(Cd_table[1][i]):
                alphas.append(Cd_table[0][0][i])
                Cds.append(Cd_table[1][i])

        i_upper = [i for i in xrange(len(alphas)) if 5.0 < alphas[i]]
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
            #return 0.02 - 0.0216*alpha + 0.400*alpha**2
            if alpha_deg < 5.0:
                return Cd_fun_lower(alpha)
            else:
                return Cd_fun_upper(alpha)
            # Cd = []
            # for ix, a in enumerate(alpha_deg):
            #     if a < 5.0:
            #         Cd.append(Cd_fun_lower(alpha[ix]))
            #     else:
            #         Cd.append(Cd_fun_upper(alpha[ix]))
            # Cd = np.array(Cd)
            # return Cd
        return Cd_fun

    @staticmethod
    def closest_Re(Re, Res_in_table):
        if type(Re) is not np.ndarray:
            try:
                Re = np.array([float(Re)])
            except TypeError:
                Re = np.array(Re)
        return np.array([min(Res_in_table, key=lambda x: abs(x-R)) for R in Re])

    @staticmethod
    def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)