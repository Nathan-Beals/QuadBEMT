import numpy as np
from lookup_table import create_table, interpolate
from scipy.interpolate import griddata
from unit_conversion import rad2deg
import resource


class Propeller(object):
    def __init__(self, twist, chord, radius, n_blades, r, y, dr, dy, solidity=None, airfoils=None):
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
            #self.solidity = n_blades*chord/(2 * np.pi * y)
            self.solidity = n_blades*chord/np.pi/radius
            # try:
            #     self.solidity = n_blades*chord/(2 * np.pi * y)
            # except ValueError:
            #     print "chord = " + str(chord)
            #     print "y = " + str(y)
            #     print "chord len = " + str(len(chord))
            #     print "y len = " + str(len(y))
            #     print "r len = " + str(len(r))
            #     raise
        else:
            # Handle case for when we just want to pass a solidity
            self.solidity = solidity
            self.chord = solidity*np.pi*radius/n_blades

        if airfoils is None:
            self.airfoils = (('simple', 0, 1),)
        else:
            self.airfoils = airfoils
            self.Cl_tables = {}
            self.Cd_tables = {}
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

    # def get_Clalpha(self, Re):
    #     """
    #     Calculate the values of lift curve slope for the airfoil sections along the span of the blade
    #     according to a spanwise distribution of Reynolds numbers
    #     :param Re: Reynolds number along the span of the blade.
    #     :return: Lift curve slopes for the airfoil sections along the span of the blade.
    #     """
    #     Clalpha = np.empty([len(self.r)])
    #     for i, airfoil in enumerate(self.airfoil_fun_r):
    #         table = self.Cl_tables[airfoil]
    #         if Re[i] > 10000:
    #             Cl0 = interpolate(table, 0.0, Re[i])
    #             Cl5deg = interpolate(table, 5.0, Re[i])
    #             Clalpha[i] = (Cl5deg - Cl0)/(5*2*np.pi/360)
    #             if np.isnan(Clalpha[i]):
    #                 print "Cl0 = " + str(Cl0)
    #                 print "Cl5deg = " + str(Cl5deg)
    #         else:
    #             Cl0 = interpolate(table, 0.0, 10000)
    #             Cl5deg = interpolate(table, 5.0, 10000)
    #             Clalpha[i] = (Cl5deg - Cl0)/(5*2*np.pi/360)
    #     return Clalpha

    def get_Clalpha_alpha0(self, Re):
        """
        Calculate the values of lift curve slope for the airfoil sections along the span of the blade
        according to a spanwise distribution of Reynolds numbers. This works only for constant airfoil section blades.
        :param Re: Reynolds number along the span of the blade.
        :return: Lift curve slopes for the airfoil sections along the span of the blade.
        """
        airfoil = self.airfoils[0][0]
        table = self.Cl_tables[airfoil]
        alfas = [0.0]*len(Re) + [5.0]*len(Re)
        points2interp = zip(alfas, np.concatenate((Re, Re)))
        points = griddata(table[0], table[1], points2interp, method='nearest')
        Cl0 = points[:len(points)/2]
        Cl5deg = points[len(points)/2:]
        Clalpha = (Cl5deg - Cl0)/(5*2*np.pi/360)
        alpha0 = -Cl0 / Clalpha
        return Clalpha, alpha0

    # def get_Cl(self, aoa, Re):
    #     """
    #     :param aoa: Angle of attack in radians, numpy array along the span
    #     :param Re: Reynolds number, numpy array along the span
    #     :return: Cl values along the span
    #     """
    #     aoa_deg = aoa * 360 / 2 / np.pi
    #     if self.airfoils == ('simple', 0, 1):
    #         return 2*np.pi * aoa
    #     else:
    #         Cl = np.empty([len(self.r)])
    #         i = 0
    #         for airfoil in self.airfoil_fun_r:
    #             table = self.Cl_tables[airfoil]
    #             if Re[i] > 10000:
    #                 Cl[i] = interpolate(table, aoa_deg[i], Re[i])
    #             else:
    #                 Cl[i] = interpolate(table, aoa_deg[i], 10000)
    #             i += 1
    #         return Cl

    def get_Cl(self, aoa, Re):
        """
        :param aoa: Angle of attack in radians, numpy array along the span
        :param Re: Reynolds number, numpy array along the span
        :return: Cl values along the span
        """
        aoa_deg = aoa * 360 / 2 / np.pi
        if self.airfoils == ('simple', 0, 1):
            return 2*np.pi * aoa
        else:
            airfoil = self.airfoils[0][0]
            table = self.Cl_tables[airfoil]
            Cl = griddata(table[0], table[1], zip(aoa_deg, Re), method='nearest')
            return Cl

    # def get_Cd(self, aoa, Re):
    #     """
    #     :param aoa: Angle of attack in radians, numpy array along the span
    #     :param Re: Reynolds number along the span
    #     :return: Cd along the span
    #     """
    #     aoa_deg = aoa * 360 / 2 / np.pi
    #     if self.airfoils == ('simple', 0, 1):
    #         return 0.02 - 0.0216*aoa + 0.400*aoa**2
    #     else:
    #         Cd = np.empty([len(self.r)])
    #         i = 0
    #         for airfoil in self.airfoil_fun_r:
    #             table = self.Cd_tables[airfoil]
    #             if Re[i] > 10000:
    #                 Cd[i] = interpolate(table, aoa_deg[i], Re[i])
    #             else:
    #                 Cd[i] = interpolate(table, aoa_deg[i], 10000)
    #             i += 1
    #         return Cd
    def get_Cd(self, aoa, Re):
        """
        :param aoa: Angle of attack in radians, numpy array along the span
        :param Re: Reynolds number along the span
        :return: Cd along the span
        """
        aoa_deg = aoa * 360 / 2 / np.pi
        if self.airfoils == ('simple', 0, 1):
            return 0.02 - 0.0216*aoa + 0.400*aoa**2
        else:
            airfoil = self.airfoils[0][0]
            table = self.Cd_tables[airfoil]
            Cd = griddata(table[0], table[1], zip(aoa_deg, Re), method='nearest')
            return Cd
