import numpy as np
from lookup_table import create_table
from scipy.interpolate import griddata
from unit_conversion import rad2deg


class Propeller(object):
    def __init__(self, twist, chord, radius, n_blades, r, y, dr, dy, Clalpha, solidity=None, airfoils=None):
        """
        :param twist: Twist distribution, numpy array
        :param chord: Chord distribution, numpy array
        :param radius: blade radius, float
        :param n_blades: number of blades, integer
        :param r: non-dimensional radial element locations
        :param y: dimensional radial element locations
        :param dr: non-dimensional distance between elements
        :param dy: dimensional distance between elements
        :param Clalpha: linear lift curve slope, usually 2*pi
        :param solidity: optional argument for use in test cases where solidity is constant along the blade
        :param airfoils: in the form (('airfoil1', r_start, r_end), ('airfoil2', r_start, r_end), etc)
        :return:
        """
        if solidity is None:
            self.chord = chord
            self.solidity = n_blades*chord/(2 * np.pi * y)
        else:
            # Handle case for when we just want to pass a solidity
            self.solidity = solidity
            self.chord = solidity*np.pi*radius/n_blades

        if airfoils is None:
            self.airfoils = ('simple', 0, 1)
        else:
            self.airfoils = airfoils
            print self.airfoils
            self.Cl_tables = {}
            self.Cd_tables = {}
            for airfoil in self.airfoils:
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
        self.Clalpha = Clalpha

        # airfoil_fun_r is a list of airfoil names where airfoil_fun_r[0] is the airfoil section at r = 0 and
        # airfoil_fun_r[-1] is the airfoil section name at r = 1
        self.airfoil_fun_r = []
        for r_loc in self.r:
            self.airfoil_fun_r.append(self.get_airfoil(r_loc))

    def get_airfoil(self, r_loc):
        for airfoil in self.airfoils:
            if airfoil[1] <= r_loc <= airfoil[2]:
                return airfoil[0]
        return 'simple'

    def get_Clalpha(self):
        return np.ones(len(self.r)) * self.Clalpha

    def get_Cl(self, aoa, Re):
        aoa = aoa * 360 / 2 / np.pi
        if self.airfoils == ('simple', 0, 1):
            return self.Clalpha * aoa
        else:
            Cl = np.empty([len(self.r)])
            i = 0
            for airfoil in self.airfoil_fun_r:
                if Re[i] > 10000:
                    Cl[i] = griddata(self.Cl_tables[airfoil][0], self.Cl_tables[airfoil][1], (aoa[i], Re[i]))
                else:
                    Cl[i] = griddata(self.Cl_tables[airfoil][0], self.Cl_tables[airfoil][1], (aoa[i], 10000))
                i += 1
            return Cl

    def get_Cd(self, aoa, Re):
        aoa = aoa * 360 / 2 / np.pi
        if self.airfoils == ('simple', 0, 1):
            return 0.02 - 0.0216*aoa + 0.400*aoa**2
        else:
            Cd = np.empty([len(self.r)])
            i = 0
            for airfoil in self.airfoil_fun_r:
                if Re[i] > 10000:
                    Cd[i] = griddata(self.Cd_tables[airfoil][0], self.Cd_tables[airfoil][1], (aoa[i], Re[i]))
                else:
                    Cd[i] = griddata(self.Cd_tables[airfoil][0], self.Cd_tables[airfoil][1], (aoa[i], 10000))
                i += 1
            return Cd
