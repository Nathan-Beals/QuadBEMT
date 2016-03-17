import numpy as np


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
            airfoil_data = []
            for airfoil in self.airfoils:
                airfoil_data.append()

        self.twist = twist
        self.radius = radius
        self.n_blades = n_blades
        self.r = r
        self.y = y
        self.dr = dr
        self.dy = dy
        self.Clalpha = Clalpha

    def get_Clalpha(self):
        if self.airfoils[0] == 'simple'
            return np.ones(len(self.r)) * self.Clalpha
        else

    def get_Cl(self, aoa, r):
        return self.Clalpha * aoa

    def get_Cd(self, aoa, r):
        return 0.02 - 0.0216*aoa + 0.400*aoa**2

    def get_airfoil(self, r):
