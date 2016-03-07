import numpy as np


class Propeller(object):
    def __init__(self, twist, chord, radius, n_blades, r, y, dr, dy, Clalpha, solidity=None):
        if solidity is None:
            self.chord = chord
            self.solidity = n_blades*chord/(np.pi * radius)
        else:
            # Handle case for when we just want to pass a solidity
            self.solidity = solidity
            self.chord = solidity*np.pi*radius/n_blades
        self.twist = twist
        self.radius = radius
        self.n_blades = n_blades
        self.r = r
        self.y = y
        self.dr = dr
        self.dy = dy
        self.Clalpha = Clalpha

    def get_Clalpha(self):
        return np.ones(len(self.r)) * self.Clalpha

    def get_Cl(self, aoa):
        return self.Clalpha * aoa

    def get_Cd(self, aoa):
        return 0.01 + 0.025*aoa + 0.65*aoa**2