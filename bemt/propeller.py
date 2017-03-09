import numpy as np
from lookup_table import create_table


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
        # Initialize propeller chord and solidity
        if solidity is None:
            self.chord = chord
            self.solidity = n_blades*chord/np.pi/radius
        else:
            # Handle case for when we just want to pass a solidity
            self.solidity = solidity
            self.chord = solidity*np.pi*radius/n_blades

        # Initialize propeller Cl and Cd tables as copies as the argument values. This can either be an empty dict or
        # a dictionary with keys corresponding to the airfoil that makes up the propeller.
        self.Cl_tables = Cl_tables
        self.Cd_tables = Cd_tables

        # Handle case where user does not specify airfoils and use simple case.
        if airfoils is None:
            self.airfoils = (('simple', 0., 1.),)
            self.Cl_tables[self.airfoils[0][0]] = ()
            self.Cd_tables[self.airfoils[0][0]] = ()

        # Handle case where the airfoil is specified as simple (Cl = 2*pi*alpha, Cd = quadratic wrt alpha)
        elif sorted(airfoils) == sorted((('simple', 0.0, 1.0),)):
            self.airfoils = airfoils
            self.Cl_tables[self.airfoils[0][0]] = ()
            self.Cd_tables[self.airfoils[0][0]] = ()

        # Handle case where the airfoil is something other than simple. For these cases we assume either the Cl and Cd
        # tables have been created outside of the optimization (i.e., the Cl_tables and Cd_tables arguments are not
        # empty) in which case the attribute variables have already been initialized to the input variables. If the
        # input variables are empty, create the tables from the airfoil data in BEMT/bemt/airfoil_data.
        else:
            self.airfoils = airfoils
            if not Cl_tables:
                for airfoil in self.airfoils:
                    # The alpha values returned from create_table are in degrees, therefore the table will need to be
                    # queried in degrees
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