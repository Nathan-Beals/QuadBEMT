from unit_conversion import lb2N
import numpy as np


class Quadrotor(object):
    def __init__(self, propeller, weight):
        self.propeller = propeller    # List of propellers
        # Quadrotor weight in Newtons. Got this figure by comparing weights of popular vehicles with approximately
        # 9 inch diameter propellers
        self.weight = weight

    def frame_drag(self, alpha, v_inf, dens):
        q = 0.5 * dens * v_inf**2
        # alpha_deg = alpha * 180 / np.pi
        # drag = lb2N((-0.0045*alpha_deg + 0.22)*q)     # Quadrotor frame drag in units of lb converted to Newtons
        # print "frame drag = " + str(drag)
        simple_drag = 0.27 * 0.092903 * q
        return simple_drag
