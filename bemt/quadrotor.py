class Quadrotor(object):
    def __init__(self, propeller, weight):
        self.propeller = propeller    # List of propellers
        # Quadrotor weight in Newtons. Got this figure by comparing weights of popular vehicles with approximately
        # 9 inch diameter propellers
        self.weight = weight

    def frame_drag(self, alpha, v_inf, dens):
        q = 0.5 * dens * v_inf**2
        drag = (-0.0045*alpha + 0.22)*q     # Quadrotor frame drag in units of lb
        return drag
