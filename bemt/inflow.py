import numpy as np


def uniform_ff(CT, alpha, mu, n_elements, tip_loss=False):
    with np.errstate(invalid='raise'):
        try:
            inflow = np.sqrt(CT/2)
        except FloatingPointError:
            print "FP error in ff inflow"
            raise
    B = 1.
    converged = False
    max_i = 100
    i = 0
    while not converged and i < max_i:
        inflow_old = inflow
        if tip_loss:
            f = inflow - mu*np.tan(alpha) - CT/(2*B**2*np.sqrt(mu**2+inflow**2))
            f_prime = 1 + CT/2*(mu**2+inflow**2)**(-float(3)/2)*inflow
        else:
            f = inflow - mu*np.tan(alpha) - CT/(2*np.sqrt(mu**2+inflow**2))
            f_prime = 1 + CT/2/B**2*(mu**2+inflow**2)**(-float(3)/2)*inflow
        inflow -= f/f_prime
        converged = abs((inflow - inflow_old)/inflow) < 0.0005
        i += 1
    inflow_span = inflow * np.ones(n_elements)
    if tip_loss:
        inflow_span[-1] *= 3.
    return inflow_span


def axial_flight(local_solidity, prop, lambda_c, local_angle, alpha0, Clalpha, v_tip, v_climb, omega, r, blade_rad, F, spd_snd,
                 mach_corr):
    u_t = omega * r * blade_rad
    local_mach = u_t / spd_snd
    this_Clalpha = np.array(Clalpha)

    # Mach number correction
    if mach_corr:
        this_Clalpha /= np.sqrt(1 - local_mach**2)
    with np.errstate(invalid='raise'):
        try:
            local_inflow = np.sqrt((local_solidity*this_Clalpha/(16*F) - lambda_c/2)**2 +
                                   (local_solidity*this_Clalpha/(8*F)*(local_angle-alpha0)*r)) - \
                           (local_solidity*this_Clalpha/(16*F)) + (lambda_c/2)
        except FloatingPointError:
            print "FP error in inflow calculation"
            raise
        if any(np.isnan(local) for local in local_inflow):
            print "Local inflow has nan value"
            raise FloatingPointError

    v_induced = local_inflow * v_tip - v_climb
    u_p = v_climb + v_induced

    u_resultant = np.sqrt(u_p**2 + u_t**2)
    rel_inflow_angle = np.arctan(u_p / u_t)
    if any(u_resultant[u_resultant > 10000.]):
        print "u_res = " + str(u_resultant)
        print "u_t = " + str(u_t)
        print "u_p = " + str(u_p)
        print "v_tip = " + str(v_tip)
        print "v_induced = " + str(v_induced)
        print "omega = " + str(omega*60/2/np.pi)
        print "r = " + str(r)
        print "R = " + str(blade_rad)

    return local_inflow, rel_inflow_angle, u_resultant