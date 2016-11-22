import numpy as np


def uniform_ff(CT_target, alpha, mu, n_elements):
    inflow = np.sqrt(CT_target/2)
    converged = False
    while not converged:
        inflow_old = inflow
        f = inflow - mu*np.tan(alpha) - CT_target/(2*np.sqrt(mu**2+inflow))
        f_prime = 1 + CT_target/2*(mu**2+inflow**2)**(-float(3)/2)*inflow
        inflow -= f/f_prime
        converged = abs((inflow - inflow_old)/inflow) < 0.0005
    return inflow * np.ones(n_elements)


def axial_flight(local_solidity, prop, lambda_c, local_angle, v_tip, v_climb, omega, r, blade_rad, F, spd_snd,
                 mach_corr):
    u_t = omega * r * blade_rad
    local_mach = u_t / spd_snd
    Clalpha = prop.get_Clalpha()
    # Mach number correction
    if mach_corr:
        Clalpha /= np.sqrt(1 - local_mach**2)
    with np.errstate(invalid='raise'):
        try:
            local_inflow = np.sqrt((local_solidity*Clalpha/(16*F) - lambda_c/2)**2 +
                               (local_solidity*Clalpha/(8*F)*local_angle*r)) - (local_solidity*Clalpha/(16*F)) + \
                               (lambda_c/2)
            if any(np.isnan(local) for local in local_inflow):
                raise FloatingPointError
        except FloatingPointError:
            print "local_solidity = %s" % str(local_solidity)
            print "Clalpha = %s" % str(Clalpha)
            print "lambda_c = %s" % str(lambda_c)
            print "theta = %s" % str(local_angle)
            print "r = %s" % str(r)
            print "F = %s" % str(F)
            raise
        # print "local_inflow = %s" % local_inflow

    v_induced = local_inflow * v_tip - v_climb
    u_p = v_climb + v_induced

    u_resultant = np.sqrt(u_p**2 + u_t**2)
    rel_inflow_angle = np.arctan(u_p / u_t)
    return local_inflow, rel_inflow_angle, u_resultant