import numpy as np


def uniform_ff(CT, alpha, mu, n_elements):
    with np.errstate(invalid='raise'):
        try:
            inflow = np.sqrt(CT/2)
        except FloatingPointError:
            print "FP error in ff inflow"
            raise
    converged = False
    max_i = 100
    i = 0
    while not converged and i < max_i:
        inflow_old = inflow
        f = inflow - mu*np.tan(alpha) - CT/(2*np.sqrt(mu**2+inflow**2))
        f_prime = 1 + CT/2**2*(mu**2+inflow**2)**(-float(3)/2)*inflow
        inflow -= f/f_prime
        converged = abs((inflow - inflow_old)/inflow) < 0.005
        i += 1
    inflow_span = inflow * np.ones(n_elements)
    return inflow_span


def new_ff(local_inflow, local_solidity, lambda_c, local_angle, alpha, mu, alpha0, Clalpha, r, F, local_mach,
           mach_corr=False):

    this_Clalpha = np.array(Clalpha)
    mu_array = np.ones(len(r))*mu

    if mach_corr:
        this_Clalpha /= np.sqrt(1 - local_mach**2)

    with np.errstate(invalid='raise'):
        try:
            hover_inflow = axial_flight(local_solidity, lambda_c, local_angle, alpha0, Clalpha, r, F, local_mach,
                                        mach_corr=mach_corr)
            local_inflow = mu_array*np.tan(alpha) + hover_inflow**2/np.sqrt(mu_array**2+local_inflow**2)
        except FloatingPointError:
            print "FP error in ff inflow calculation"
            raise
        if any(np.isnan(local) for local in local_inflow):
            print "Local inflow has nan value"
            raise FloatingPointError

    return local_inflow


def axial_flight(local_solidity, lambda_c, local_angle, alpha0, Clalpha, r, F, local_mach, mach_corr=False):

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
            raise
        if any(np.isnan(local) for local in local_inflow):
            raise FloatingPointError
    return local_inflow
