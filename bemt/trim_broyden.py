import numpy as np
from bemt import bemt_forward_flight


def get_dens(alt):
    RAD_EARTH = 6371000    # Radius of Earth in meters
    GAMMA = 1.4    # Heat capacity ratio for air
    R = 287    # For air R = 287 J/kg K
    A = -0.0065    # lapse rate in /m
    Gs = -9.8    # Acceleration due to gravity in m/s**2
    TEMP_SSL = 288.16    # Standard sea level temp in K
    PRESS_SSL = 1.01325 * 10**5    # Standard sea level pressure in N/m**2
    kine_visc = 1.460 * 10**-5  # Kinematic viscosity of air
    alt_geop = (RAD_EARTH*alt)/(RAD_EARTH+alt)
    temp = TEMP_SSL - A * alt_geop
    press = PRESS_SSL * (temp/TEMP_SSL)**(-Gs/(A*R))
    dens = press/(R*temp)
    spd_snd = np.sqrt(GAMMA * R * temp)
    return dens


def equilibrium(x, v_inf, quadrotor, kwargs, called_from='trim'):
        try:
            # If x is a 2x1 np.matrix this will handle it
            alpha = x[0, 0]
            omega = x[1, 0]
        except (IndexError, TypeError):
            # If x is a list or tuple this will handle it
            alpha, omega = x
        n_azi_elements = kwargs['n_azi_elements']
        alt = kwargs['alt']
        allowable_Re = kwargs['allowable_Re']
        Cl_funs = kwargs['Cl_funs']
        Cd_funs = kwargs['Cd_funs']
        dens = get_dens(alt)
        pitch = kwargs['pitch']
        tip_loss = kwargs['tip_loss']
        mach_corr = kwargs['mach_corr']
        try:
            thrust, rotor_drag, rotor_power, bemt_converged = bemt_forward_flight(quadrotor, pitch, omega, alpha, v_inf,
                                                                                  n_azi_elements, alt=alt,
                                                                                  tip_loss=tip_loss, mach_corr=mach_corr,
                                                                                  allowable_Re=allowable_Re,
                                                                                  Cl_funs=Cl_funs, Cd_funs=Cd_funs)
        except Exception as e:
            print "{} in ff bemt".format(type(e).__name__)
            raise
        f1 = thrust*np.sin(alpha) - rotor_drag*np.cos(alpha) - quadrotor.frame_drag(alpha, v_inf, dens)/4
        f2 = thrust*np.cos(alpha) + rotor_drag*np.sin(alpha) - quadrotor.weight/4
        return f1, f2, bemt_converged


def jacobian(fun, x, vinf, quadrotor, fun_x, kwargs):
        n = len(x)
        jac = np.zeros([n, len(fun_x)])
        eps = 1e-8
        x_perturb = [[x[0, 0]+eps, x[1, 0]], [x[0, 0], x[1, 0]+eps]]
        for i in xrange(n):
            perturb_eval = fun(x_perturb[i], vinf, quadrotor, kwargs, called_from='jacobian')[:-1]
            jac[i] = (np.array(perturb_eval) - np.array(fun_x.T).reshape(-1,)) / eps
        return np.matrix(np.transpose(jac))


def trim_broyden(quadrotor, v_inf, x0, kwargs):
    e = 0.0005
    converged = False
    max_i = 20
    i = 1

    x_old = np.matrix(x0).T         # 2x1 matrix
    fx_old = np.matrix(equilibrium(x_old, v_inf, quadrotor, kwargs)[:-1]).T    # 2x1 matrix
    jac0 = jacobian(equilibrium, x_old, v_inf, quadrotor, fx_old, kwargs)   # 2x2 np.matrix
    print "jacobian = [[%f, %f], [%f, %f]]" % (jac0[0, 0], jac0[0, 1], jac0[1, 0], jac0[1, 1])
    a_inv_old = np.linalg.inv(jac0)     # 2x2 np.matrix
    dx = -1 * a_inv_old*fx_old     # 2x1 np.matrix
    while not converged:
        if i > max_i:
            break
        x = x_old + dx      # 2x1 np.matrix
        print "xnew = (%f, %f)" % (x[0, 0], x[1, 0])
        print "xold = (%f, %f)" % (x_old[0, 0], x_old[1, 0])
        print "dx   = (%f, %f)" % (dx[0, 0], dx[1, 0])
        fx = np.matrix(equilibrium(x, v_inf, quadrotor, kwargs)).T  # 3x1 matrix
        if not fx[2, 0]:
            break
        fx = fx[:-1, 0]        # 2x1 matrix
        print "(f1, f2) = (%f, %f)" % (fx[0, 0], fx[1, 0])
        y = fx - fx_old    # 2x1 matrix
        s = x - x_old     # 2x1 matrix
        sT = s.T                        # 1x2 matrix
        a_inv = a_inv_old + ((s - a_inv_old*y)/(sT*a_inv_old*y))*sT*a_inv_old   # 2x2 matrix
        print "a_inv = [[%f, %f], [%f, %f]]" % (a_inv[0, 0], a_inv[0, 1], a_inv[1, 0], a_inv[1, 1])
        dx = -1 * a_inv*fx    # 2x1 matrix
        alpha_conv = (x[0, 0] - x_old[0, 0])/x[0, 0]
        omega_conv = (x[1, 0] - x_old[1, 0])/x[1, 0]
        print "alpha_conv = " + str(alpha_conv)
        print "omega_conv = " + str(omega_conv)
        converged = abs((x[0, 0] - x_old[0, 0])/x[0, 0]) < e and abs((x[1, 0] - x_old[1, 0])/x[1, 0]) < e
        x_old = x
        i += 1
    if not converged:
        print "trim did not converge"
    else:
        print str(i-1) + " iterations to converge"
    return x_old[0, 0], x_old[1, 0], converged
