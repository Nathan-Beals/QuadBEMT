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


def equilibrium(x, v_inf, quadrotor, kwargs):
        print "entering equilibrium"
        alpha, omega = x
        n_azi_elements = kwargs['n_azi_elements']
        alt = kwargs['alt']
        dens = get_dens(alt)
        try:
            thrust, rotor_drag, rotor_power, bemt_converged = bemt_forward_flight(quadrotor, pitch=0.,
                                                                                  omega=omega, alpha=alpha, v_inf=v_inf,
                                                                                  n_azi_elements=n_azi_elements)
            print "thrust = " + str(thrust)
        except Exception as e:
            print "{} in ff bemt".format(type(e).__name__)
            raise
        f1 = thrust*np.sin(alpha) - rotor_drag*np.cos(alpha) - quadrotor.frame_drag(alpha, v_inf, dens)
        f2 = thrust*np.cos(alpha) + rotor_drag*np.sin(alpha) - quadrotor.weight/4
        print "leaving equilibrium"
        return f1, f2, bemt_converged


def jacobian(fun, x, vinf, quadrotor, kwargs):
        n = len(x)
        fun_x = np.array(fun(x, vinf, quadrotor, kwargs)[:-1])
        print "type fun_x is " + str(type(fun_x))
        jac = np.zeros([n, len(fun_x)])
        print jac
        eps = 1e-8
        x_perturb = x
        for i in xrange(n):
            x_perturb[i] += eps
            perturb_eval = np.array(fun(x_perturb, vinf, quadrotor, kwargs)[:-1])
            jac[i] = (perturb_eval - fun_x) / eps
            x_perturb[i] = x[i]
        return np.transpose(jac)


def trim(quadrotor, v_inf, x0, kwargs):
    e = 0.0005
    alpha0, omega0 = x0     # tilt angle in radians, rotational speed in rad/s
    x = np.array([alpha0, omega0])
    converged = False
    max_i = 10
    i = 0
    while not converged and i <= max_i:
        if i == max_i:
            return x[0], x[1], converged
        equilibrium_eval = np.array(equilibrium(x, v_inf, quadrotor, kwargs))
        print "equilibrium_eval = " + str(equilibrium_eval)
        if not equilibrium_eval[-1]:
            return x, converged
        equilibrium_jac = np.array(jacobian(equilibrium, x, v_inf, quadrotor, kwargs))
        print "equilibrium_jac = " + str(equilibrium_jac)
        xnew = x - np.linalg.solve(equilibrium_jac, equilibrium_eval[:-1])
        converged = all(abs((xnew[i] - x[i])/xnew[i]) < e for i in xrange(len(xnew)))
        x = xnew
        print "(alpha, omega) = (%f, %f)" % (x[0], x[1])
        i += 1
    if not converged:
        print "trim did not converge"
    return x[0], x[1], converged

