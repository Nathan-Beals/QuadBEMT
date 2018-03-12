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
        alpha, omega = x
        n_azi_elements = kwargs['n_azi_elements']
        alt = kwargs['alt']
        dens = get_dens(alt)
        thrust, rotor_drag, rotor_power = bemt_forward_flight(quadrotor.propeller, pitch=0., omega=omega,
                                                              alpha=alpha, v_inf=v_inf, n_azi_elements=n_azi_elements)
        f1 = thrust*np.sin(alpha) - rotor_drag*np.cos(alpha) - quadrotor.frame_drag(alpha, v_inf, dens)
        f2 = thrust*np.cos(alpha) + rotor_drag*np.sin(alpha) - quadrotor.weight/4
        return f1, f2


def jacobian(fun, x):
        n = len(x)
        fun_x = fun(x)
        jac = np.zeros([n, len(fun_x)])
        eps = 1e-8
        x_perturb = x
        for i in xrange(n):
            x_perturb[i] += eps
            jac[i] = (fun(x_perturb) - fun_x) / eps
            x_perturb[i] = x[i]
        return np.transpose(jac)


def trim(quadrotor, v_inf, x0, kwargs):
    e = 0.0005
    alpha0, omega0 = x0     # tilt angle in radians, rotational speed in rad/s
    x = np.array([alpha0, omega0])
    converged = False
    max_i = 20
    i = 0
    while not converged and i <= max_i:
        if i == max_i:
            return []
        i += 1
        equilibrium_eval = np.array(equilibrium(x, v_inf, quadrotor, kwargs))
        equilibrium_jac = np.array(jacobian(equilibrium, x))
        xnew = np.subtract(x, np.linalg.solve(equilibrium_jac, equilibrium_eval))
        converged = all(abs((xnew[i] - x[i])/xnew[i]) < e for i in xrange(len(xnew)))
        x = xnew
    return x

