import numpy as np
from bemt import bemt_forward_flight


def equilibrium(x, v_inf, quadrotor):
        alpha, omega = x
        thrust, rotor_drag, rotor_power = bemt_forward_flight(quadrotor.propeller, pitch=0., omega=omega,
                                                              alpha=alpha, v_inf=v_inf, n_azi_elements=20)
        f1 = thrust*np.sin(alpha) - rotor_drag*np.cos(alpha) - quadrotor.frame_drag(alpha, v_inf)
        f2 = thrust*np.cos(alpha) + rotor_drag*np.sin(alpha) - quadrotor.weight
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


def trim(quadrotor, v_inf, x0):
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
        equilibrium_eval = np.array(equilibrium(x, v_inf, quadrotor))
        equilibrium_jac = np.array(jacobian(equilibrium, x))
        xnew = np.subtract(x, np.linalg.solve(equilibrium_jac, equilibrium_eval))
        converged = all(abs((xnew[i] - x[i])/xnew[i]) < e for i in xrange(len(xnew)))
        x = xnew
    return x

