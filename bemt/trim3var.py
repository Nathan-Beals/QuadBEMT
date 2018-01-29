import numpy as np
from bemt import bemt_forward_flight


def equilibrium(x, v_inf, quadrotor):
        alpha, omega1, omega2 = x
        prop1, prop2 = quadrotor.propellers
        quad_length = quadrotor.length
        frame_drag = quadrotor.frame_drag(alpha)
        frame_lift = quadrotor.frame_lift(alpha)
        thrust1, rotor_drag1, rotor_power1 = bemt_forward_flight(prop1, pitch=0., omega=omega1, alpha=alpha,
                                                                 v_inf=v_inf, n_azi_elements=20)
        thrust2, rotor_drag2, rotor_power2 = bemt_forward_flight(prop2, pitch=0., omega=omega2, alpha=alpha,
                                                                 v_inf=v_inf, n_azi_elements=20)
        frame_normal = frame_drag*np.sin(alpha) + frame_lift*np.cos(alpha)    # Frame lift is positive down frame_normal is positive down
        f1 = thrust1*quad_length/2 + frame_normal*quad_length/4 - thrust2*quad_length/2
        f2 = thrust1*np.sin(alpha) + thrust2*np.sin(alpha) - rotor_drag1*np.cos(alpha) - rotor_drag2*np.cos(alpha) - \
             frame_drag
        f3 = thrust1*np.cos(alpha) + thrust2*np.cos(alpha) + rotor_drag1*np.sin(alpha) + rotor_drag2*np.sin(alpha) - \
             quadrotor.weight
        return f1, f2, f3


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
    alpha0, omega10, omega20 = x0     # tilt angle in radians, rotational speed in rad/s
    x = np.array([alpha0, omega10, omega20])
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
