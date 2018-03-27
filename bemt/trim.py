import numpy as np
import bemt


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
        lift_curve_info_dict = kwargs['lift_curve_info_dict']
        bemt_converged = True
        try:
            thrust, rotor_drag, rotor_power = bemt.bemt_forward_flight(quadrotor, pitch, omega, alpha, v_inf,
                                                                       n_azi_elements, alt=alt, tip_loss=tip_loss,
                                                                       mach_corr=mach_corr, allowable_Re=allowable_Re,
                                                                       Cl_funs=Cl_funs, Cd_funs=Cd_funs,
                                                                       lift_curve_info_dict=lift_curve_info_dict)
        except Exception as e:
            #print "{} in ff bemt".format(type(e).__name__)
            raise
        f1 = thrust*np.sin(alpha) - rotor_drag*np.cos(alpha) - quadrotor.frame_drag(alpha, v_inf, dens)/4
        f2 = thrust*np.cos(alpha) + rotor_drag*np.sin(alpha) - quadrotor.weight/4
        return f1, f2, bemt_converged


def jacobian(fun, x, vinf, quadrotor, fun_x, kwargs):
        n = len(x)
        jac = np.zeros([n, len(fun_x)])
        eps = 1e-8
        x_perturb = [[x[0]+eps, x[1]], [x[0], x[1]+eps]]
        for i in xrange(n):
            perturb_eval = fun(x_perturb[i], vinf, quadrotor, kwargs, called_from='jacobian')[:-1]
            jac[i] = (perturb_eval - fun_x) / eps
        return np.transpose(jac)


def trim(quadrotor, v_inf, x0, kwargs):
    e = 0.005
    alpha0, omega0 = x0     # tilt angle in radians, rotational speed in rad/s
    x = np.array([alpha0, omega0])
    converged = False
    max_i = 10
    i = 0
    while not converged and i <= max_i:
        if i == max_i:
            break
        equilibrium_x = np.array(equilibrium(x, v_inf, quadrotor, kwargs))
        if not equilibrium_x[-1]:
            break
        equilibrium_jac = jacobian(equilibrium, x, v_inf, quadrotor, equilibrium_x[:-1], kwargs)
        dx = np.linalg.solve(equilibrium_jac, -equilibrium_x[:-1])
        xnew = x + dx
        converged = all(abs((xnew[i] - x[i])/xnew[i]) < e for i in xrange(len(xnew)))
        x = xnew
        i += 1
    return x[0], x[1], converged

