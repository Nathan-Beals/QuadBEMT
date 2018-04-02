import numpy as np
import inflow
import aero_coeffs

# Define some constants
RAD_EARTH = 6371000    # Radius of Earth in meters
GAMMA = 1.4    # Heat capacity ratio for air
R = 287    # For air R = 287 J/kg K
A = -0.0065    # lapse rate in /m
Gs = -9.8    # Acceleration due to gravity in m/s**2
TEMP_SSL = 288.16    # Standard sea level temp in K
PRESS_SSL = 1.01325 * 10**5    # Standard sea level pressure in N/m**2
kine_visc = 1.460 * 10**-5  # Kinematic viscosity of air


def bemt_forward_flight(quadrotor, pitch, omega, alpha, v_inf, n_azi_elements, alt=0, tip_loss=True, mach_corr=False,
                        allowable_Re=[], Cl_funs={}, Cd_funs={}, lift_curve_info_dict={}):

    alt_geop = (RAD_EARTH*alt)/(RAD_EARTH+alt)
    temp = TEMP_SSL - A * alt_geop
    press = PRESS_SSL * (temp/TEMP_SSL)**(-Gs/(A*R))
    dens = press/(R*temp)
    spd_snd = np.sqrt(GAMMA * R * temp)

    # Define blade geometry parameters. Pitch, chord, and r are all lists of the same length defining the blade
    # geometry at a specific span location r
    propeller = quadrotor.propeller
    n_blades = propeller.n_blades
    blade_rad = propeller.radius
    twist = np.array(propeller.twist)
    local_angle = pitch + twist
    chord = np.array(propeller.chord)
    dy = propeller.dy
    dr = propeller.dr
    r = np.array(propeller.r)
    y = np.array(propeller.y)
    n_elements = len(r)
    local_solidity = np.array(propeller.solidity)
    airfoil = propeller.airfoils[0][0]
    Cl_table = propeller.Cl_tables[airfoil]
    Cd_table = propeller.Cd_tables[airfoil]
    allowable_Re = allowable_Re
    Cl_funs = Cl_funs
    Cd_funs = Cd_funs
    lift_curve_info_dict = lift_curve_info_dict
    lambda_c = 0.

    # Define some other parameters for use in calculations
    v_tip = blade_rad * omega   # Blade tip speed
    mu = v_inf*np.cos(alpha)/(omega*blade_rad)  # Advance ratio

    # Now do the forward flight case. Since the velocity components normal and in plane with the TPP are now a
    # function of the blade azimuth angle, psi, we need to average forces over the entire range of psi.
    psi = np.linspace(0, 2*np.pi, n_azi_elements+1)[:-1]
    dpsi = 2 * np.pi * y / n_azi_elements   # size of d_psi for each annulus

    dT_mat = np.empty([n_azi_elements, n_elements], dtype=float)
    dH_mat = np.empty([n_azi_elements, n_elements], dtype=float)
    dP_mat = np.empty([n_azi_elements, n_elements], dtype=float)
    # dPo_mat = np.empty([n_azi_elements, n_elements], dtype=float)
    # dPi_mat = np.empty([n_azi_elements, n_elements], dtype=float)
    # inflow_mat = np.empty([n_azi_elements, n_elements], dtype=float)

    for i_azi, azi_ang in enumerate(psi):
        u_t = omega*y + v_inf*np.cos(alpha)*np.sin(azi_ang)
        local_mach = u_t / spd_snd
        Re_approx = u_t * chord / kine_visc
        closest_Re = aero_coeffs.closest_Re(Re_approx, allowable_Re)
        Clalpha = np.zeros(len(r))
        Cl0 = np.zeros(len(r))
        alpha0 = np.zeros(len(r))
        for i, Re in enumerate(closest_Re):
            Clalpha[i], Cl0[i], alpha0[i] = lift_curve_info_dict[Re]

        ############################################################################################################
        #                                  Begin inflow calculation                                                #
        ############################################################################################################

        # First calculate the hover inflow value without tip loss for the current value of omega
        F = 1
        hover_inflow = inflow.axial_flight(local_solidity, lambda_c, local_angle, alpha0, Clalpha, r, F, local_mach,
                                           mach_corr=mach_corr)
        # Now use the guess to calculate a guess of the forward flight inflow without tip loss
        local_inflow = mu*np.tan(alpha)*np.ones(len(r)) + hover_inflow**2/np.sqrt(mu**2+hover_inflow**2)
        rel_inflow_angle = np.arctan(local_inflow*v_tip / u_t)

        converged = np.array([False]*n_elements)
        i = 0
        max_i = 25
        eps = 1e-8
        while not all(converged) and i < max_i:

            local_inflow_old = np.array(local_inflow)

            if tip_loss:
                f_tip = n_blades/2. * ((1 - r)/(r * rel_inflow_angle))
                f = f_tip
                f[-1] = 0.0000000000001
                F = (2/np.pi) * np.arccos(np.exp(-f))

            this_Clalpha = np.array(Clalpha)
            mu_array = np.ones(len(r))*mu

            if mach_corr:
                this_Clalpha /= np.sqrt(1 - local_mach**2)

            with np.errstate(invalid='raise'):
                try:
                    hover_inflow = np.sqrt((local_solidity*this_Clalpha/(16*F) - lambda_c/2)**2 +
                                           (local_solidity*this_Clalpha/(8*F)*(local_angle-alpha0)*r)) - \
                                           (local_solidity*this_Clalpha/(16*F)) + (lambda_c/2)
                    ff_inflow = mu_array*np.tan(alpha) + hover_inflow**2/np.sqrt(mu_array**2+local_inflow**2)
                    ff_inflow_perturb = mu_array*np.tan(alpha) + hover_inflow**2/np.sqrt(mu_array**2+(local_inflow+eps)**2)
                except FloatingPointError:
                    print "FP error in ff inflow calculation"
                    raise
                if any(np.isnan(local) for local in local_inflow):
                    print "Local inflow has nan value"
                    raise FloatingPointError

            f_lambda = local_inflow - ff_inflow
            f_perturb = (local_inflow + eps) - ff_inflow_perturb
            f_lambda_prime = (f_perturb - f_lambda) / eps
            local_inflow -= f_lambda/f_lambda_prime

            converged = abs((local_inflow - local_inflow_old)/local_inflow) < 0.005

            u_p = local_inflow * v_tip
            rel_inflow_angle = np.arctan(u_p / u_t)
            u_resultant = np.sqrt(u_p**2 + u_t**2)

            i += 1
        if not all(converged):
            raise FloatingPointError
        #############################################################################################################
        #inflow_mat[i_azi:] = local_inflow
        eff_aoa = local_angle - rel_inflow_angle
        # if np.any(eff_aoa < -5.*np.pi/180):
        #     raise IndexError

        # Calculate Reynolds number along the span of the blade
        Re = u_resultant * chord / kine_visc
        if allowable_Re:
            Re = aero_coeffs.closest_Re(Re, allowable_Re)

        # If we are using only a selection of discrete Reynolds numbers for the sake of calculating lift and
        # drag coefficients, calculate using the linearized Cl functions and polynomial Cd functions found in
        # aero_coeffs.py.
        Cl = []
        Cd = []
        if Cl_funs:
            for i in xrange(len(eff_aoa)):
                Cl.append(float(Cl_funs[Re[i]](eff_aoa[i])))
                Cd.append(Cd_funs[Re[i]](eff_aoa[i]))
            Cl = np.array(Cl)
            Cd = np.array(Cd)
        # Otherwise use the full table interpolation
        else:
            Cl = np.nan_to_num(np.array(aero_coeffs.get_Cl(eff_aoa, Re, Cl_table)))
            Cd = np.nan_to_num(np.array(aero_coeffs.get_Cd(eff_aoa, Re, Cd_table)))

        # Calculate sectional lift and drag
        dL = 0.5 * dens * u_resultant**2 * chord * Cl * dy
        if tip_loss:
            dL[-1] = 0
        dD = 0.5 * dens * u_resultant**2 * chord * Cd * dy

        # Calculate sectional (wrt blade reference frame) normal and in-plane forces
        dFz = dL * np.cos(rel_inflow_angle) - dD * np.sin(rel_inflow_angle)
        dFx = dL * np.sin(rel_inflow_angle) + dD * np.cos(rel_inflow_angle)

        # Calculate rotor thrust and drag and torque with respect to the rotor frame
        # (i.e., drag, H, is rotor in-plane)
        dT = dFz
        dH = dFx * np.sin(azi_ang)
        dP = y * dFx * omega
        dT_mat[i_azi:] = dT
        dH_mat[i_azi:] = dH
        dP_mat[i_azi:] = dP

        # dPo = y * dD * np.cos(rel_inflow_angle) * omega
        # dPi = y * dL * np.sin(rel_inflow_angle) * omega
        # dPo_mat[i_azi:] = dPo
        # dPi_mat[i_azi:] = dPi

    # Calculate total thrust and thrust coefficient
    dT = np.mean(dT_mat, axis=0) * n_blades
    T = sum(dT)
    CT = T / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)

    # Calculate total rotor frame drag, H
    dH = np.mean(dH_mat, axis=0) * n_blades
    H = sum(dH)
    CH = H / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)

    # Calculate rotor torque, and power
    dP = np.mean(dP_mat, axis=0) * n_blades
    P = sum(dP)

    # dPo = np.mean(dPo_mat, axis=0) * n_blades
    # dPi = np.mean(dPi_mat, axis=0) * n_blades
    # print "Profile power dist = " + str(dPo)
    # print "Total profile power = " + str(sum(dPo))
    # print "Induced power dist = " + str(dPi)
    # print "Total induced power = " + str(sum(dPi))

    # Q = sum(dQ)
    # CP = P / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**3)
    # CQ = Q / (dens * np.pi * blade_rad**3 * (omega*blade_rad)**2)
    return T, H, P#, inflow_mat


def bemt_axial(propeller, pitch, omega, allowable_Re=[], Cl_funs={}, Cd_funs={}, v_climb=0, alt=0, tip_loss=True,
               mach_corr=False, output='short'):
    # Calculate geopotential altitude
    alt_geop = (RAD_EARTH*alt)/(RAD_EARTH+alt)

    # Calculate atmospheric conditions
    temp = TEMP_SSL - A * alt_geop
    press = PRESS_SSL * (temp/TEMP_SSL)**(-Gs/(A*R))
    dens = press/(R*temp)
    if alt != 0:
        dens = 0.8194
    spd_snd = np.sqrt(GAMMA * R * temp)

    # Define blade geometry parameters. Pitch, chord, and r are all lists of the same length defining the blade
    # geometry at a specific span location r
    n_blades = propeller.n_blades
    blade_rad = propeller.radius
    twist = np.array(propeller.twist)
    local_angle = pitch + twist
    chord = np.array(propeller.chord)
    dy = propeller.dy
    dr = propeller.dr
    r = np.array(propeller.r)
    y = np.array(propeller.y)
    n_elements = len(r)
    local_solidity = np.array(propeller.solidity)
    airfoil = propeller.airfoils[0][0]
    Cl_table = propeller.Cl_tables[airfoil]
    Cd_table = propeller.Cd_tables[airfoil]

    # Due to the possible camber of the airfoils along the span, we need to correct the local angle to include the zero
    # lift angle of attack. For positively cambered airfoils this will be a negative angle (all values of alpha0 will be
    # negative. Also find the lift curve slope along the span of the blade. Both quantities are calculated using an
    # approximate Reynolds number which calculates Re using only the in-plane portion of the freestream velocity.
    u_t = omega * r * blade_rad
    local_mach = u_t / spd_snd
    Re_approx = u_t * chord / kine_visc
    Clalpha, Cl0, alpha0 = aero_coeffs.get_liftCurveInfo(Re_approx, Cl_table)
    v_tip = blade_rad * omega
    lambda_c = v_climb/v_tip

    # Now handle hover and vertical flight cases
    # First calculate inflow along span by using F = 1 to get initial value not including tip loss
    F = 1
    local_inflow = inflow.axial_flight(local_solidity, lambda_c, local_angle, alpha0, Clalpha, r, F, local_mach,
                                       mach_corr=mach_corr)
    v_induced = local_inflow * v_tip - v_climb
    u_p = v_climb + v_induced
    rel_inflow_angle = np.arctan(u_p / u_t)
    u_resultant = np.sqrt(u_p**2 + u_t**2)

    # Now if tip_loss correction is desired, use the F = 1 solution as a starting guess to find the inflow
    if tip_loss:
        converged = np.array([False]*n_elements)
        i = 0
        max_i = 25
        while not all(converged) and i < max_i:
            local_inflow_old = local_inflow
            f_tip = n_blades/2. * ((1 - r)/(r * rel_inflow_angle))
            f = f_tip
            f[-1] = 0.0000000000001
            F = (2/np.pi) * np.arccos(np.exp(-f))
            #F[np.isnan(F)] = 0.999999999
            try:
                local_inflow = inflow.axial_flight(local_solidity, lambda_c, local_angle, alpha0, Clalpha, r, F,
                                                   local_mach, mach_corr=mach_corr)
                v_induced = local_inflow * v_tip - v_climb
                u_p = v_climb + v_induced
                u_resultant = np.sqrt(u_p**2 + u_t**2)
                rel_inflow_angle = np.arctan(u_p / u_t)
            except FloatingPointError:
                raise
            converged = abs((local_inflow - local_inflow_old)/local_inflow) < 0.0005
            i += 1

    # Calculate Reynolds number along the span of the blade
    Re = u_resultant * chord / kine_visc
    if allowable_Re:
        Re = aero_coeffs.closest_Re(Re, allowable_Re)

    # Now calculate the effective angle of attack at the blade stations.
    eff_aoa = local_angle - rel_inflow_angle

    # If we are using only a selection of discrete Reynolds numbers for the sake of calculating lift and drag
    # coefficients, calculate using the linearized Cl functions and polynomial Cd functions found in aero_coeffs.py.
    Cl = []
    Cd = []
    if Cl_funs:
        for i in xrange(len(eff_aoa)):
            Cl.append(float(Cl_funs[Re[i]](eff_aoa[i])))
            Cd.append(Cd_funs[Re[i]](eff_aoa[i]))
        Cl = np.array(Cl)
        Cd = np.array(Cd)
    # Otherwise use the full table interpolation
    else:
        Cl = np.nan_to_num(np.array(aero_coeffs.get_Cl(eff_aoa, Re, Cl_table)))
        Cd = np.nan_to_num(np.array(aero_coeffs.get_Cd(eff_aoa, Re, Cd_table)))

    # Calculate forces
    dL = 0.5*dens*u_resultant**2*chord*Cl*dy
    # if tip_loss:
    #     dL[-1] = 0.
    dD = 0.5*dens*u_resultant**2*chord*Cd*dy

    dFz = dL*np.cos(rel_inflow_angle) - dD*np.sin(rel_inflow_angle)
    dFx = dD*np.cos(rel_inflow_angle) + dL*np.sin(rel_inflow_angle)

    dT = n_blades * dFz
    # dQ = n_blades * dFx * y
    dP = n_blades * dFx * omega * y

    #dPp = n_blades * dD*np.cos(rel_inflow_angle) * omega * y
    #dPo = n_blades * dL*np.sin(rel_inflow_angle) * omega * y

    #T = sum(dT)
    #Q = sum(dQ)
    P = sum(dP)
    # Pp = sum(dPp)
    # Po = sum(dPo)
    #
    # CT = T / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)
    # CP = P / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**3)
    # CQ = Q / (dens * np.pi * blade_rad**3 * (omega*blade_rad)**2)
    #
    # prop_CT = T / (dens * (omega/2/np.pi)**2 * (blade_rad*2)**4)
    # prop_CP = P / (dens * (omega/2/np.pi)**3 * (blade_rad*2)**5)
    #FM = prop_CT**(3./2)/np.sqrt(2)/prop_CP

    if output == 'short':
        return dT, P
    return dT, P





