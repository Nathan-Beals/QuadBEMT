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


def bemt_forward_flight(quadrotor, pitch, omega, alpha, v_inf, n_azi_elements, alt=0, tip_loss=True,
                        mach_corr=False, inflow_model='uniform', allowable_Re=[], Cl_funs={}, Cd_funs={}):

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


    # Define some other parameters for use in calculations
    v_tip = blade_rad * omega   # Blade tip speed
    mu = v_inf*np.cos(alpha)/(omega*blade_rad)  # Advance ratio
    # Define initial guess of thrust required using the hover value for the current omega
    dT = np.empty(n_elements)
    dP = np.empty(n_elements)
    P = 0.0
    T = 0.0
    H = 0.0

    # Now do the forward flight case. Since the velocity components normal and in plane with the TPP are now a
    # function of the blade azimuth angle, psi, we need to average forces over the entire range of psi.
    psi = np.linspace(0, 2*np.pi, n_azi_elements)
    dpsi = 2 * np.pi * y / n_azi_elements   # size of d_psi for each annulus

    def bemt_eval(T_old):
        CT_guess = T_old / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)
        if inflow_model == 'uniform':
            local_inflow = inflow.uniform_ff(CT_guess, alpha, mu, n_elements, tip_loss=tip_loss)
        else:
            local_inflow = 0

        dT_mat = np.empty([n_azi_elements, n_elements], dtype=float)
        dH_mat = np.empty([n_azi_elements, n_elements], dtype=float)
        dQ_mat = np.empty([n_azi_elements, n_elements], dtype=float)
        for i_azi, azi_ang in enumerate(psi):
            u_p = local_inflow * v_tip
            u_t = omega*y + v_inf*np.cos(alpha)*np.sin(azi_ang)
            U = (u_t**2 + u_p**2)**0.5
            local_mach = U / spd_snd
            rel_inflow_angle = np.arctan(u_p / u_t)
            eff_aoa = local_angle - rel_inflow_angle

            # Calculate Reynolds number along the span of the blade
            Re = U * chord / kine_visc
            # Re_actual = np.array(Re)
            if allowable_Re:
                Re = np.array([min(allowable_Re, key=lambda x: abs(x-rn)) for rn in Re])

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

            # Mach number correction to the lift coefficient
            if mach_corr:
                Cl /= np.sqrt(1 - local_mach**2)

            # Calculate sectional lift and drag
            dL = 0.5 * dens * U**2 * chord * Cl * dy
            dD = 0.5 * dens * U**2 * chord * Cd * dy

            # Calculate sectional (wrt blade reference frame) normal and in-plane forces
            dFz = dL * np.cos(rel_inflow_angle) - dD * np.sin(rel_inflow_angle)
            dFx = dL * np.sin(rel_inflow_angle) + dD * np.cos(rel_inflow_angle)

            # Calculate rotor thrust and drag and torque with respect to the rotor frame
            # (i.e., drag, H, is rotor in-plane)
            dT = dFz
            dH = dFx * np.sin(azi_ang)
            dQ = y * dFx
            dT_mat[i_azi:] = dT
            dH_mat[i_azi:] = dH
            dQ_mat[i_azi:] = dQ

            # # Calculate the lift and drag in the observer frame for a flight path angle of zero.
            # dL_observer = dH * np.sin(alpha) + dT * np.cos(alpha)
            # dD_observer = dT * np.sin(alpha) - dH * np.cos(alpha)
            # dL_observer_mat[i_azi:] = dL_observer
            # dD_observer_mat[i_azi:] = dD_observer

        # Calculate total thrust and thrust coefficient
        dT = np.mean(dT_mat, axis=0) * n_blades
        T = sum(dT)
        CT = T / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)

        # Calculate total rotor frame drag, H
        dH = np.mean(dH_mat, axis=0) * n_blades
        H = sum(dH)
        CH = H / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)

        # # Calculate observer frame lift and drag
        # dL_observer = np.mean(dL_observer_mat, axis=0) * dy * psi * n_blades
        # dD_observer = np.mean(dD_observer_mat, axis=0) * dy * psi * n_blades
        # L_observer = sum(dL_observer)
        # D_observer = sum(dD_observer)

        # Calculate rotor torque, and power
        dQ = np.mean(dQ_mat, axis=0) * n_blades
        dP = dQ * omega
        P = sum(dP)
        # Q = sum(dQ)
        # CP = P / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**3)
        # CQ = Q / (dens * np.pi * blade_rad**3 * (omega*blade_rad)**2)
        return T, H, P

    t = sum(bemt_axial(propeller, pitch, omega, allowable_Re=allowable_Re, Cl_funs=Cl_funs, Cd_funs=Cd_funs)[0])
    bemt_longoutput = []
    converged = False
    i = 0
    max_i = 10
    while not converged and i < max_i:
        tp = 1e-8
        bemt_longoutput = bemt_eval(t)
        ft = t - bemt_longoutput[0]
        ftp = (t + tp) - bemt_eval(t + tp)[0]
        ftprime = (ftp - ft) / tp
        tnew = t - ft / ftprime
        converged = abs((tnew - t)/tnew) < 0.0005
        t = tnew
        i += 1
    return bemt_longoutput[0], bemt_longoutput[1], bemt_longoutput[2], converged


def bemt_axial(propeller, pitch, omega, allowable_Re=[], Cl_funs={}, Cd_funs={}, v_climb=0, alt=0, tip_loss=True, mach_corr=False, output='short'):
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
    Re_approx = u_t * chord / kine_visc
    Clalpha, Cl0, alpha0 = aero_coeffs.get_liftCurveInfo(Re_approx, Cl_table)
    v_tip = blade_rad * omega
    lambda_c = v_climb/v_tip

    # Now handle hover and vertical flight cases
    # First calculate inflow along span by using F = 1 to get initial value not including tip loss
    F = 1
    local_inflow, rel_inflow_angle, u_resultant = inflow.axial_flight(local_solidity, propeller, lambda_c, local_angle,
                                                                      alpha0, Clalpha, v_tip, v_climb, omega, r,
                                                                      blade_rad, F, spd_snd, mach_corr=mach_corr)
    #print "inflow no tip loss = " + str(local_inflow)
    inflow_notl = local_inflow
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
            try:
                local_inflow, rel_inflow_angle, u_resultant = inflow.axial_flight(local_solidity, propeller, lambda_c,
                                                                                  local_angle, alpha0, Clalpha, v_tip,
                                                                                  v_climb, omega, r, blade_rad, F,
                                                                                  spd_snd, mach_corr=mach_corr)
            except FloatingPointError:
                raise
            converged = abs((local_inflow - local_inflow_old)/local_inflow) < 0.0005
            i += 1
    #print "inflow tip loss = " + str(local_inflow)
    #print "tl ratio = " + str(local_inflow / inflow_notl)

    # Calculate Reynolds number along the span of the blade
    Re = u_resultant * chord / kine_visc
    Re_actual = np.array(Re)
    if allowable_Re:
        Re = np.array([min(allowable_Re, key=lambda x: abs(x-rn)) for rn in Re])

    # Now calculate the effective angle of attack at the blade stations.
    eff_aoa = local_angle - rel_inflow_angle
    #print "hover effaoa = " + str(eff_aoa)

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
    #print "hover Cl = " + str(Cl)
    # Calculate forces
    dL = 0.5*dens*u_resultant**2*chord*Cl*dy
    dD = 0.5*dens*u_resultant**2*chord*Cd*dy

    dFz = dL*np.cos(rel_inflow_angle) - dD*np.sin(rel_inflow_angle)
    dFx = dD*np.cos(rel_inflow_angle) + dL*np.sin(rel_inflow_angle)

    dT = n_blades * dFz
    dQ = n_blades * dFx * y
    dP = n_blades * dFx * omega * y

    dPp = n_blades * dD*np.cos(rel_inflow_angle) * omega * y
    dPo = n_blades * dL*np.sin(rel_inflow_angle) * omega * y

    T = sum(dT)
    Q = sum(dQ)
    P = sum(dP)
    Pp = sum(dPp)
    Po = sum(dPo)

    CT = T / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)
    CP = P / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**3)
    CQ = Q / (dens * np.pi * blade_rad**3 * (omega*blade_rad)**2)

    prop_CT = T / (dens * (omega/2/np.pi)**2 * (blade_rad*2)**4)
    prop_CP = P / (dens * (omega/2/np.pi)**3 * (blade_rad*2)**5)
    FM = prop_CT**(3./2)/np.sqrt(2)/prop_CP

    if output == 'short':
        return dT, P
    return dT, dP, Cd, Cl, u_resultant, chord, dL, local_inflow, rel_inflow_angle, eff_aoa, dFx, dFz, Re, prop_CT, prop_CP, Pp, Po, Re_actual





