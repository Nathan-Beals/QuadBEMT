import numpy as np
import inflow
from unit_conversion import rad2deg
import resource

# Define some constants
RAD_EARTH = 6371000    # Radius of Earth in meters
GAMMA = 1.4    # Heat capacity ratio for air
R = 287    # For air R = 287 J/kg K
A = -0.0065    # lapse rate in /m
Gs = -9.8    # Acceleration due to gravity in m/s**2
TEMP_SSL = 288.16    # Standard sea level temp in K
PRESS_SSL = 1.01325 * 10**5    # Standard sea level pressure in N/m**2


def bemt_forward_flight(propeller, pitch, omega, alpha, T_req, v_inf, n_azi_elements, v_climb=0, alt=0, tip_loss=True,
                        mach_corr=True, inflow_model='uniform'):
    # Calculate geopotential altitude
    alt_geop = (RAD_EARTH*alt)/(RAD_EARTH+alt)

    # Calculate atmospheric conditions
    temp = TEMP_SSL - A * alt_geop
    press = PRESS_SSL * (temp/TEMP_SSL)**(-Gs/(A*R))
    dens = press/(R*temp)
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

    # Define some other parameters for use in calculations
    v_tip = blade_rad * omega   # Blade tip speed
    mu = v_inf*np.cos(alpha)/(omega*blade_rad)  # Advance ratio
    CT_req = T_req / (dens * np.pi * blade_rad * (omega*blade_rad)**2)

    # Now do the forward flight case. Since the velocity components normal and in plane with the TPP are now a
    # function of the blade azimuth angle, psi, we need to average forces over the entire range of psi.
    psi = np.linspace(0, 2*np.pi, n_azi_elements)
    dpsi = 2 * np.pi * y / n_azi_elements   # size of d_psi for each annulus

    if inflow_model == 'uniform':
        local_inflow = inflow.uniform_ff(CT_req, alpha, mu, n_elements)
    else:
        local_inflow = 0

    dT_mat = np.empty([n_azi_elements, n_elements], dtype=float)
    dH_mat = np.empty([n_azi_elements, n_elements], dtype=float)
    dQ_mat = np.empty([n_azi_elements, n_elements], dtype=float)
    dL_observer_mat = np.empty([n_azi_elements, n_elements], dtype=float)
    dD_observer_mat = np.empty([n_azi_elements, n_elements], dtype=float)
    i_azi = 0
    for azi_ang in psi:
        u_p_nd = local_inflow
        u_t_nd = r + mu * np.sin(azi_ang)
        u_r_nd = mu * np.cos(azi_ang)
        U_nd = (u_t_nd**2 + u_p_nd**2)**0.5
        local_mach = omega*blade_rad*U_nd / spd_snd
        rel_inflow_angle = np.arctan(u_p_nd / u_t_nd)
        eff_aoa = local_angle - rel_inflow_angle
        Cl = propeller.get_Cl(eff_aoa)
        Cd = propeller.get_Cd(eff_aoa)
        # Mach number correction to the lift coefficient
        if mach_corr:
            Cl /= np.sqrt(1 - local_mach**2)
        # Calculate sectional lift and drag
        dL = 0.5 * U_nd**2 * chord * Cl
        dD = 0.5 * U_nd**2 * chord * Cd
        # Calculate sectional (wrt blade reference frame) normal and in-plane forces
        dFz = dL * np.cos(rel_inflow_angle) - dD * np.sin(rel_inflow_angle)
        dFx = dL * np.sin(rel_inflow_angle) + dD * np.cos(rel_inflow_angle)
        # Calculate rotor thrust and drag and torque with respect to the rotor frame
        # (i.e., drag, H, is rotor in-plane)
        dT = n_blades * dFz
        dH = n_blades * dFx * np.sin(azi_ang)
        dQ = n_blades * y * dFx
        dT_mat[i_azi:] = dT
        dH_mat[i_azi:] = dH
        dQ_mat[i_azi:] = dQ
        # Calculate the lift and drag in the observer frame for a flight path angle of zero.
        dL_observer = dH * np.sin(alpha) + dT * np.cos(alpha)
        dD_observer = dT * np.sin(alpha) - dH * np.cos(alpha)
        dL_observer_mat[i_azi:] = dL_observer
        dD_observer_mat[i_azi:] = dD_observer

    # Calculate total thrust and thrust coefficient
    dT = np.mean(dT_mat, axis=0) * dy * dpsi
    T = sum(dT)
    CT = T / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)

    # Calculate total rotor frame drag, H
    dH = np.mean(dH_mat, axis=0) * dy * dpsi
    H = sum(dH)
    CH = H / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)

    # Calculate observer frame lift and drag
    dL_observer = np.mean(dL_observer_mat, axis=0) * dy * psi
    dD_observer = np.mean(dD_observer_mat, axis=0) * dy * psi
    L_observer = sum(dL_observer)
    D_observer = sum(dD_observer)

    # Calculate rotor torque, and power
    dQ = np.mean(dQ_mat, axis=0) * dy * psi
    Q = sum(dQ)
    CQ = Q / (dens * np.pi * blade_rad**3 * (omega*blade_rad)**2)
    CP = CQ
    P = CP * dens * np.pi * blade_rad * (omega*blade_rad)**3

    return T, CT, H, CH, L_observer, D_observer, P, CP


def bemt_axial(propeller, pitch, omega, v_climb=0, alt=0, tip_loss=True, mach_corr=True, output='short'):
    # Calculate geopotential altitude
    alt_geop = (RAD_EARTH*alt)/(RAD_EARTH+alt)

    # Calculate atmospheric conditions
    temp = TEMP_SSL - A * alt_geop
    press = PRESS_SSL * (temp/TEMP_SSL)**(-Gs/(A*R))
    dens = press/(R*temp)
    spd_snd = np.sqrt(GAMMA * R * temp)
    kine_visc = 1.460 * 10**-5

    # Define blade geometry parameters. Pitch, chord, and r are all lists of the same length defining the blade
    # geometry at a specific span location r
    n_blades = propeller.n_blades
    blade_rad = propeller.radius
    twist = np.array(propeller.twist)

    chord = np.array(propeller.chord)
    dy = propeller.dy
    dr = propeller.dr
    r = np.array(propeller.r)
    y = np.array(propeller.y)
    n_elements = len(r)
    local_solidity = np.array(propeller.solidity)

    # Due to the possible camber of the airfoils along the span, we need to correct the local angle to include the zero
    # lift angle of attack. For positively cambered airfoils this will be a negative angle (all values of alpha0 will be
    # negative. Also find the lift curve slope along the span of the blade. Both quantities are calculated using an
    # approximate Reynolds number which calculates Re using only the in-plane portion of the freestream velocity.
    u_t = omega * r * blade_rad
    Re_approx = u_t * chord / kine_visc
    # Clalpha = propeller.get_Clalpha(Re_approx)
    # alpha0 = propeller.get_alpha0(Re_approx)
    Clalpha, alpha0 = propeller.get_Clalpha_alpha0(Re_approx)
    local_angle = pitch + twist

    # Define some other parameters for use in calculations
    v_tip = blade_rad * omega   # Blade tip speed
    lambda_c = v_climb/v_tip    # Climb inflow ratio

    # Now handle hover and vertical flight cases
    # First calculate inflow along span by using F = 1 to get initial value not including tip loss
    F = 1
    local_inflow, rel_inflow_angle, u_resultant = inflow.axial_flight(local_solidity, propeller, lambda_c, local_angle,
                                                                      alpha0, Clalpha, v_tip, v_climb, omega, r,
                                                                      blade_rad, F, spd_snd, mach_corr=mach_corr)
    # Now if tip_loss correction is desired, use the F = 1 solution as a starting guess to find the inflow
    if tip_loss:
        converged = np.array([False]*n_elements)
        while not all(converged):
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

    # Calculate Reynolds number along the span of the blade
    Re = u_resultant * chord / kine_visc

    # Now calculate the effective angle of attack at the blade stations
    eff_aoa = local_angle - rel_inflow_angle

    # Retrieve Cl and Cd values according to effective angle of attack along the blades. This will return NaN toward
    # the root
    Cl = np.nan_to_num(np.array(propeller.get_Cl(eff_aoa, Re)))
    Cd = np.nan_to_num(np.array(propeller.get_Cd(eff_aoa, Re)))

    # Calculate forces
    dL = 0.5*dens*u_resultant**2*chord*Cl*dy
    dD = 0.5*dens*u_resultant**2*chord*Cd*dy

    dFz = dL*np.cos(rel_inflow_angle) - dD*np.sin(rel_inflow_angle)
    dFx = dD*np.cos(rel_inflow_angle) + dL*np.sin(rel_inflow_angle)

    dT = n_blades * dFz
    dQ = n_blades * dFx * y
    dP = n_blades * dFx * omega * y

    T = sum(dT)
    Q = sum(dQ)
    P = sum(dP)

    CT = T / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)
    CP = P / (dens * np.pi * blade_rad**2 * (omega*blade_rad)**3)
    CQ = Q / (dens * np.pi * blade_rad**3 * (omega*blade_rad)**2)

    prop_CT = T / (dens * (omega/2/np.pi)**2 * (blade_rad*2)**4)
    prop_CP = P / (dens * (omega/2/np.pi)**3 * (blade_rad*2)**5)

    if output == 'short':
        return dT, P
    return dT, dP, P, Cd, Cl, u_resultant, chord, dL, local_inflow, rel_inflow_angle, eff_aoa, dFx, dFz, Re


# def bemt_axial_alt(propeller, pitch, omega, v_climb=0, alt=0, tip_loss=True, mach_corr=True, output='short'):
#     alt_geop = (RAD_EARTH*alt)/(RAD_EARTH+alt)
#
#     # Calculate atmospheric conditions
#     temp = TEMP_SSL - A * alt_geop
#     press = PRESS_SSL * (temp/TEMP_SSL)**(-Gs/(A*R))
#     dens = press/(R*temp)
#     spd_snd = np.sqrt(GAMMA * R * temp)
#     kine_visc = 1.460 * 10**-5
#
#     # Define blade geometry parameters. Pitch, chord, and r are all lists of the same length defining the blade
#     # geometry at a specific span location r
#     n_blades = propeller.n_blades
#     blade_rad = propeller.radius
#     twist = np.array(propeller.twist)
#
#     chord = np.array(propeller.chord)
#     dy = propeller.dy
#     dr = propeller.dr
#     r = np.array(propeller.r)
#     y = np.array(propeller.y)
#     n_elements = len(r)
#     local_solidity = np.array(propeller.solidity)
#
#     # Due to the possible camber of the airfoils along the span, we need to correct the local angle to include the zero
#     # lift angle of attack. For positively cambered airfoils this will be a negative angle (all values of alpha0 will be
#     # negative. Also find the lift curve slope along the span of the blade. Both quantities are calculated using an
#     # approximate Reynolds number which calculates Re using only the in-plane portion of the freestream velocity.
#     u_t = omega * r * blade_rad
#     Re_approx = u_t * chord / kine_visc
#     # Clalpha = propeller.get_Clalpha(Re_approx)
#     # alpha0 = propeller.get_alpha0(Re_approx)
#     Clalpha, alpha0 = propeller.get_Clalpha_alpha0(Re_approx)
#     local_angle = pitch + twist
#     Cd_approx = propeller.get_Cd(local_angle, Re_approx)
#
#     # Define some other parameters for use in calculations
#     v_tip = blade_rad * omega   # Blade tip speed
#     lambda_c = v_climb/v_tip    # Climb inflow ratio
#
#     # Now handle hover and vertical flight cases
#     # First calculate inflow along span by using F = 1 to get initial value not including tip loss
#     F = 1
#     local_inflow, rel_inflow_angle, u_resultant = inflow.axial_flight(local_solidity, propeller, lambda_c,
#                                                                                 local_angle, alpha0, Clalpha, v_tip,
#                                                                                 v_climb, omega, r, blade_rad, F,
#                                                                                 spd_snd, mach_corr=mach_corr)
#     # Now if tip_loss correction is desired, use the F = 1 solution as a starting guess to find the inflow
#     if tip_loss:
#         converged = np.array([False]*n_elements)
#         while not all(converged):
#             local_inflow_old = local_inflow
#             f_tip = n_blades/2. * ((1 - r)/(r * rel_inflow_angle))
#             f = f_tip
#             f[-1] = 0.0000000000001
#             F = (2/np.pi) * np.arccos(np.exp(-f))
#             try:
#                 local_inflow, rel_inflow_angle, u_resultant = inflow.axial_flight(local_solidity, propeller, lambda_c,
#                                                                                   local_angle, alpha0, Clalpha, v_tip,
#                                                                                   v_climb, omega, r, blade_rad, F,
#                                                                                   spd_snd, mach_corr=mach_corr)
#             except FloatingPointError:
#                 raise
#             converged = abs((local_inflow - local_inflow_old)/local_inflow) < 0.0005
#
#     # Calculate Reynolds number along the span of the blade
#     Re = u_resultant * chord / kine_visc
#
#     # Now calculate the effective angle of attack at the blade stations
#     eff_aoa = local_angle - rel_inflow_angle
#
#     # Retrieve Cl and Cd values according to effective angle of attack along the blades. This will return NaN toward
#     # the root
#     Cl = np.nan_to_num(np.array(propeller.get_Cl(eff_aoa, Re)))
#     Cd = np.nan_to_num(np.array(propeller.get_Cd(eff_aoa, Re)))
#
#     dL = 0.5*dens*u_resultant**2*chord*Cl*dy
#     dD = 0.5*dens*u_resultant**2*chord*Cd*dy
#
#     dFz = dL*np.cos(rel_inflow_angle) - dD*np.sin(rel_inflow_angle)
#     dFx = dL*np.sin(rel_inflow_angle) + dD*np.cos(rel_inflow_angle)
#
#     dT = n_blades * dFz
#     dP_i = n_blades * dL*np.sin(rel_inflow_angle) * omega * y
#     dP_o = n_blades * dD*np.cos(rel_inflow_angle) * omega * y
#     dP = n_blades * dFx * omega * y
#
#     dCT = 0.5 * n_blades * u_resultant**2 * chord * (Cl*np.cos(rel_inflow_angle) - Cd*np.sin(rel_inflow_angle)) * dy \
#           / (np.pi*blade_rad**2) / (omega*blade_rad)**2
#
#     dCP = n_blades*chord/np.pi/blade_rad/2*(rel_inflow_angle*Cl+Cd)*r**3*dy/blade_rad
#
#     dP_alt = dCP * (dens * np.pi * blade_rad**2 * (omega*blade_rad)**3)
#     dT_alt = dCT * (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)
#
#     P_alt = sum(dP_alt)
#     T_alt = sum(dT_alt)
#
#     T = sum(dT)
#     P = sum(dP)
#
#     # dCT = 0.5 * local_solidity * u_resultant**2 * (Cl*np.cos(rel_inflow_angle) - Cd*np.sin(rel_inflow_angle)) * dy / \
#     #       (omega**2 * blade_rad**3)
#     # dCP = 0.5 * local_solidity * u_resultant**2 * (Cl*np.sin(rel_inflow_angle) + Cd*np.cos(rel_inflow_angle)) * y * dy / \
#     #       (omega**2 * blade_rad**4)
#     #
#     # CT = sum(dCT)
#     # CP = sum(dCP)
#     #
#     # T = CT * (dens * np.pi * blade_rad**2 * (omega*blade_rad)**2)
#     # P = CP * (dens * np.pi * blade_rad**2 * (omega*blade_rad)**3)
#
#     return dT, dP, dP_i, dP_o





