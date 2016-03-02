import numpy as np
import inflow

# Define some constants
RAD_EARTH = 6371000    # Radius of Earth in meters
GAMMA = 1.4    # Heat capacity ratio for air
R = 287    # For air R = 287 J/kg K
A = -0.0065    # lapse rate in /m
Gs = -9.8    # Acceleration due to gravity in m/s**2
TEMP_SSL = 288.16    # Standard sea level temp in K
PRESS_SSL = 1.01325 * 10**5    # Standard sea level pressure in N/m**2


def bemt(propeller, pitch, omega, alpha, v_climb=0, v_inf=0, alt=0, tip_loss=True, CT_target=0, mach_corr=True,
         ff_inflow='uniform'):

    if v_climb == 0 and v_inf == 0:
        mode = 'hover'
    elif v_climb != 0 and v_inf == 0:
        mode = 'vertical'
    else:
        mode = 'forward'

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
    n_azi_elements = 60
    local_solidity = np.array(propeller.solidity)

    # Define some other parameters for use in calculations
    v_tip = blade_rad * omega   # Blade tip speed
    lambda_c = v_climb/v_tip    # Climb inflow ratio
    mu = v_inf*np.cos(alpha)/(omega*blade_rad)  # Advance ratio

    # Now handle hover and vertical flight cases
    if mode in ['hover', 'vertical']:
        # First calculate inflow along span by using F = 1 to get initial value not including tip loss
        F = 1
        local_inflow, rel_inflow_angle, u_resultant = inflow.axial_flight(local_solidity, propeller, lambda_c,
                                                                          local_angle, v_tip, v_climb, omega, r,
                                                                          blade_rad, F, spd_snd, mach_corr=mach_corr)
        # Now if tip_loss correction is desired, use the F = 1 solution as a starting guess to find the inflow
        if tip_loss:
            converged = np.array([False]*n_elements)
            c = 0
            while not all(converged):
                local_inflow_old = local_inflow
                f_tip = n_blades/2. * ((1 - r)/(r * rel_inflow_angle))
                f = f_tip
                F = (2/np.pi) * np.arccos(np.exp(-f))
                local_inflow, rel_inflow_angle, u_resultant = inflow.axial_flight(local_solidity, propeller, lambda_c,
                                                                                  local_angle, v_tip, v_climb, omega, r,
                                                                                  blade_rad, F, spd_snd,
                                                                                  mach_corr=mach_corr)
                converged = abs((local_inflow - local_inflow_old)/local_inflow) < 0.0005
                c += 1

        # Now calculate the effective angle of attack at the blade stations
        eff_aoa = local_angle - rel_inflow_angle

        # Retrieve Cl and Cd values according to effective angle of attack along the blades
        Cl = np.array(propeller.get_Cl(eff_aoa))
        Cd = np.array(propeller.get_Cd(eff_aoa))
        dCl = Cl

        # Calculate forces
        dL = 0.5*dens*u_resultant**2*chord*Cl*dy
        dD = 0.5*dens*u_resultant**2*chord*Cd*dy

        dT = n_blades*(dL*np.cos(rel_inflow_angle)-dD*np.sin(rel_inflow_angle))
        dQ = n_blades*(dD*np.cos(rel_inflow_angle)+dL*np.sin(rel_inflow_angle))
        T = sum(dT)
        Q = sum(dQ)

        dPi = n_blades*dL*np.sin(rel_inflow_angle)*y*omega
        dPp = n_blades*(dD*np.cos(rel_inflow_angle))*y*omega
        dPtot = dPi + dPp
        Pi = sum(dPi)
        Pp = sum(dPp)
        Ptot = sum(dPtot)

        dCt = dT/(dens*np.pi*blade_rad**2*v_tip**2)
        dCp = dPtot/(dens*np.pi*blade_rad**2*v_tip**3)
        CT = sum(dCt)
        CP = sum(dCp)
    else:
        # Now do the forward flight case. Since the velocity components normal and in plane with the TPP are now a
        # function of the blade azimuth angle, psi, we need to average forces over the entire range of psi.
        psi = np.linspace(0, 2*np.pi, n_azi_elements)
        dpsi = 2 * np.pi * y / n_azi_elements   # size of d_psi for each annulus

        if ff_inflow == 'uniform':
            local_inflow = inflow.uniform_ff(CT_target, alpha, mu, n_elements)

        dT_azi_total = np.zeros(n_elements)
        dPi_azi_total = np.zeros(n_elements)
        dCpo_azi_total = np.zeros(n_elements)
        dCl_azi_total = np.zeros(n_elements)
        for azi_ang in psi:
            u_p = local_inflow * omega * blade_rad * np.sin(azi_ang)
            u_t = omega*y + mu*omega*blade_rad*np.sin(azi_ang)
            u_r = mu * omega * blade_rad * np.cos(azi_ang)
            local_mach = u_t / spd_snd
            rel_inflow_angle = np.arctan(u_p/u_t)
            eff_aoa = local_angle - rel_inflow_angle
            Cl, Cd = propeller.aero_coeffs(eff_aoa)
            # Mach number correction to the lift coefficient
            if mach_corr:
                Cl /= np.sqrt(1 - local_mach**2)
            dCl_azi_total += Cl
            # Calculate thrust
            dT_azi_total += 0.5 * dens * (u_t**2+u_p**2) * chord * \
                            (Cl*np.cos(rel_inflow_angle)-Cd*np.sin(rel_inflow_angle)) * dpsi
            # Calculate induced power
            dPi_azi_total += 0.5 * dens * (u_t**2+u_p**2) * chord * \
                             (Cd*np.cos(rel_inflow_angle)+Cl*np.sin(rel_inflow_angle)) * omega * dpsi
            # Calculate profile power coefficient
            coslamb = abs(u_t)/(u_t**2+u_r**2)**0.5
            Cd3d = Cd * (alpha*coslamb)/coslamb
            D = 0.5 * u_t * abs(u_t) * chord * Cd3d
            Fr = D * (u_r/u_t)
            Fx = D * np.cos(rel_inflow_angle)
            dCpo_azi_total += local_solidity * (u_t*Fx/chord + u_r*Fr/chord) * dpsi

        # Calculate forces
        dCl = dCl_azi_total / n_azi_elements
        dT = n_blades * dy * dT_azi_total / n_azi_elements
        dCt = dT/(dens * np.pi * blade_rad**2 * v_tip**2)
        T = sum(dT)
        CT = sum(dCt)

        dPi = n_blades * y * dy * dPi_azi_total / n_azi_elements
        dCpi = dPi/(dens * np.pi * blade_rad**2 * v_tip**3)
        dCpo = n_blades * dr * dCpo_azi_total / n_azi_elements
        dPo = dens * np.pi * blade_rad**2 * v_tip**3 * dCpo
        Pi = sum(dPi)
        CPi = sum(dCpi)
        Po = sum(dPo)
        CPo = sum(dCpo)
        CP = CPi + CPo
        Ptot = Pi + Po

    return CT, CP, dCt, Ptot, local_inflow, rel_inflow_angle, dCl




