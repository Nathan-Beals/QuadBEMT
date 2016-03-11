import numpy as np
import matplotlib.pyplot as plt
import propeller
import bemt
from unit_conversion import rpm2radps, radps2rpm, rad2deg, deg2rad, in2m


n_elements = 100
n_azi_elements = 100
n_blades = 2
blade_radius = in2m(8.90)
root_cutout = 0.1 * blade_radius
chord = np.array([0.127, 0.135, 0.158, 0.178, 0.195, 0.209, 0.219, 0.225, 0.227, 0.226, 0.221, 0.212, 0.199, 0.182,
                  0.161, 0.135, 0.097, 0.058]) * blade_radius
twist = np.array([27.54, 25.28, 26.08, 25.47, 24.07, 22.18, 20.00, 18.18, 16.38, 14.83, 13.63, 12.56, 11.56, 10.65,
                  9.68, 8.51, 6.72, 4.89]) * 2 * np.pi / 360
pitch = 0
dy = float(blade_radius-root_cutout)/n_elements
dr = float(1)/n_elements
y = root_cutout + dy*np.arange(1, n_elements+1)
r = y/blade_radius
Clalpha = 2 * np.pi

slow_fly_9 = propeller.Propeller(twist, chord, blade_radius, n_blades, r, y, dr, dy, Clalpha)

# Let's say we want to calculate the power required to fly forward with a required thrust of T = 3 N at a forward speed
# of 4 m/s.
L_req = 3
v_inf = 4

# Set up iteration parameters
omega0_RPM = 1000
omega_max_RPM = 7000
omega_inc_RPM = 1
alpha0_deg = 0.01
alpha_max_deg = 10
alpha_inc_deg = 0.01

omega0 = rpm2radps(omega0_RPM)
omega_max = rpm2radps(omega_max_RPM)
omega_inc = rpm2radps(omega_inc_RPM)
alpha0 = deg2rad(alpha0_deg)
alpha_max = deg2rad(alpha_max_deg)
alpha_inc = deg2rad(alpha_inc_deg)

converged = False
while not converged:
    for omega in xrange(omega0, omega_max, omega_inc):
        for alpha in xrange(alpha0, alpha_max, alpha_inc):
            T, CT, H, CH, L_observer, D_observer, P, CP = bemt.bemt_forward_flight(slow_fly_9, pitch, omega, alpha,
                                                                                   T_req, v_inf, n_azi_elements)
            alpha_requirement = np.arctan(D_observer/)