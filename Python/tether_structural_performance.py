# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 12:26:00 2019

@author: thoma
"""

import numpy as np
import matplotlib.pyplot as plt

"""
inputs:
    - tether mass per m
    - operating altitude
    - altitude angle
    - tension on the ground
    - tension angle on the ground
"""

g_gravity = 9.81 #[m/s2]
m_tether_per_length = 0.12 #[kg/m]

w_tether_per_length = m_tether_per_length*g_gravity #[N]

h_operating_altitude = 4800 #[m]

theta_altitude_angle_deg = 20 #[deg]
theta_altitude_angle = theta_altitude_angle_deg*np.pi/180 #[rad]

#t_tension_tether_ground = 42000 #[N]
#theta_tension_angle_ground_deg = 10 #[deg]
#theta_tension_angle_ground = theta_tension_angle_ground_deg * np.pi/180 #[rad]
#H_horizontal_tension_force = t_tension_tether_ground*np.cos(theta_tension_angle_ground) #[N]

H_horizontal_tension_force = 25000 #[N]

L_horizontal_operating_distance = h_operating_altitude/np.tan(theta_altitude_angle) #[m]

L_star_tether_chord_length = h_operating_altitude/np.sin(theta_altitude_angle) #[m]

K1_sag_constant = np.arcsinh(w_tether_per_length*h_operating_altitude/(2*H_horizontal_tension_force*np.sinh(w_tether_per_length*L_horizontal_operating_distance/(2*H_horizontal_tension_force))))-w_tether_per_length*L_horizontal_operating_distance/(2*H_horizontal_tension_force)

K2_sag_constant = - H_horizontal_tension_force/w_tether_per_length*np.cosh(K1_sag_constant)

C_length_sag_cable = H_horizontal_tension_force/w_tether_per_length*(np.sinh(w_tether_per_length*L_horizontal_operating_distance/H_horizontal_tension_force+K1_sag_constant)-np.sinh(K1_sag_constant)) #[m]

x_cable_coordinate = [0] #[m]
y_cable_coordinate = [0] #[m]

for x in range(0,int(L_horizontal_operating_distance),1):
    x_cable_coordinate.append(x)
    y = H_horizontal_tension_force/w_tether_per_length*np.cosh(w_tether_per_length/H_horizontal_tension_force*x + K1_sag_constant) + K2_sag_constant
    y_cable_coordinate.append(y)
    
plt.plot(x_cable_coordinate,y_cable_coordinate)
print('total cable length', C_length_sag_cable, 'm')
print('straight cable length', L_star_tether_chord_length, 'm')
print('percentage difference = ', (C_length_sag_cable-L_star_tether_chord_length)/L_star_tether_chord_length*100, '%')

theta_tension_angle_ground = np.arctan(np.sinh(K1_sag_constant))
theta_tension_angle_glider = np.arctan(np.sinh(w_tether_per_length*L_horizontal_operating_distance/L_horizontal_operating_distance+K1_sag_constant))

print('glider rope angle = ', theta_tension_angle_glider*180/np.pi)
if theta_tension_angle_ground < 0:
    print('ERROR, tether touches ground, theta =', theta_tension_angle_ground*np.pi/180, 'degrees')
else:
    print('Phew, we are OK, theta=', theta_tension_angle_ground*180/np.pi, 'degrees')
