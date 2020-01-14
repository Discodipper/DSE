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
    - horizontal tension force
"""
g_gravity = 9.80665 #[m/s2]



###
def cable_sag_calculator(m_tether_per_length,h_operating_altitude,theta_altitude_angle_deg, H_horizontal_tension_force):
    
    w_tether_per_length = m_tether_per_length*g_gravity #[N] Weight per unit length of the tether
    
    theta_altitude_angle = theta_altitude_angle_deg*np.pi/180 #[rad] Angle between ground station and glider (straight cable)
    
    L_horizontal_operating_distance = h_operating_altitude/np.tan(theta_altitude_angle) #[m] Ground distance
    
    L_star_tether_chord_length = h_operating_altitude/np.sin(theta_altitude_angle) #[m] Chord length of cable, straight line cable length
    
    K1_sag_constant = np.arcsinh(w_tether_per_length*h_operating_altitude/(2*H_horizontal_tension_force*np.sinh(w_tether_per_length*L_horizontal_operating_distance/(2*H_horizontal_tension_force))))-w_tether_per_length*L_horizontal_operating_distance/(2*H_horizontal_tension_force)
    #print(K1_sag_constant)
    K2_sag_constant = - H_horizontal_tension_force/w_tether_per_length*np.cosh(K1_sag_constant)
    #print(K2_sag_constant)
    C_length_sag_cable = H_horizontal_tension_force/w_tether_per_length*(np.sinh(w_tether_per_length*L_horizontal_operating_distance/H_horizontal_tension_force+K1_sag_constant)-np.sinh(K1_sag_constant)) #[m]
    
    cable_angle_calculator(K1_sag_constant, w_tether_per_length, L_horizontal_operating_distance,H_horizontal_tension_force)
    
    cable_coordinates_calculator(L_horizontal_operating_distance, H_horizontal_tension_force, w_tether_per_length,K1_sag_constant,K2_sag_constant)
    
    #print('total cable length', C_length_sag_cable, 'm')
    #print('straight cable length', L_star_tether_chord_length, 'm')
    #print('percentage difference = ', (C_length_sag_cable-L_star_tether_chord_length)/L_star_tether_chord_length*100, '%')
    
    return C_length_sag_cable;    
###        


###
def cable_angle_calculator(K1_sag_constant, w_tether_per_length, L_horizontal_operating_distance,H_horizontal_tension_force):
    
    theta_tension_angle_ground = np.arctan(np.sinh(K1_sag_constant))
    theta_tension_angle_glider = np.arctan(np.sinh(w_tether_per_length*L_horizontal_operating_distance/H_horizontal_tension_force+K1_sag_constant))
    
    #print('glider rope angle = ', theta_tension_angle_glider*180/np.pi, 'degrees')
    if np.all(theta_tension_angle_ground) < 0:
        #print('ERROR, tether touches ground, theta =', theta_tension_angle_ground*np.pi/180, 'degrees')
    # else:
    #     #print('Phew, we are OK, theta=', theta_tension_angle_ground*180/np.pi, 'degrees')
        return(theta_tension_angle_glider, theta_tension_angle_ground)
###

###
def cable_coordinates_calculator(L_horizontal_operating_distance, H_horizontal_tension_force, w_tether_per_length,K1_sag_constant,K2_sag_constant):
    x_cable_coordinate = [] #[m]
    y_cable_coordinate = [] #[m]
    for x in range(0,int(L_horizontal_operating_distance),1):
        x_cable_coordinate.append(x)
        y = H_horizontal_tension_force/w_tether_per_length*np.cosh(w_tether_per_length/H_horizontal_tension_force*x + K1_sag_constant) + K2_sag_constant
        y_cable_coordinate.append(y)
        
    plt.plot(x_cable_coordinate,y_cable_coordinate)
    plt.ylim(0, L_horizontal_operating_distance+0.1*L_horizontal_operating_distance)
    plt.xlim(0, L_horizontal_operating_distance+0.1*L_horizontal_operating_distance)
    plt.gca().set_aspect('equal', adjustable='box')
    
    return x_cable_coordinate, y_cable_coordinate
###

###
def cable_dimensions_calculator(tension_force_cable, ultimate_tensile_strength, density_tether):   
    safety_factor = 5
    A_crosssectional_area_cable = safety_factor*tension_force_cable/ultimate_tensile_strength #[m^2]
    r_radius_cable = np.sqrt(A_crosssectional_area_cable/np.pi) #[m]
    
    m_tether_per_length = density_tether*A_crosssectional_area_cable
    return(m_tether_per_length, r_radius_cable)
###


#inputs:
m_tether_per_length, r_radius_cable = cable_dimensions_calculator(42000,3000*10**6, 1000) #[kg/m], [m]
d_diameter_cable_mm = 2*r_radius_cable*1000 #[mm]
h_operating_altitude = 2500 #[m]
theta_altitude_angle_deg = 25 #[deg]
H_horizontal_tension_force = 33000 #[N]

C_length_sag_cable = cable_sag_calculator(m_tether_per_length,h_operating_altitude,theta_altitude_angle_deg, H_horizontal_tension_force)