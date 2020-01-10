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
#-----------------------------------------------------------------------------

######## INPUTS #########

g_gravity = 9.80665 #[m/s2]


#Banana fibre: ultimate_strength = 790*10**6; material_density = 1300 # kg/m3
lift_force = 400000
ultimate_strength = 790*10**6 # yield or ultimate stress
material_density = 1300 # kg/m3
#E_youngs_modulus = 26*10**9

h_operating_altitude = 2000 #[m]
theta1_ideal = 45.
H_horizontal_tension_force = lift_force*np.cos(theta1_ideal*np.pi/180) #[N]

theta = []
theta1 = []
theta0 = []
sagged_cable_length = []
theta1comp = []
sag_list = []

#-----------------------------------------------------------------------------

########## DEFINITIONS ##########

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
    
    theta_tension_angle_glider, theta_tension_angle_ground = cable_angle_calculator(K1_sag_constant, w_tether_per_length, L_horizontal_operating_distance,H_horizontal_tension_force)
    sag = cable_coordinates_calculator(L_horizontal_operating_distance, H_horizontal_tension_force, w_tether_per_length,K1_sag_constant,K2_sag_constant, theta_altitude_angle)
    
    #print('total cable length', C_length_sag_cable, 'm')
    #print('straight cable length', L_star_tether_chord_length, 'm')
    #print('percentage difference = ', (C_length_sag_cable-L_star_tether_chord_length)/L_star_tether_chord_length*100, '%')
    
    return C_length_sag_cable, theta_tension_angle_glider, theta_tension_angle_ground, sag;    
###        


###
def cable_angle_calculator(K1_sag_constant, w_tether_per_length, L_horizontal_operating_distance,H_horizontal_tension_force):
    
    theta_tension_angle_ground = np.arctan(np.sinh(K1_sag_constant))
    theta_tension_angle_glider = np.arctan(np.sinh(w_tether_per_length*L_horizontal_operating_distance/H_horizontal_tension_force+K1_sag_constant))
    
    #print('glider rope angle = ', theta_tension_angle_glider*180/np.pi, 'degrees')
    #if theta_tension_angle_ground < 0:
    #    print('ERROR, tether touches ground, theta =', theta_tension_angle_ground*np.pi/180, 'degrees')
    #else:
    #    print('Phew, we are OK, theta=', theta_tension_angle_ground*180/np.pi, 'degrees')
    return theta_tension_angle_glider, theta_tension_angle_ground;
###

###
def cable_coordinates_calculator(L_horizontal_operating_distance, H_horizontal_tension_force, w_tether_per_length,K1_sag_constant,K2_sag_constant,theta_altitude_angle):
    x_cable_coordinate = [] #[m]
    y_cable_coordinate = [] #[m]
    y_chord_coordinate = [] #[m]
    
    for x in range(0,int(L_horizontal_operating_distance),1):
        x_cable_coordinate.append(x)
        y = H_horizontal_tension_force/w_tether_per_length*np.cosh(w_tether_per_length/H_horizontal_tension_force*x + K1_sag_constant) + K2_sag_constant
        y_cable_coordinate.append(y)
        y_chord_coordinate.append(x*np.tan(theta_altitude_angle))
        
    plt.plot(x_cable_coordinate,y_cable_coordinate)
    plt.plot(x_cable_coordinate,y_chord_coordinate)
    plt.ylim(0, L_horizontal_operating_distance+0.1*L_horizontal_operating_distance)
    plt.xlim(0, L_horizontal_operating_distance+0.1*L_horizontal_operating_distance)
    plt.gca().set_aspect('equal', adjustable='box')
    
    sag = max(np.subtract(y_chord_coordinate,y_cable_coordinate))
    return sag, x_cable_coordinate, y_cable_coordinate
###

###
def cable_dimensions_calculator(tension_force_cable, ultimate_tensile_strength, density_tether):   
    safety_factor = 5
#    A_crosssectional_area_cable = safety_factor*tension_force_cable/ultimate_tensile_strength #[m^2]
#    r_radius_cable = np.sqrt(A_crosssectional_area_cable/np.pi) #[m]
    pure_fibre_area = safety_factor*tension_force_cable/ultimate_tensile_strength #[m^2]
    m_tether_per_length = density_tether*pure_fibre_area*(1 + pure_fibre_area/0.0005037987746501848) #[kg/m]
    rope_area = m_tether_per_length/density_tether #[m^2]
#    m_tether_per_length = density_tether*A_crosssectional_area_cable #[kg/m]
    r_radius_cable = np.sqrt(rope_area/np.pi)
    return m_tether_per_length, r_radius_cable;
###

###

#-----------------------------------------------------------------------------

######## COMPUTATIONS #########

# Calculate properties of tether: 
m_tether_per_length, r_radius_cable = cable_dimensions_calculator(lift_force,ultimate_strength, material_density) #[kg/m], [m]
d_diameter_cable_mm = 2*r_radius_cable*1000 #[mm]

# Calculate sag properties of tether:

# Call functions in for loop for different angles between glider and ground (theta)
# for the closest match between the angle of lift (theta1, theta1_ideal).
for theta_altitude_angle_deg in np.arange(15.,theta1_ideal,0.2):
    C_length_sag_cable, theta_tension_angle_glider, theta_tension_angle_ground, sag = cable_sag_calculator(m_tether_per_length,h_operating_altitude,theta_altitude_angle_deg, H_horizontal_tension_force)
    theta.append(theta_altitude_angle_deg)
    theta1.append(theta_tension_angle_glider*180/np.pi)
    theta1comp.append(abs(theta_tension_angle_glider*180/np.pi-theta1_ideal))
    theta0.append(theta_tension_angle_ground*180/np.pi)
    sagged_cable_length.append(C_length_sag_cable)
    sag_list.append(sag[0])

    
# Following are the final angles:    
operating_angle_final = theta[theta1comp.index(min(theta1comp))]
theta0_final = theta0[theta1comp.index(min(theta1comp))]
sag_final = sag_list[theta1comp.index(min(theta1comp))]

if theta0_final < 0:
    print('Whoops, the cable will hit the ground, theta0 = ', theta0_final)
theta_tension_angle_glider_final = theta1[theta1comp.index(min(theta1comp))]


# Plot the cable sag:
plt.clf()
C_length_sag_cable_final = cable_sag_calculator(m_tether_per_length,h_operating_altitude,operating_angle_final, H_horizontal_tension_force)   
print('total (sagged) cable length =', C_length_sag_cable_final[0])
total_mass_cable = m_tether_per_length*C_length_sag_cable_final[0]   #kg
print('total mass tether =',total_mass_cable)
total_tension_ground = H_horizontal_tension_force/np.cos(theta0_final*np.pi/180)
print('total tesion force at ground =', total_tension_ground) 
print('maximum sag = ', sag_final)   
    
    
    