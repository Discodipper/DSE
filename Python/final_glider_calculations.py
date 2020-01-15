#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 12:22:02 2019

@author: pranav
"""

from math import *
import matplotlib.pyplot as plt
import numpy as np
import time
from ISA_calculator import isa
from flight_dynamics_tether_weight import cable_dimensions_calculator
from flight_dynamics_tether_weight import cable_sag_calculator
from flight_dynamics_tether_weight import cable_angle_calculator
from flightpath import reduced_lift_by_roll
from flightpath import flight_radius
import scipy.linalg as la
start = time.time()

reel_factor_array = np.arange(0.3, 4, 0.1) #m/s
altitude_array = np.arange(3000, 3100, 100)

"""theta is polar angle, phi is azimuth angle, beta is elevation angle (operation angle)."""
azimuth_angle = 0*np.pi/180 #rad
hi = 0*pi/180 #rad
g_gravity = 9.80665
gliders_per_tether = 2

lift_coefficient = 0.995 + 0.092285#1.005 + 0.092285#1.005 + 0.0464#0.067+0.0546
drag_coefficient = 0.03388 + 0.0509 + 0.006473#0.034+0.092285
glide_ratio = lift_coefficient/drag_coefficient
wing_area = 60 #m^2


def reel_speed(f, V_wind):
    V_reel = V_wind * f
    return(V_reel)

def tangential_velocity_factor(theta, phi, hi, glide_ratio, f):
    a = np.cos(theta) * np.cos(phi) * np.cos(hi) - np.sin(phi)*np.sin(hi)
    b = np.sin(theta) * np.cos(phi)
    tvf = a + np.sqrt(a**2 + b**2 -1 + (glide_ratio**2)*(b - f)**2)
    return(tvf, a, b)

def azimuth_constraint(polar_angle, phi, glide_ratio, f):
    x = np.sin(polar_angle) * np.cos(phi)
    y = (np.sqrt(1 + glide_ratio**2 * (1-f**2)) + f*glide_ratio**2)/(1 + glide_ratio**2)
    diff = x - y
    return(diff)

def spherical_to_cartesian(x, theta, phi):
    #A = la.inv(np.matrix([[np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)], [np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)], [-np.sin(phi), np.cos(phi), 0]]))
    A = np.matrix([[np.sin(theta) * np.cos(phi), np.cos(theta) * np.cos(phi), -np.sin(phi)], [np.sin(theta) * np.sin(phi), np.cos(theta) * np.sin(phi), np.cos(phi)], [np.cos(theta), -np.sin(theta), 0]])
    B = np.dot(A, x)
    return(B)

def resultant_force_coefficient(lift_coefficient, drag_coefficient):
    C_R = np.sqrt(lift_coefficient**2 + drag_coefficient**2)
    return(C_R)

def max_elevation_angle(glide_ratio, reel_factor):
    beta_max = np.arccos((reel_factor*glide_ratio**2 + np.sqrt(1 + glide_ratio**2 * (1-reel_factor**2)))/(1 + glide_ratio**2))
    return(beta_max)

def apparent_wind_speed_values(windspeed, reelspeed, Lambda, polar_angle, azimuth_angle, hi):
    kite_speed_tangential = Lambda * windspeed
    wind_component = np.matrix([[np.sin(polar_angle) * np.cos(azimuth_angle)], [np.cos(polar_angle) * np.cos(azimuth_angle)], [-np.sin(azimuth_angle)]]) * windspeed
    kite_radial_component = np.matrix([[1], [0], [0]]) * reelspeed
    kite_tangential_component = np.matrix([[0], [np.cos(hi)], [np.sin(hi)]]) * kite_speed_tangential
    V_a_spherical = wind_component - kite_radial_component - kite_tangential_component
    V_a_cartesian = spherical_to_cartesian(V_a_spherical, polar_angle, azimuth_angle)
    return(V_a_spherical, V_a_cartesian)

def apparent_wind_speed_magnitude(apparent_wind_speed_cartesian):
    V_a_magnitude = la.norm(np.array(apparent_wind_speed_cartesian))
    return(V_a_magnitude)
    
def tether_force_max(V_a, air_density, lift_coefficient, drag_coefficient, wing_area, elevation_angle, gliders_per_tether):
    F_t_max = (0.5*air_density*resultant_force_coefficient(lift_coefficient, drag_coefficient)* V_a**2 * wing_area) * gliders_per_tether
    F_t_max_horizontal = F_t_max * np.sin(elevation_angle)
    return(F_t_max, F_t_max_horizontal)

def tether_mass_guess(tether_diameter_guess, tether_density, altitude, elevation_angle, F_t_max_horizontal):
    tethermass_per_length = pi*(tether_diameter_guess/2)**2 * tether_density
    tether_length = cable_sag_calculator(tethermass_per_length, altitude, elevation_angle, F_t_max_horizontal)[0]
    K1_sag_constant = cable_sag_calculator(tethermass_per_length, altitude, elevation_angle, F_t_max_horizontal)[1]
    L_horizontal_operating_distance = cable_sag_calculator(tethermass_per_length, altitude, elevation_angle, F_t_max_horizontal)[2]
    tethermass_guess = tethermass_per_length * tether_length
    return(tethermass_guess, tether_length, K1_sag_constant, L_horizontal_operating_distance)

def force_z_direction(V_a_z, air_density, wing_area, drag_coefficient):
    #V_a_z = apparent_wind_speed_cartesian.item(2)
    F_z = 0.5 * air_density * V_a_z**2 * wing_area * drag_coefficient
    return(F_z)

def net_tether_force(V_a, air_density, lift_coefficient, drag_coefficient, wing_area, elevation_angle, V_a_z, glidermass, g_gravity, F_z, tethermass, gliders_per_tether):
    F_t_net = tether_force_max(V_a, air_density, lift_coefficient, drag_coefficient, wing_area, elevation_angle, gliders_per_tether)[0] - (glidermass*g_gravity*np.cos(elevation_angle)) - (tethermass*g_gravity) - (F_z*np.cos(elevation_angle))
    #F_t_net_z = F_t_net * np.sin(elevation_angle)
    return(F_t_net)#, F_t_net_z)

def net_tether_force_z(F_t_net, elevation_angle):
    F_t_net_z = F_t_net * np.sin(elevation_angle)
    return(F_t_net_z)

def tether_diameter_new(F_t_net_z, ultimate_tensile_strength, density_tether):
    tethermass_per_length_new = cable_dimensions_calculator(F_t_net_z, ultimate_tensile_strength, density_tether)[0]
    tether_diameter_new = cable_dimensions_calculator(F_t_net_z, ultimate_tensile_strength, density_tether)[1] * 2
    return(tether_diameter_new)    

def polarangle(elevation_angle):
    #Make sure the input elevation angle is in radians
    THETA = (pi/2) - elevation_angle
    return(THETA)

def generator_power(tether_force, reelspeed):
    P_gen = tether_force * reelspeed
    #Still to be completed
    return(P_gen)

def reel_out_effective_angle_of_attack(V_a):
    V_a_x = V_a.item(0)
    V_a_y = abs(V_a.item(1))
    aoa = np.arctan2(V_a_x/V_a_y)
    return(aoa)

def added_velocity_z_direction(hi, elevation_angle, radius_flight_path, g_gravity):
    if np.sin(hi) == 0:
        added_velocity = 0
    else:
        added_velocity = np.sqrt(2 * g_gravity * radius_flight_path * abs(np.sin(hi)) * np.cos(elevation_angle)) * -(np.sin(hi) / abs(np.sin(hi)))
    return(added_velocity)

apparent_wind_speed_lst = []
apparent_wind_speed_magnitude_lst = []
reelspeed_lst = []
altitude_lst = []
operation_angle_lst = []
azimuth_angle_lst = []
Lambda_lst = []
a_lst = []
b_lst = []
b_f_lst = []
wind_component_lst = []
kite_radial_component_lst = []
kite_tangential_component_lst = []
V_w_lst = []
generator_power_lst = []
tether_force_lst = []
reel_factor_lst = []
lift_lst = []


max_apparent_wind_speed_magnitude = 0
corresponding_apparent_wind_speed = 0
corresponding_altitude = 0
corresponding_beta = 0
corresponding_reel_speed = 0
corresponding_wind_speed = 0
corresponding_tether_force = 0
corresponding_generator_power = 0
corresponding_tether_mass = 0
corresponding_generator_power = 0
corresponding_total_mass = 0
corresponding_Lambda = 0
corresponding_lift = 0
corresponding_flight_radius = 0
corresponding_F_z = 0
#generator_angle = 0

d_tether = 0
l_tether = 0
corresponding_reel_factor = 0

tether_guess_lst = []
tether_needed_lst = []
tether_diff_lst = []
count_lst = []
flight_radius_lst = []

for reel_factor in reel_factor_array:
    for altitude in altitude_array:
        
        temperature, pressure, air_density, windspeed = isa(altitude)
        reelspeed = reel_speed(reel_factor, windspeed)
        
        beta_max = max_elevation_angle(glide_ratio, reel_factor) * 180/pi
        if beta_max >= 60:
            beta = (np.linspace(20, beta_max, 20))*pi/180 #rad
            for operation_angle in beta:
                polar_angle = polarangle(operation_angle) #rad
                    
                Lambda = tangential_velocity_factor(polar_angle, azimuth_angle, hi, glide_ratio, reel_factor)[0]
                
                apparent_wind_speed_spherical, apparent_wind_speed_cartesian = apparent_wind_speed_values(windspeed, reelspeed, Lambda, polar_angle, azimuth_angle, hi)
                magnitude_apparent_wind_speed = apparent_wind_speed_magnitude(apparent_wind_speed_cartesian)
                V_a_z = apparent_wind_speed_cartesian.item(2)
                tether_diameter_initial_guess = 0.07 #m
                tether_diameter_difference = 1
                ultimate_tensile_strength = 3.6*10**9 # yield or ultimate stress
                tether_density = 975 # kg/m3
                glidermass = 3080*gliders_per_tether #kg
                roll_angle = 30 * pi/180 #operation_angle
                count = 0
                radius_flight_path = 2000 #m
                total_mass = 0 #kg
                lift = 0
                tether_mass_final = 0
                F_z = 0
                reduced_lift = 0
                tether_length = 0
                cable_angle = 0
               
                while tether_diameter_difference > 0.005:
                    total_tether_force, total_tether_force_horizontal = (tether_force_max(magnitude_apparent_wind_speed, air_density, lift_coefficient, drag_coefficient, wing_area, operation_angle, gliders_per_tether))
                    #print("total_tether_force = ", total_tether_force)
                    F_z = force_z_direction(V_a_z, air_density, wing_area, drag_coefficient)
                    tether_mass_guess_value = tether_mass_guess(tether_diameter_initial_guess, tether_density, altitude, operation_angle*180/np.pi, total_tether_force_horizontal)[0]
                    #print("tether mass = ", tether_mass_guess_value)
                    tether_force_net = net_tether_force(magnitude_apparent_wind_speed, air_density, lift_coefficient, drag_coefficient, wing_area, operation_angle, V_a_z, glidermass, g_gravity, F_z, tether_mass_guess_value, gliders_per_tether)
                    #print("net_tether_force = ", tether_force_net)
                    tether_force_net_z = net_tether_force_z(tether_force_net, operation_angle)
                    
                    #reduced_force_by_roll = reduced_lift_by_roll(corresponding_tether_force, roll_angle)[0]
                    reduced_force_by_roll = reduced_lift_by_roll(tether_force_net, roll_angle)#, wing_area, air_density, lift_coefficient, radius_flight_path)
                    tether_force_net = tether_force_net - reduced_force_by_roll
                    final_generator_power = reduced_force_by_roll * reelspeed
                    
                    
                    tether_diameter_needed = tether_diameter_new(tether_force_net_z, ultimate_tensile_strength, tether_density)
                    tether_diameter_difference_absolute = abs(tether_diameter_initial_guess - tether_diameter_needed)
                    tether_diameter_difference = tether_diameter_difference_absolute
                    tether_mass_final = tether_mass_guess(tether_diameter_initial_guess, tether_density, altitude, operation_angle*180/np.pi, total_tether_force_horizontal)[0]
                    tether_length__final = tether_mass_guess(tether_diameter_initial_guess, tether_density, altitude, operation_angle*180/np.pi, total_tether_force_horizontal)[1]
                    # K1 = tether_mass_guess(tether_diameter_initial_guess, tether_density, altitude, operation_angle*180/np.pi, total_tether_force_horizontal)[2]
                    # L_h = tether_mass_guess(tether_diameter_initial_guess, tether_density, altitude, operation_angle*180/np.pi, total_tether_force_horizontal)[3]
                    # angle_tether = cable_angle_calculator(K1, tether_mass_final*g_gravity, L_h, tether_force_net*np.sin(operation_angle))
                    
                    total_mass = tether_mass_final + glidermass*gliders_per_tether
                    total_mass_per_glider = glidermass + tether_mass_final/gliders_per_tether
                    kite_speed = np.sqrt(apparent_wind_speed_cartesian.item(1)**2 + apparent_wind_speed_cartesian.item(2)**2)
                    final_flight_radius = flight_radius(roll_angle, total_tether_force/gliders_per_tether, total_mass_per_glider, magnitude_apparent_wind_speed)
                    
                    
                    #print("d_guess = ", tether_diameter_initial_guess)
                    tether_diameter_initial_guess = tether_diameter_needed
                    tether_length = tether_length__final
                    #print("diff = ", tether_diameter_difference_absolute)
                    #print("d_needed = ", tether_diameter_needed)
                    count = count + 1
                    radius_flight_path = final_flight_radius
                    lift = total_tether_force
                    reduced_lift = reduced_force_by_roll
                    #cable_angle = angle_tether
                    
                    tether_guess_lst.append(tether_diameter_initial_guess)
                    tether_needed_lst.append(tether_diameter_needed)
                    tether_diff_lst.append(tether_diameter_difference_absolute)
                    count_lst.append(count)
                if np.isnan(tether_diameter_difference)==True:
                        break
                
                # plt.scatter(count_lst, tether_guess_lst, tether_needed_lst, tether_diff_lst)
                # plt.show()
                
                apparent_wind_speed_lst.append(apparent_wind_speed_cartesian)
                apparent_wind_speed_magnitude_lst.append(apparent_wind_speed_magnitude)
                reelspeed_lst.append(reelspeed)
                altitude_lst.append(altitude)
                operation_angle_lst.append(operation_angle*180/pi)
                azimuth_angle_lst.append(azimuth_angle)
                Lambda_lst.append(Lambda)
                V_w_lst.append(windspeed)
                tether_force_lst.append(tether_force_net)
                generator_power_lst.append(generator_power(tether_force_net, reelspeed))
                reel_factor_lst.append(reel_factor)
                lift_lst.append(lift)
                flight_radius_lst.append(radius_flight_path)
                
                # if magnitude_apparent_wind_speed > max_apparent_wind_speed_magnitude:
                if generator_power(tether_force_net, reelspeed) > corresponding_generator_power:
                    max_apparent_wind_speed_magnitude = magnitude_apparent_wind_speed
                    corresponding_altitude = altitude
                    corresponding_reel_speed = reelspeed
                    corresponding_wind_speed = windspeed
                    corresponding_reel_factor = reel_factor
                    corresponding_generator_power = generator_power(tether_force_net, reelspeed)
                    corresponding_beta = operation_angle
                    corresponding_apparent_wind_speed = apparent_wind_speed_cartesian
                    corresponding_apparent_wind_speed_spherical = apparent_wind_speed_spherical
                    corresponding_tether_force = tether_force_net
                    corresponding_tether_mass = tether_mass_final
                    # corresponding_generator_power = final_generator_power
                    corresponding_total_mass = total_mass + glidermass
                    corresponding_flight_radius =  final_flight_radius
                    corresponding_reduced_lift_roll = reduced_lift 
                    corresponding_total_mass = total_mass
                    corresponding_Lambda = Lambda
                    corresponding_lift = lift
                    corresponding_F_z = F_z
                    d_tether = tether_diameter_initial_guess
                    l_tether = tether_length
                    #generator_angle = cable_angle
                

            





#print("V_a = ", apparent_wind_speed_lst[7149])
print("Maximum magnitude of V_a = ", max_apparent_wind_speed_magnitude)
print("V_a vector = ", corresponding_apparent_wind_speed)
print("alitutude = ", corresponding_altitude)
print("V_reel = ", corresponding_reel_speed)
print("f = ", corresponding_reel_factor)
print("wind speed = ", corresponding_wind_speed)
print("Elevation angle = ", corresponding_beta*180/pi)
print("Tether force = ", corresponding_tether_force)
print("Final generator power", corresponding_generator_power)
print("Reduced lift due to roll", corresponding_reduced_lift_roll)
print("Flight radius", corresponding_flight_radius)
print("Total mass = :", corresponding_total_mass)
print("Tether mass", corresponding_tether_mass)
print("Lambda = ", corresponding_Lambda)
print("Lift per glider = ", corresponding_lift/gliders_per_tether)
print("F_z", F_z)
print("Tether diamter", d_tether)
print("Tether length: ", l_tether)
print("Generators needed: ", 8710000/corresponding_generator_power)
print("Total gliders needed: ", 8710000*gliders_per_tether/corresponding_generator_power)
print("Minimum height: ", corresponding_altitude + corresponding_flight_radius*np.sin(roll_angle))
#print("Generator angle: ", generator_angle)


# print("a = ", a_lst[110])
# print("b = ", b_lst[110])
# print("wind component = ", wind_component_lst[110])
# print("kite radial component = ", kite_radial_component_lst[110])
# print("kite tangential compoent = ", kite_tangential_component_lst[110])



                    
                
    


end = time.time()
print("Elapsed time: ", end-start)