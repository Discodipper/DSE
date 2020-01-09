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
import scipy.linalg as la
start = time.time()

reel_speed_array = np.arange(3, 16, 1) #m/s
altitude_array = np.arange(2000, 3025, 25)

"""theta is polar angle, phi is azimuth angle, beta is elevation angle (operation angle)."""
beta = (np.arange(20, 90, 2))*pi/180 #rad
#phi = (np.arange(0, 95, 5))*pi/180 #rad
azimuth_angle = 0 #rad
hi = 270*pi/180 #rad
g_gravity = 9.80665

drag_coefficient = 0.2
lift_coefficient = 2.0
glide_ratio = lift_coefficient/drag_coefficient
wing_area = 40 #m^2


def reeling_factor(V_reel, V_wind):
    f = V_reel/V_wind
    return(f)

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
    A = la.inv(np.matrix([[np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)], [np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)], [-np.sin(phi), np.cos(phi), 0]]))
    B = A*x
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
    
def tether_force_max(V_a, air_density, lift_coefficient, drag_coefficient, wing_area, elevation_angle):
    F_t_max = (0.5*air_density*resultant_force_coefficient(lift_coefficient, drag_coefficient)* V_a**2 * wing_area)
    F_t_max_horizontal = F_t_max * np.sin(elevation_angle)
    return(F_t_max, F_t_max_horizontal)

def tether_mass_guess(tether_diameter_guess, tether_density, altitude, elevation_angle, F_t_max_horizontal)
    tethermass_per_length = pi*(tether_diameter_guess/2)**2 * tether_density
    tether_length = cable_sag_calculator(tethermass_per_length, altitude, elevation_angle, F_t_max_horizontal)
    tethermass_guess = tethermass_per_length * tether_length
    return(tethermass_guess)

def force_z_direction(V_a_z, air_density, wing_area, drag_coefficient):
    #V_a_z = apparent_wind_speed_cartesian.item(2)
    F_z = 0.5 * air_density * V_a_z**2 * wing_area * drag_coefficient
    return(F_z)

def net_tether_force(V_a, air_density, lift_coefficient, drag_coefficient, wing_area, elevation_angle, V_a_z, glidermass, g_gravity, elevation_angle, F_z):
    F_t_net = tether_force_max(V_a, air_density, lift_coefficient, drag_coefficient, wing_area, elevation_angle, V_a_z) - (glidermass*g_gravity*np.cos(elevation_angle)) - (tethermass*g_gravity) - (F_z*np.cos(elevation_angle))
    F_t_net_z = F_t_net * np.sin(elevation_angle)
    return(F_t_net, F_t_net_z)

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


max_apparent_wind_speed_magnitude = 0
corresponding_apparent_wind_speed = 0
corresponding_altitude = 0
corresponding_beta = 0
corresponding_reel_speed = 0
corresponding_wind_speed = 0
corresponding_tether_force = 0
corresponding_generator_power = 0



operation_angle_fail_lst = []
azimuth_angle_fail_lst = []
for reelspeed in reel_speed_array:
    for altitude in altitude_array:
        for operation_angle in beta:
            temperature, pressure, air_density, windspeed = isa(altitude)
            reel_factor = reeling_factor(reelspeed, windspeed)
            polar_angle = polarangle(operation_angle) #rad
            
            
            if max_elevation_angle(glide_ratio, reel_factor) <= operation_angle and azimuth_constraint(polar_angle, azimuth_angle, glide_ratio, reel_factor):
                
                Lambda = tangential_velocity_factor(polar_angle, azimuth_angle, hi, glide_ratio, reel_factor)[0]
                
                apparent_wind_speed_spherical, apparent_wind_speed_cartesian = apparent_wind_speed_values(windspeed, reelspeed, Lambda, polar_angle, azimuth_angle, hi)

                if isnan(Lambda)==False:
                    magnitude_apparent_wind_speed = apparent_wind_speed_magnitude(apparent_wind_speed_cartesian)
                    V_a_z = apparent_wind_speed_cartesian.item(2)
                    
                    
                    tether_diameter_initial_guess = 0.07 #m
                    total_tether_force, total_tether_force_horizontal = tether_force_max(V_a, air_density, lift_coefficient, drag_coefficient, wing_area, operation_angle, V_a_z)
                    force_in_z = force_z_direction(magnitude_apparent_wind_speed, air_density, wing_area, drag_coefficient)
                    tether_mass_guess_value = tether_mass_guess(tether_diameter_initial_guess, tether_density, altitude, operation_angle, force_in_z)
                    tether_force_net = net_tether_force(magnitude_apparent_wind_speed, air_density, lift_coefficient, drag_coefficient, wing_area, operation_angle)
                    
                    
                    
                    
                    apparent_wind_speed_lst.append(apparent_wind_speed_cartesian)
                    apparent_wind_speed_magnitude_lst.append(apparent_wind_speed_magnitude)
                    reelspeed_lst.append(reelspeed)
                    altitude_lst.append(altitude)
                    operation_angle_lst.append(operation_angle)
                    azimuth_angle_lst.append(azimuth_angle)
                    Lambda_lst.append(Lambda)
                    V_w_lst.append(windspeed)
                    #tether_force_lst.append(tether_force)
                    #generator_power_lst.append(generator_power(tether_force, reelspeed))
                    reel_factor_lst.append(reel_factor)
                    
                    if magnitude_apparent_wind_speed > max_apparent_wind_speed_magnitude:
                        max_apparent_wind_speed_magnitude = magnitude_apparent_wind_speed
                        corresponding_altitude = altitude
                        corresponding_reel_speed = reelspeed
                        corresponding_wind_speed = windspeed
                        #corresponding_tether_force = tether_force
                        #corresponding_generator_power = generator_power(tether_force, reelspeed)
                        corresponding_beta = operation_angle
                        corresponding_apparent_wind_speed = apparent_wind_speed_cartesian
                        corresponding_apparent_wind_speed_spherical = apparent_wind_speed_spherical
                        
            
        else:
            operation_angle_fail_lst.append(operation_angle)
            azimuth_angle_fail_lst.append(azimuth_angle)
            corresponding_theta = polar_angle*180/pi     
            





#print("V_a = ", apparent_wind_speed_lst[7149])
print("Maximum magnitude of V_a = ", max_apparent_wind_speed_magnitude)
print("V_a vector = ", corresponding_apparent_wind_speed)
print("alitutude = ", corresponding_altitude)
print("V_reel = ", corresponding_reel_speed)
print("wind speed = ", corresponding_wind_speed)
print("tether force = ", corresponding_tether_force)
print("generator power = ", corresponding_generator_power)
print("Elevation angle = ", corresponding_beta*180/pi)
# print("a = ", a_lst[110])
# print("b = ", b_lst[110])
# print("wind component = ", wind_component_lst[110])
# print("kite radial component = ", kite_radial_component_lst[110])
# print("kite tangential compoent = ", kite_tangential_component_lst[110])



                    
                
    


end = time.time()
print("Elapsed time: ", end-start)