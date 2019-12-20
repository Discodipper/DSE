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
import scipy.linalg as la
start = time.time()

generator_power = 800*1000 #W
reel_speed_array = np.arange(3, 9, 1) #m/s
altitude_array = np.arange(2000, 2525, 25)

"""theta is polar angle, phi is azimuth angle, beta is elevation angle (operation angle)."""
beta = (np.arange(20, 50, 5))*pi/180 #rad
phi = (np.arange(0, 90, 1))*pi/180 #rad
hi = 90*pi/180 #rad


drag_coefficient = 0.2
lift_coefficient = 2.0
glide_ratio = lift_coefficient/drag_coefficient
wing_area = 40 #m^2
resultant_force_coefficient = np.sqrt(drag_coefficient**2 + lift_coefficient**2)

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

operation_angle_fail_lst = []
azimuth_angle_fail_lst = []
for reelspeed in reel_speed_array:
    for altitude in altitude_array:
        for operation_angle in beta:
            for azimuth_angle in phi:
                temperature, pressure, air_density, windspeed = isa(altitude)
                polar_angle = (pi/2) - operation_angle #rad
                
                reel_factor = reeling_factor(reelspeed, windspeed)
                if azimuth_constraint(polar_angle, azimuth_angle, glide_ratio, reel_factor) < 0:
                    Lambda = tangential_velocity_factor(polar_angle, azimuth_angle, hi, glide_ratio, reel_factor)[0]
                    A = tangential_velocity_factor(polar_angle, azimuth_angle, hi, glide_ratio, reel_factor)[1]
                    B = tangential_velocity_factor(polar_angle, azimuth_angle, hi, glide_ratio, reel_factor)[2]
                    
                    kite_speed_tangential = Lambda * windspeed
                    
                    
                    wind_component = np.matrix([[np.sin(polar_angle) * np.cos(azimuth_angle)], [np.cos(polar_angle) * np.cos(azimuth_angle)], [-np.sin(azimuth_angle)]]) * windspeed
                    kite_radial_component = np.matrix([[1], [0], [0]]) * reelspeed
                    kite_tangential_component = np.matrix([[0], [np.cos(hi)], [np.sin(hi)]]) * kite_speed_tangential
                    
                    apparent_wind_speed = wind_component + kite_radial_component - kite_tangential_component
                    
                    #tether_force = generator_power/reelspeed #N
                    if isnan(Lambda)==False:
                        apparent_wind_speed_magnitude = la.norm(np.array(apparent_wind_speed))
                        if apparent_wind_speed_magnitude > 10:
                            apparent_wind_speed_lst.append(apparent_wind_speed)
                            apparent_wind_speed_magnitude_lst.append(apparent_wind_speed_magnitude)
                            reelspeed_lst.append(reelspeed)
                            altitude_lst.append(altitude)
                            operation_angle_lst.append(operation_angle)
                            azimuth_angle_lst.append(azimuth_angle)
                            Lambda_lst.append(Lambda)
                            a_lst.append(A)
                            b_lst.append(B)
                            b_f_lst.append(B - reel_factor)
                            wind_component_lst.append(wind_component)
                            kite_radial_component_lst.append(kite_radial_component)
                            kite_tangential_component_lst.append(kite_tangential_component)
                            V_w_lst.append(windspeed)
                
            else:
                operation_angle_fail_lst.append(operation_angle)
                azimuth_angle_fail_lst.append(azimuth_angle)
                    

                    
                
    


end = time.time()
print("Elapsed time: ", end-start)