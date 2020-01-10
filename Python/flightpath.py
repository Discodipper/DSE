# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 11:31:43 2020

@author: Xander
"""

from math import *
import matplotlib.pyplot as plt
import numpy as np
import time
from ISA_calculator import isa
import scipy.linalg as la
start = time.time()

#Assumption: CL/CD >> 1 and therefore crosswind kitespeed = apparent wind speed

altitude = 2000 #m
temperature, pressure, air_density, windspeed = isa(altitude)

drag_coefficient = 0.2
lift_coefficient = 2.0
glide_ratio = lift_coefficient/drag_coefficient
resultant_force_coefficient = np.sqrt(drag_coefficient**2 + lift_coefficient**2)
wing_area = 60 #m^2
mass_glider = 1800 #kg
mass_area_array = np.arange(0,42,2) #kg/m^2
flight_radius_array = np.arange(25,125,25) #m
reel_speed = 8 #m/s
roll_angle = 90 * pi/180 #rad


def normalised_power(mass_area, flight_radius, air_density, lift_coefficient):
    P_normalised = (1-(2*mass_area / (flight_radius * air_density * lift_coefficient))**2)**(3/2)
    return(P_normalised)



def flight_radius(roll_angle, net_lift, total_mass, kite_velocity):
    flight_radius = total_mass * kite_velocity**2 / net_lift / np.sin(roll_angle)
    return(flight_radius)

def reduced_lift_by_roll(net_lift, roll_angle, total_mass, wing_area, air_density, lift_coefficient, radius_flight_path):
    L_reduced = net_lift * np.cos(roll_angle)
    L_reduced_2 = net_lift * (1-((minimum_radius_flight_path(total_mass, wing_area, air_density, lift_coefficient)) / radius_flight_path)**2)**0.5
    return(L_reduced, L_reduced_2)

def minimum_radius_flight_path(total_mass, wing_area, air_density, lift_coefficient):
    r_min = 2 * total_mass / wing_area / air_density / lift_coefficient
    return(r_min)

def power_produced():
    power = 4/27 * normalised_power1 * 0.5 * air_density * windspeed**3 * wing_area * lift_coefficient**3 / drag_coefficient**2

for flight_radius_iteration in flight_radius_array:
    normalised_power_lst = []
    flight_radius_lst = []    
    mass_area_lst = []
    power_lst = []
    for mass_area in mass_area_array:
        #power = 1 / flight_radius / lift_coefficient * resultant_force_coefficient * (3*windspeed)**2 * wing_area * reel_speed
        normalised_power1 = normalised_power(mass_area, flight_radius_iteration, air_density, lift_coefficient)
        normalised_power_lst.append(normalised_power1)
        power = 4/27 * normalised_power1 * 0.5 * air_density * windspeed**3 * wing_area * lift_coefficient**3 / drag_coefficient**2
        power_lst.append(power)
        flight_radius_lst.append(flight_radius_iteration)
        mass_area_lst.append(mass_area)
    
     
# =============================================================================
#     plt.plot(mass_area_lst,power_lst)
#     plt.plot(mass_area_lst,normalised_power_lst)
#plt.show()
# =============================================================================

    # plt.figure(1)
    plt.subplot(211)
    plt.plot(mass_area_lst,power_lst)
    plt.subplot(212)
    plt.plot(mass_area_lst,normalised_power_lst)






