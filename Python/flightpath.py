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

# =============================================================================
# power_lst = []
# =============================================================================
# =============================================================================
# flight_radius_lst = []
# mass_area_lst = []
# =============================================================================
# =============================================================================
# normalised_power_lst = []
# =============================================================================


for flight_radius in flight_radius_array:
    normalised_power_lst = []
    flight_radius_lst = []    
    mass_area_lst = []
    power_lst = []
    for mass_area in mass_area_array:
        #power = 1 / flight_radius / lift_coefficient * resultant_force_coefficient * (3*windspeed)**2 * wing_area * reel_speed
        normalised_power = (1-(2*mass_area / (flight_radius * air_density * lift_coefficient))**2)**(3/2)
        normalised_power_lst.append(normalised_power)
        power = 4/27 * normalised_power * 0.5 * air_density * windspeed**3 * wing_area * lift_coefficient**3 / drag_coefficient**2
        power_lst.append(power)
        flight_radius_lst.append(flight_radius)
        mass_area_lst.append(mass_area)
    
    
    plt.plot(mass_area_lst,power_lst)
    plt.plot(mass_area_lst,normalised_power_lst)
plt.show()





