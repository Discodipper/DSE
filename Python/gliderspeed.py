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
start = time.time()

generator_power = 800*1000 #W
reel_speed_array = np.arange(3, 26, 1) #m/s
altitude_array = np.arange(2000, 2510, 10)

"""theta is polar angle, phi is azimuth angle, beta is elevation angle (operation angle)."""
operation_angle = np.arange(20, 46, 1) #deg
theta = 90- operation_angle
phi = np.arange(0, 92, 2) #deg


drag_coefficient = 0.05
lift_coefficient = 2.0
wing_area = 40 #m^2
resultant_force_coefficient = np.sqrt(drag_coefficient**2 + lift_coefficient**2)


for reelspeed in reel_speed_array:
    for altitude in altitude_array:
        for 
        temperature, pressure, air_density, windspeed = isa(altitude)
        
        
        tether_force = generator_power/reelspeed #N
        apparent_speed = np.sqrt(tether_force/(0.5*air_density*resultant_force_coefficient*wing_area))
    


end = time.time()
print("Elapsed time: ", end-start)