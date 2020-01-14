# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 13:57:02 2020

@author: Thomas Arblaster
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#diameter_of_ropes = [2.5,3,4,5,6,7,8,9,10,11,12,13,15,17] #[mm]

#strength_of_ropes = [9.9,14.8,24.3,31.4,44.8,73.7,92.1,106,124,148,177,207,259,323] #[kN]

#weight_of_ropes = [4.5,6.8,11.1,15.6,22.3,35.6,44.5,54,63,75.5,90,107,134,184] #[g/m]
#weight_of_ropes = np.asarray(weight_of_ropes)*0.001 #[kg/m]

diameter_of_ropes = [5,6,7,8,9,10,11,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,48,52,56,60,64,72,80,88,96] #[mm]

strength_of_ropes = [3.3,4,6.1,8.1,9,10.9,13.9,17.8,22,26.1,36,41,50.5,55,65,70,78,84.5,94,110,133,145,155,170,185,220,250,280,325,400,450,525,625] #[1000kg]
strength_of_ropes = np.asarray(strength_of_ropes)*9.81 #[kN]

weight_of_ropes = [1.4,2.1,2.7,3.5,4.2,4.6,6.2,9.0,12.0,14.0,19.0,22,25.5,30,32.5,39,43,46,53,60,68,73,83,100,125,143,175,200,230,300,350,430,500] #[kg/100m]
weight_of_ropes = np.asarray(weight_of_ropes)*0.01 #[kg/m]

x_data_rope = weight_of_ropes
y_data_rope_1 = strength_of_ropes

def f(x, A): # this is your 'straight line' y=f(x)
    return A*x

popt1, pcov1 = curve_fit(f, x_data_rope, y_data_rope_1) 
#popt2, pcov2 = curve_fit(f, x_data_rope, y_data_rope_2) 

x_value_range = np.empty(100) #[kN]
b_range_of_strengths = np.arange(0, max(x_data_rope), max(x_data_rope)/100)
ind_i_still_dk_what_this_does = np.arange(len(x_value_range))
np.put(x_value_range, ind_i_still_dk_what_this_does, b_range_of_strengths)

strength_of_fibre = 3.6 #[GPa]
density_of_fibre = 975 #[kg/m^3]

a_strength_array = np.empty(100) #[kN]
b_range_of_strengths = np.arange(1, round(max(strength_of_ropes)), max(strength_of_ropes)/100)
ind_idk_what_this_does = np.arange(len(a_strength_array))
np.put(a_strength_array, ind_idk_what_this_does, b_range_of_strengths)

area_pure_fibre = a_strength_array*1000/(strength_of_fibre*10**9) #[m^2]
weight_pure_fibre = density_of_fibre * area_pure_fibre #[kg/m]

x_data_rope_2 = weight_pure_fibre
y_data_rope_2 = a_strength_array

#plt.subplot(1, 2, 1)
plt.xlabel('weight (kg/m)')
plt.ylabel('tensile strength (kN)')
plt.plot(x_data_rope_2,y_data_rope_2,'k--',label='Dyneema SK75 fibre')
plt.plot(x_data_rope,y_data_rope_1,'ko',label='Dyneema ropes of different thicknesses')
plt.plot(x_value_range,f(np.asarray(x_value_range),popt1[0]),'k-.',label='Rope approximation')
plt.legend()

#a_first_graph = popt1[0]*1000
#b_second_graph = strength_of_fibre*10**9/density_of_fibre
#c_difference_graph = (b_second_graph - a_first_graph)/(1 + a_first_graph*b_second_graph)
#
#b_for_weight_per_area = (b_second_graph/c_difference_graph)*(strength_of_fibre*10**9)

x_value_range = x_value_range*11

weight_fibre_per_diameter = [] #kg/m^2
weight_rope_per_diameter = [] #kg/m^2

for x in x_value_range: #AREA OF PURE FIBRE
   weight_fibre_per_diameter.append(density_of_fibre*x)
   weight_rope_per_diameter.append(strength_of_fibre*10**9*x/(popt1[0]*1000)) 
   
x_value_range = x_value_range[2:]
weight_fibre_per_diameter = weight_fibre_per_diameter[2:]
weight_rope_per_diameter = weight_rope_per_diameter[2:]

y_values_weight_difference = [b/a for a,b in zip(weight_fibre_per_diameter,weight_rope_per_diameter)]
   
#plt.subplot(1, 2, 2)
#plt.xlabel('diameter of pure fibre (m^2)')
#plt.ylabel('weight increase (%)')
#plt.plot(x_value_range,y_values_weight_difference,'k-',label='Rope performance') #,label='Rope performance'
##plt.plot(x_value_range,weight_fibre_per_diameter,'k-.',label='fibre') #,label='Rope performance'
##plt.plot(x_value_range,weight_rope_per_diameter,'k-',label='Rope') #,label='Rope performance'
#plt.legend()

#plt.subplot(1, 2, 2)
#plt.xlabel('area of pure fibre (10^-2 m^2)')
#plt.ylabel('weight increase (%)')
#plt.plot(x_value_range*10**2,weight_change_per_area_percentages,'k-') #,label='Rope performance'
#plt.legend()

plt.show()

def tether_rope_calculator(required_tether_force_newton,tether_material_tension_strength_pa,tether_material_density):
    pure_fibre_area = required_tether_force_newton/tether_material_tension_strength_pa
    rope_weight_per_m = tether_material_density*pure_fibre_area*(1 + pure_fibre_area/0.0005037987746501848)
    rope_area = rope_weight_per_m/tether_material_density
    return rope_weight_per_m, rope_area
