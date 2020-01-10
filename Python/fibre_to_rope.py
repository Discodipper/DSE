# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 13:57:02 2020

@author: Thomas Arblaster
"""

import numpy as np
import matplotlib.pyplot as plt
import statistics as stat
from scipy.optimize import curve_fit

diameter_of_ropes = [2.5,3,4,5,6,7,8,9,10,11,12,13,15,17] #[mm]

strength_of_ropes = [9.9,14.8,24.3,31.4,44.8,73.7,92.1,106,124,148,177,207,259,323] #[kN]
8
weight_of_ropes = [4.5,6.8,11.1,15.6,22.3,35.6,44.5,54,63,75.5,90,107,134,184] #[g/m]
weight_of_ropes = np.asarray(weight_of_ropes)*0.001 #[kg/m]

area_of_ropes = [] #[m^2]

for d in diameter_of_ropes:
   area_of_ropes.append((d*0.001)**2*3.1415)

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

strength_of_fibre = (3.3 + 3.9)/2 #[GPa]
density_of_fibre = 975 #[kg/m^3]

a_strength_array = np.empty(100) #[kN]
b_range_of_strengths = np.arange(1, round(max(strength_of_ropes)), max(strength_of_ropes)/100)
ind_idk_what_this_does = np.arange(len(a_strength_array))
np.put(a_strength_array, ind_idk_what_this_does, b_range_of_strengths)

area_pure_fibre = a_strength_array*1000/(strength_of_fibre*10**9) #[m^2]
weight_pure_fibre = density_of_fibre * area_pure_fibre #[kg/m]

x_data_rope_2 = weight_pure_fibre
y_data_rope_2 = a_strength_array

plt.subplot(1, 2, 1)
plt.xlabel('weight (kg/m)')
plt.ylabel('tensile strength (kN)')
plt.plot(x_data_rope_2,y_data_rope_2,'k--',label='Dyneema SK75 fibre')
plt.plot(x_data_rope,y_data_rope_1,'ko',label='Dyneema D12 MAX rope')
plt.plot(x_value_range,f(np.asarray(x_value_range),popt1[0]),'k-.',label='Rope approximation')
plt.legend()

x_value_range = x_value_range*3

b_for_weight_per_area = (strength_of_fibre*10**6/density_of_fibre-popt1[0])*1000/(strength_of_fibre*10**9)

plt.subplot(1, 2, 2)
plt.xlabel('weight increase (%)')
plt.ylabel('area of pure fibre (10^-5 m^2)')
plt.plot(x_value_range*100,f(np.asarray(x_value_range),b_for_weight_per_area*10**5),'k-') #,label='Rope performance'
#plt.legend()

plt.show()

def tether_rope_calculator(required_tether_force_newton,tether_material_tension_strength_pa,tether_material_density):
    pure_fibre_area = required_tether_force_newton/tether_material_tension_strength_pa
    rope_weight_per_m = tether_material_density*pure_fibre_area*(1 + pure_fibre_area/0.0005037987746501848)
    rope_area = rope_weight_per_m/tether_material_density
    return rope_weight_per_m, rope_area
