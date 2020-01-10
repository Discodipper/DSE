# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 16:13:23 2020

@author: Thomas Arblaster
"""

import numpy as np
import matplotlib.pyplot as plt
import statistics as stat
from scipy.optimize import curve_fit

diameter_of_ropes = [2.5,3,4,5,6,7,8,9,10,11,12,13,15,17] #[mm]

strength_of_ropes = [9.9,14.8,24.3,31.4,44.8,73.7,92.1,106,124,148,177,207,259,323] #[kN]

weight_of_ropes = [4.5,6.8,11.1,15.6,22.3,35.6,44.5,54,63,75.5,90,107,134,184] #[g/m]

strength_of_fibre = (3.3 + 3.9)/2 #[GPa]

density_of_fibre = 975 #[kg/m^3]

area_of_ropes = [] #[m^2]

for d in diameter_of_ropes:
   area_of_ropes.append((d*0.001)**2*3.1415)
   
density_of_ropes = [] #[kg/m^2]

i = 0
for w in weight_of_ropes:
    density_of_ropes.append(w*0.001/area_of_ropes[i])
    i = i+1

ult_tension_of_ropes = [] #[GPa]

i = 0
for s in strength_of_ropes:
    ult_tension_of_ropes.append((s*1000/area_of_ropes[i])/(10**9))
    i = i+1
    
def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B

x_data_rope = diameter_of_ropes
y_data_rope_1 = ult_tension_of_ropes
y_data_rope_2 = density_of_ropes

x_values_range = range(round(min(x_data_rope)),max(x_data_rope))

popt1, pcov1 = curve_fit(f, x_data_rope, y_data_rope_1) 
popt2, pcov2 = curve_fit(f, x_data_rope, y_data_rope_2) 

fig, plt_tensile_strengths_1 = plt.subplots()

#plt.subplot(2, 2, 1)8
plt_tensile_strengths_1.set_xlabel('diameter (mm)')
plt_tensile_strengths_1.set_ylabel('tensile strength (GPa)')
plt_tensile_strengths_1.plot(diameter_of_ropes,ult_tension_of_ropes,'ko-')
#plt_tensile_strengths_1.plot(x_values_range,f(np.asarray(x_values_range),popt1[0],popt1[1]))

y_data_rope_3 = np.asarray(ult_tension_of_ropes)/strength_of_fibre
y_data_rope_4 = np.asarray(density_of_ropes)/density_of_fibre

popt3, pcov3 = curve_fit(f, x_data_rope, y_data_rope_3) 
popt4, pcov4 = curve_fit(f, x_data_rope, y_data_rope_4) 

plt_tensile_strengths_2 = plt_tensile_strengths_1.twinx()

#plt.subplot(2, 2, 3)
plt_tensile_strengths_2.set_ylabel('tensile strength (%fibre)')
#plt_tensile_strengths_2.plot(diameter_of_ropes,y_data_rope_3*100,'ko-')
plt_tensile_strengths_2.plot(x_values_range,(f(np.asarray(x_values_range),popt3[0],popt3[1]))*100)

#specific_strength = [b/m for b,m in zip(ult_tension_of_ropes,density_of_ropes)]
#popt3, pcov3 = curve_fit(f, x_data_rope, specific_strength) 
#
#plt.subplot(1, 3, 3)
#plt.xlabel('rope diameter (mm)')
#plt.ylabel('rope tensile strength over density (GPa*m^3/kg)')
#plt.plot(diameter_of_ropes,specific_strength,'ko-')
#plt.plot(x_values_range,f(np.asarray(x_values_range),popt3[0],popt3[1]))

plt.tight_layout()
plt.show()

fig, plt_density_1 = plt.subplots()
plt_density_2 = plt_density_1.twinx()

#plt.subplot(2, 2, 2)
plt_density_1.set_xlabel('diameter (mm)')
plt_density_1.set_ylabel('density (kg/m^3)')
plt_density_1.plot(diameter_of_ropes,density_of_ropes,'ko-')
#plt_density_1.plot(x_values_range,f(np.asarray(x_values_range),popt2[0],popt2[1]))

#plt.subplot(2, 2, 4)
#plt.xlabel('diameter (mm)')
plt_density_2.set_ylabel('density (%fibre)')
#plt_density_2.plot(diameter_of_ropes,y_data_rope_4*100,'ko-')
plt_density_2.plot(x_values_range,(f(np.asarray(x_values_range),popt4[0],popt4[1]))*100)

print('average density = ',stat.mean(density_of_ropes),' kg/m^3')
print('average strength = ',stat.mean(ult_tension_of_ropes),' GPa')

print('equation for density change = ',popt4[0], '*x + ',popt4[1])
print('equation for ultimate strength change = ',popt3[0], '*x + ',popt3[1])