# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 21:36:16 2020

@author: Thomas Arblaster
"""

import numpy as np
import matplotlib.pyplot as plt

L_den_noise_allowable = 47 #dBA

noise_from_take_off_1 = 50 #dBA
noise_from_take_off_2 = 57 #dBA
noise_from_take_off_3 = 64 #dBA
noise_from_take_off_4 = 70 #dBA
noise_from_landing = 56 #dBA

noise_ambient = 30 #dBA

variable_minutes = 3

T_day = 12*60*variable_minutes #min
T_evening = 3*60*variable_minutes #min
T_night = 9*60*variable_minutes #min

#10**(L_den_noise_allowable/10)*24

possible_l_den_noises_1 = []
possible_l_den_noises_2 = []
possible_l_den_noises_3 = []
possible_l_den_noises_4 = []
x_values = []

range_m = 80

for m in range(range_m):
#    noises_day = noise_ambient + np.zeros(T_day)
#    noises_day[:m] = noise_from_take_off
#    noises_evening = noise_ambient + np.zeros(T_evening)
#    noises_evening[:m] = noise_from_take_off
#    noises_night = noise_ambient + np.zeros(T_night)
#    noises_night[:m] = noise_from_take_off
    
    noises_day = noise_ambient + np.zeros(T_day)
    noises_day[:(4*m)] = noise_from_take_off_1
    noises_evening = noise_ambient + np.zeros(T_evening)
    noises_evening[:m] = noise_from_take_off_1
    noises_night = noise_ambient + np.zeros(T_night)
    noises_night[:(3*m)] = noise_from_take_off_1
    
    l_a_day_sum = []
    l_a_evening_sum = []
    l_a_night_sum = []
    
    for n in noises_day:
        l_a_day_sum.append(10**(0.1*n))
        
    for n in noises_evening:
        l_a_evening_sum.append(10**(0.1*n))
        
    for n in noises_night:
        l_a_night_sum.append(10**(0.1*n))
    
    l_a_eq_t_day = 10*np.log10((1/T_day)*(sum(l_a_day_sum)))
    l_a_eq_t_evening = 10*np.log10((1/T_evening)*(sum(l_a_evening_sum)))
    l_a_eq_t_night = 10*np.log10((1/T_night)*(sum(l_a_night_sum)))
    
    l_den = 10*np.log10(1/24 * (12*(10**(0.1*l_a_eq_t_day)) + 3*(10**(0.1*(l_a_eq_t_evening + 5))) + 9*(10**(0.1*(l_a_eq_t_night + 10)))))
    possible_l_den_noises_1.append(l_den)
    x_values.append(m*8/variable_minutes)
    
for m in range(range_m):
    noises_day = noise_ambient + np.zeros(T_day)
    noises_day[:(4*m)] = noise_from_take_off_2
    noises_evening = noise_ambient + np.zeros(T_evening)
    noises_evening[:m] = noise_from_take_off_2
    noises_night = noise_ambient + np.zeros(T_night)
    noises_night[:(3*m)] = noise_from_take_off_2
    
    l_a_day_sum = []
    l_a_evening_sum = []
    l_a_night_sum = []
    
    for n in noises_day:
        l_a_day_sum.append(10**(0.1*n))
        
    for n in noises_evening:
        l_a_evening_sum.append(10**(0.1*n))
        
    for n in noises_night:
        l_a_night_sum.append(10**(0.1*n))
    
    l_a_eq_t_day = 10*np.log10((1/T_day)*(sum(l_a_day_sum)))
    l_a_eq_t_evening = 10*np.log10((1/T_evening)*(sum(l_a_evening_sum)))
    l_a_eq_t_night = 10*np.log10((1/T_night)*(sum(l_a_night_sum)))
    
    l_den = 10*np.log10(1/24 * (12*(10**(0.1*l_a_eq_t_day)) + 3*(10**(0.1*(l_a_eq_t_evening + 5))) + 9*(10**(0.1*(l_a_eq_t_night + 10)))))
    possible_l_den_noises_2.append(l_den)
    
for m in range(range_m):
    noises_day = noise_ambient + np.zeros(T_day)
    noises_day[:(4*m)] = noise_from_take_off_3
    noises_evening = noise_ambient + np.zeros(T_evening)
    noises_evening[:m] = noise_from_take_off_3
    noises_night = noise_ambient + np.zeros(T_night)
    noises_night[:(3*m)] = noise_from_take_off_3
    
    l_a_day_sum = []
    l_a_evening_sum = []
    l_a_night_sum = []
    
    for n in noises_day:
        l_a_day_sum.append(10**(0.1*n))
        
    for n in noises_evening:
        l_a_evening_sum.append(10**(0.1*n))
        
    for n in noises_night:
        l_a_night_sum.append(10**(0.1*n))
    
    l_a_eq_t_day = 10*np.log10((1/T_day)*(sum(l_a_day_sum)))
    l_a_eq_t_evening = 10*np.log10((1/T_evening)*(sum(l_a_evening_sum)))
    l_a_eq_t_night = 10*np.log10((1/T_night)*(sum(l_a_night_sum)))
    
    l_den = 10*np.log10(1/24 * (12*(10**(0.1*l_a_eq_t_day)) + 3*(10**(0.1*(l_a_eq_t_evening + 5))) + 9*(10**(0.1*(l_a_eq_t_night + 10)))))
    possible_l_den_noises_3.append(l_den)
    
for m in range(range_m):
    noises_day = noise_ambient + np.zeros(T_day)
    noises_day[:(4*m)] = noise_from_take_off_4
    noises_evening = noise_ambient + np.zeros(T_evening)
    noises_evening[:m] = noise_from_take_off_4
    noises_night = noise_ambient + np.zeros(T_night)
    noises_night[:(3*m)] = noise_from_take_off_4
    
    l_a_day_sum = []
    l_a_evening_sum = []
    l_a_night_sum = []
    
    for n in noises_day:
        l_a_day_sum.append(10**(0.1*n))
        
    for n in noises_evening:
        l_a_evening_sum.append(10**(0.1*n))
        
    for n in noises_night:
        l_a_night_sum.append(10**(0.1*n))
    
    l_a_eq_t_day = 10*np.log10((1/T_day)*(sum(l_a_day_sum)))
    l_a_eq_t_evening = 10*np.log10((1/T_evening)*(sum(l_a_evening_sum)))
    l_a_eq_t_night = 10*np.log10((1/T_night)*(sum(l_a_night_sum)))
    
    l_den = 10*np.log10(1/24 * (12*(10**(0.1*l_a_eq_t_day)) + 3*(10**(0.1*(l_a_eq_t_evening + 5))) + 9*(10**(0.1*(l_a_eq_t_night + 10)))))
    possible_l_den_noises_4.append(l_den)
    
#print(possible_l_den_noises)

plt.xlabel('daily exposure to sound of active propellers (min)')
plt.ylabel('L_den (dBA)')
plt.plot(x_values,possible_l_den_noises_1,'k.',label = str(noise_from_take_off_1)+' dBA')
plt.plot(x_values,possible_l_den_noises_2,'k:',label = str(noise_from_take_off_2)+' dBA') 
plt.plot(x_values,possible_l_den_noises_3,'k--',label = str(noise_from_take_off_3)+' dBA')
plt.plot(x_values,possible_l_den_noises_4,'k-.',label = str(noise_from_take_off_4)+' dBA')
plt.plot(x_values,np.zeros(len(x_values))+L_den_noise_allowable,'k-',label = 'maximum allowable L_den')
plt.legend()


