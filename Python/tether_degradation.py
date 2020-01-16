# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 22:37:42 2020

@author: Thomas Arblaster
"""

import numpy as np
import matplotlib.pyplot as plt

joint_efficiency = 0.8

abrasion_factor = 0.9
abrasion_time = 500 #hr

fatigue_factor = 0.75
fatigue_time = 2200 #hr

moisture_factor = 1

uv_factor = 0.8
uv_time = 1050 #hr

dyneema_safety = 1.5
dyneema_dynamic = 1.1
dyneema_asymetric = 1.2

step_size = 100

hours_of_degradation = np.arange(0, 2100, step_size).tolist()

joint_array = np.zeros(len(hours_of_degradation)) + joint_efficiency
abrasion_array = 1 + np.zeros(len(hours_of_degradation)) - np.asarray(hours_of_degradation)*(1 - abrasion_factor)/abrasion_time
fatigue_array = 1 + np.zeros(len(hours_of_degradation)) - np.asarray(hours_of_degradation)*(1 - fatigue_factor)/fatigue_time
moisture_array = np.zeros(len(hours_of_degradation)) + moisture_factor
uv_array = 1 + np.zeros(len(hours_of_degradation)) - np.asarray(hours_of_degradation)*(1 - uv_factor)/uv_time

safety_factor = 1.5
dynamic_factor = 1.1
asymetric_factor = 1.2

design_factor_array = np.zeros(len(hours_of_degradation)) + safety_factor*dynamic_factor*asymetric_factor
i = 0
for (j, a, f, m, u) in zip(joint_array, abrasion_array, fatigue_array, moisture_array, uv_array):
    design_factor_array[i] = design_factor_array[i]/(j*a*f*m*u)
    i = i + 1

fig, ax1 = plt.subplots()

ax1.set_xlabel('time (hrs)')
ax1.set_ylabel('degradation factor (-)')
ax1.plot(hours_of_degradation,joint_array,'k--',label = 'joint efficiency')
ax1.plot(hours_of_degradation,abrasion_array,'k-.',label = 'abrasion')
ax1.plot(hours_of_degradation,fatigue_array,'k:',label = 'fatigue')
ax1.plot(hours_of_degradation,moisture_array,'k-.<',label = 'moisture')
ax1.plot(hours_of_degradation,uv_array,'k:>',label = 'UV')
ax1.plot(0,1,'k-',label = 'design factor')
axes = plt.gca()
axes.set_ylim([0.5,1.1])

plt.legend()

ax2 = ax1.twinx()

ax2.set_xlabel('design factor (-)')
ax2.plot(hours_of_degradation,design_factor_array,'k-',label = 'design factor')

fig.tight_layout()
plt.show()
