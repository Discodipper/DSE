# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 11:52:22 2020

@author: Mart

This script describes the main part of defining material and properties of the wing box
"""

import math as m
import matplotlib.pyplot as plt
import scipy as sp
import copy
import sys, os
sys.path.clear()
sys.path.append(os.path.realpath('..\\wing'))
sys.path.append(os.path.realpath('..\\Flight_Performance\\Shape\\weight_and_cg_estimation'))

from shape_parameters import sweep_angle_deg, wing_area, wing_span, chord_root, chord_tip
from wing_box_stress_calculations import bending_stress_wing_box, thickness_required_bending
from wing_loads_and_deflection_calculations import wing_box_segmentation, wing_segmentation, determine_Kq

wing_box_chord_root = 0.5 #m
wing_box_chord_tip = 0.4 #m
wing_box_height_root = 0.2 #m
wing_box_height_tip = 0.15 #m
number_of_wing_segments = 1000 #-
span = wing_span #m
lift = 250000 #N
load_factor = 1 #-
wing_weight = 0 #N
wing_surface_area = wing_area #m^2
chord_root = chord_root #m
chord_tip = chord_tip #m
E = 70*10**9 #Pa
drag_wing = lift/30 #N
moment_of_inertia_y = 0.0003333333 #dit moet nog flink aangepast en geparametriseerd worden
moment_of_inertia_x = 0.001
safety_factor = 3
yield_stress_plate = 1000000
density_plates = 975
#comment: the wing surface area has to be coherent with the span and such, otherwise it messes up
#print(wing_area)
aspect_ratio = 12
x = 0.2
chord_root = m.sqrt(wing_area/aspect_ratio)*(0.5/(0.35+0.15*x))
#print(chord_root)
chord_tip = 0.4*chord_root
wing_area_total = 60
wing_area_outside_fuselage = wing_area-x*m.sqrt(wing_area)*m.sqrt(aspect_ratio)*chord_root
#print(wing_area_outside_fuselage)
wing_area = wing_area_outside_fuselage
span = (1-x)*m.sqrt(wing_area*aspect_ratio)
#print(span)
fuselage_weight = 100

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

chord_lengths, spanwise_locations = wing_segmentation(chord_root, chord_tip, number_of_wing_segments, span)
#plt.plot(spanwise_locations, chord_lengths)

Kq = determine_Kq(lift, load_factor, wing_weight, wing_surface_area)
local_load_distribution = [0]
local_drag_distribution = [0]
local_shear_distribution = [0]
shear_force_previous = 0
local_drag_shear_distribution = [0]
drag_shear_force_previous = 0
local_load_previous = 0
bending_moment_distribution = [0]
bending_moment_previous = 0
drag_bending_moment_distribution = [0]
drag_bending_moment_previous = 0 
plate_thickness_list = [0] 
local_drag_previous = 0
local_shear_previous = 0
local_drag_shear_previous = 0
wing_box_chord_height_list = [0]
bending_stress_list = [0]
wing_box_chord_width_list = [0]

for n in range(1, number_of_wing_segments):
    chord_length = chord_lengths[-1-n]
    wing_box_chord_height = chord_length*0.19
    wing_box_chord_width = chord_length*0.45
    local_load = Kq*chord_length
    local_load_distribution.append(local_load)
    local_drag = Kq*chord_length*drag_wing/lift
    local_drag_distribution.append(local_drag)
    local_shear = shear_force_previous -(local_load+local_load_previous)/2*(spanwise_locations[-n-1]-spanwise_locations[-n])
    local_drag_shear = drag_shear_force_previous -(local_drag+local_drag_previous)/2*(spanwise_locations[-n-1]-spanwise_locations[-n])
    local_bending_moment = bending_moment_previous - (local_shear+local_shear_previous)/2*(spanwise_locations[n]-spanwise_locations[n-1])
    local_drag_bending_moment = drag_bending_moment_previous - (local_drag_shear+local_drag_shear_previous)/2*(spanwise_locations[-n-1]-spanwise_locations[-n])
    thickness_required_bending = 0.001
    running = True
    while running == True:
        #bending_stress = -local_bending_moment*0.5*wing_box_chord_height/2/(thickness_required_bending**(3)/12*wing_box_chord_width+wing_box_chord_width*thickness_required_bending*(0.5*wing_box_chord_height-thickness_required_bending/2)**(2))-local_drag_bending_moment*0.5*wing_box_chord_width/2/(wing_box_chord_width**(3)/12*thickness_required_bending)
        bending_stress = -0.5*local_bending_moment*(0.5*wing_box_chord_height-thickness_required_bending/2)/(1/12*thickness_required_bending**3*wing_box_chord_width+thickness_required_bending*wing_box_chord_width*(wing_box_chord_height/2-thickness_required_bending/2))+0.5*local_drag_bending_moment*0.5*wing_box_chord_width/(1/12*wing_box_chord_width**3*thickness_required_bending)
        thickness_required_bending = thickness_required_bending + 0.001
        if bending_stress < yield_stress_plate*safety_factor: 
            running = False
    
    drag_bending_moment_distribution.append(local_drag_bending_moment)
    drag_bending_moment_previous = local_drag_bending_moment
    shear_force_previous = local_shear
    local_load_previous = local_load
    local_drag_previous = local_drag
    local_shear_previous = local_shear
    local_drag_shear_previous = local_drag_shear
    local_shear_distribution.append(local_shear)
    drag_shear_force_previous = local_drag_shear
    local_drag_shear_distribution.append(local_drag_shear)
    bending_moment_previous = local_bending_moment
    bending_moment_distribution.append(local_bending_moment)
    plate_thickness_list.append(thickness_required_bending)
    wing_box_chord_height_list.append(wing_box_chord_height)
    wing_box_chord_width_list.append(wing_box_chord_width)
    bending_stress_list.append(bending_stress)

local_load_distribution.reverse()
local_drag_distribution.reverse()
local_shear_distribution.reverse()
local_drag_shear_distribution.reverse()
bending_moment_distribution.reverse()
drag_bending_moment_distribution.reverse()
plate_thickness_list.reverse()
wing_box_chord_height_list.reverse()
bending_stress_list.reverse()

mass_wing_box = plate_thickness_list[0]*(wing_box_chord_width_list[0]+wing_box_chord_width_list[-2])/2*spanwise_locations[-1]*density_plates*2
print('plate thickness') 
print(plate_thickness_list[0])
print('mass wing box outside fuselage')
print( mass_wing_box*2)
#print(wing_box_chord_height_list[0])
#plt.plot(spanwise_locations, local_load_distribution) 
#plt.plot(spanwise_locations, local_shear_distribution)   
plt.plot(spanwise_locations, bending_moment_distribution)
#plt.plot(spanwise_locations, plate_thickness_list)  
#plt.plot(spanwise_locations, wing_box_chord_height_list)  
#plt.plot(spanwise_locations, bending_stress_list)
#plt.plot(spanwise_locations, drag_bending_moment_distribution)

#plate_thickness = thickness_required_bending(bending_moment_y, bending_moment_x, yield_stress_plate, safety_factor, wing_box_chord_height, wing_box_chord_width)

lift_half_plane_middle_section = lift*(wing_area_total-wing_area_outside_fuselage)/wing_area_total*0.5
bending_moment_middle_section = -x/2/2*m.sqrt(wing_area*aspect_ratio)*lift_half_plane_middle_section+x/2*m.sqrt(wing_area*aspect_ratio)*fuselage_weight*9.81+bending_moment_distribution[0]+local_shear_distribution[0]*-x/2*m.sqrt(wing_area*aspect_ratio)
drag_bending_moment_middle_section = bending_moment_middle_section*drag_wing/lift
print(bending_moment_middle_section)
thickness_required_bending_middle_section = 0.001
running = True
while running == True:
       #bending_stress = -local_bending_moment*0.5*wing_box_chord_height/2/(thickness_required_bending**(3)/12*wing_box_chord_width+wing_box_chord_width*thickness_required_bending*(0.5*wing_box_chord_height-thickness_required_bending/2)**(2))-local_drag_bending_moment*0.5*wing_box_chord_width/2/(wing_box_chord_width**(3)/12*thickness_required_bending)
       #bending_stress = -0.5*bending_moment_middle_section*(0.5*wing_box_chord_height_list[0]-thickness_required_bending_middle_section/2)/(1/12*thickness_required_bending_middle_section**3*wing_box_chord_width_list[0]+thickness_required_bending_middle_section*wing_box_chord_width_list[0]*(wing_box_chord_height_list[0]/2-thickness_required_bending_middle_section/2))+0.5*drag_bending_moment_middle_section*0.5*wing_box_chord_width_list[0]/(1/12*wing_box_chord_width_list[0]**3*thickness_required_bending_middle_section)
       bending_stress = -0.5*bending_moment_middle_section*(0.5*wing_box_chord_height-thickness_required_bending_middle_section/2)/(1/12*thickness_required_bending_middle_section**3*wing_box_chord_width+thickness_required_bending_middle_section*wing_box_chord_width*(wing_box_chord_height/2-thickness_required_bending_middle_section/2))+0.5*drag_bending_moment_middle_section*0.5*wing_box_chord_width/(1/12*wing_box_chord_width**3*thickness_required_bending_middle_section)
       thickness_required_bending_middle_section = thickness_required_bending_middle_section + 0.001
       if bending_stress < yield_stress_plate*safety_factor: 
           running = False
print('plate thickness middle section')
print(thickness_required_bending_middle_section)
mass_middle_section = thickness_required_bending_middle_section*x*m.sqrt(wing_area*aspect_ratio)*wing_box_chord_width*2*density_plates
print('mass middle section')
print(mass_middle_section)
print(wing_box_chord_width)