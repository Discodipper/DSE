# -*- coding: utf-8 -*-
"""
Created on Mon Jan 06 13:33:76 2020

@author: Mart
"""

import matplotlib.pyplot as plt
import scipy as sp
import copy
import sys
sys.path.append(r'C:\Users\Mart\Documents\GitHub\DSE\Python\Flight Performance\Shape\weight and cg estimation')

"""
This section enters the variables needed in the definitions.
""" 

from shape_parameters import sweep_angle_deg, wing_area, wing_span, chord_root, chord_tip

wing_box_chord_root = 0.5 #m
wing_box_chord_tip = 0.4 #m
wing_box_height_root = 0.2 #m
wing_box_height_tip = 0.15 #m
number_of_wing_segments = 1000 #-
span = wing_span #m
lift = 42000 #N
load_factor = 1 #-
wing_weight = 0 #N
wing_surface_area = wing_area #m^2
chord_root = chord_root #m
chord_tip = chord_tip #m
E = 70*10**9 #Pa
drag_wing = 1000 #N
moment_of_inertia_y = 0.0003333333 #dit moet nog flink aangepast en geparametriseerd worden
moment_of_inertia_x = 0.001
#comment: the wing surface area has to be coherent with the span and such, otherwise it messes up

"""
This piece segments the wing box in a number of bits, calculating the chords of wing box 
"""
def wing_box_segmentation(wing_box_chord_root, wing_box_chord_tip, wing_box_height_root, wing_box_height_tip, number_of_wing_segments, span):
    wing_box_chord_lengths = []
    spanwise_locations = []
    wing_box_heights = []
    for i in range(0,number_of_wing_segments):
        local_chord = wing_box_chord_root - i*(wing_box_chord_root-wing_box_chord_tip)/number_of_wing_segments
        wing_box_chord_lengths.append(local_chord)
        local_spanwise_location = i*span/2/number_of_wing_segments
        spanwise_locations.append(local_spanwise_location)
        local_height = wing_box_height_root - i*(wing_box_height_root-wing_box_height_tip)/number_of_wing_segments
        wing_box_heights.append(local_height)
    return wing_box_chord_lengths, spanwise_locations, wing_box_heights

wing_box_chord_lengths, spanwise_locations, wing_box_heights = wing_box_segmentation(wing_box_chord_root, wing_box_chord_tip, wing_box_height_root, wing_box_height_tip, number_of_wing_segments, span)
#plt.plot(wing_box_segmentation_outcomes[1], wing_box_segmentation_outcomes[2])
#plt.show()

"""
This part enters the segmentation of the wing, calculting the different cords. It is tested to see if it adds up as lift
"""
def wing_segmentation(chord_root, chord_tip, number_of_wing_segments, span):
    chord_lengths = []
    spanwise_locations = []
    for i in range(0,number_of_wing_segments):
        local_chord = chord_root - i*(chord_root-chord_tip)/number_of_wing_segments
        chord_lengths.append(local_chord)
        local_spanwise_location = i*span/2/number_of_wing_segments
        spanwise_locations.append(local_spanwise_location)
    return chord_lengths,spanwise_locations

chord_lengths, spanwise_locations = wing_segmentation(chord_root, chord_tip, number_of_wing_segments, span)
#plt.plot(spanwise_locations, chord_lengths)
#plt.show

testarea = 0
width = span/2/len(chord_lengths)
for n in range(0, len(chord_lengths)):
    testarea = testarea + width*chord_lengths[n]
    
#print(testarea)
    
"""
This part calculates the integrating constant for the wing loadings
"""
def determine_Kq(lift, load_factor, wing_weight, wing_surface_area):
    Kq = (lift - load_factor*wing_weight)/wing_surface_area
    return Kq

Kq = determine_Kq(lift, load_factor, wing_weight, wing_surface_area)
#print(Kq)

"""
This part calculates the distribution of loads over the wing, with a specific laod for every subpart of the wing. It is tested on total lift.
"""
def load_distribution_over_wing(Kq, chord_lengths, drag_wing, lift):
    local_load_distribution = []
    local_drag_distribution = []
    for i in chord_lengths:
        local_load = Kq*i
        local_load_distribution.append(local_load)
        local_drag = Kq*i*drag_wing/lift
        local_drag_distribution.append(local_drag)
    return local_load_distribution, local_drag_distribution

local_load_distribution, local_drag_distribution = load_distribution_over_wing(Kq, chord_lengths, drag_wing, lift)
#plt.plot(spanwise_locations, local_load_distribution)
#plt.plot(spanwise_locations, local_drag_distribution)

#test_total_lift = 0
#for n in range(0, len(local_load_distribution)):
#    test_total_lift = test_total_lift + local_load_distribution[n]*width
    
#print(test_total_lift)

#test_total_drag = 0
#for n in range(0, len(local_load_distribution)):
#    test_total_drag = test_total_drag + local_drag_distribution[n]*width
#print(test_total_drag)

"""
The shear forces are calculated in this part
"""
def shear_force_wing(local_load_distribution, spanwise_locations, local_drag_distribution):
    local_shear_distribution = [0]
    shear_force_previous = 0
    local_drag_shear_distribution = [0]
    drag_shear_force_previous = 0
    spanwise_locations.reverse()
    local_load_distribution.reverse()
    local_drag_distribution.reverse()
    for n in range(1,number_of_wing_segments):
        local_shear = shear_force_previous -(local_load_distribution[n]+local_load_distribution[n-1])/2*(spanwise_locations[n]-spanwise_locations[n-1])
        shear_force_previous = local_shear
        local_shear_distribution.append(local_shear)
        local_drag_shear = drag_shear_force_previous -(local_drag_distribution[n]+local_drag_distribution[n-1])/2*(spanwise_locations[n]-spanwise_locations[n-1])
        drag_shear_force_previous = local_drag_shear
        local_drag_shear_distribution.append(local_drag_shear)
    spanwise_locations.reverse()
    local_load_distribution.reverse()
    local_shear_distribution.reverse()
    local_drag_shear_distribution.reverse()
    return local_shear_distribution, local_drag_shear_distribution

local_shear_distribution, local_drag_shear_distribution = shear_force_wing(local_load_distribution, spanwise_locations, local_drag_distribution)
#plt.plot(spanwise_locations, local_shear_distribution)
#plt.plot(spanwise_locations, local_drag_shear_distribution)

"""
This part calculates the bending moments in the wing caused by the shear
"""    
def bending_moment_distribution(local_shear_distribution, spanwise_locations, local_drag_shear_distribution):
    bending_moment_distribution = [0]
    bending_moment_previous = 0
    drag_bending_moment_distribution = [0]
    drag_bending_moment_previous = 0
    local_drag_shear_distribution.reverse()
    spanwise_locations.reverse()
    local_shear_distribution.reverse()
    for n in range(1,number_of_wing_segments):
        local_bending_moment = bending_moment_previous - (local_shear_distribution[n]+local_shear_distribution[n-1])/2*(spanwise_locations[n]-spanwise_locations[n-1])
        bending_moment_distribution.append(local_bending_moment)
        bending_moment_previous = local_bending_moment
        local_drag_bending_moment = drag_bending_moment_previous - (local_drag_shear_distribution[n]+local_drag_shear_distribution[n-1])/2*(spanwise_locations[n]-spanwise_locations[n-1])
        drag_bending_moment_distribution.append(local_drag_bending_moment)
        drag_bending_moment_previous = local_drag_bending_moment
    spanwise_locations.reverse()
    local_shear_distribution.reverse()
    bending_moment_distribution.reverse()
    drag_bending_moment_distribution.reverse()
    return bending_moment_distribution, drag_bending_moment_distribution

bending_moment_distribution, drag_bending_moment_distribution = bending_moment_distribution(local_shear_distribution, spanwise_locations, local_drag_shear_distribution)
#plt.plot(spanwise_locations, bending_moment_distribution)
#plt.plot(spanwise_locations, drag_bending_moment_distribution)

"""
This part calculates the deflection angle at every part
"""
def local_deflection_angle(spanwise_locations, bending_moment_distribution, E, moment_of_inertia_y, drag_bending_moment_distribution, moment_of_inertia_x):
    deflection_angles = [0]
    deflection_angle_previous = 0
    drag_deflection_angles = [0]
    drag_deflection_angle_previous = 0
    for n in range(1,number_of_wing_segments):
        local_deflection_angle = deflection_angle_previous + 0.5*(bending_moment_distribution[n]/E/moment_of_inertia_y+bending_moment_distribution[n-1]/E/moment_of_inertia_y)*(spanwise_locations[n]-spanwise_locations[n-1])
        deflection_angles.append(local_deflection_angle)
        deflection_angle_previous = local_deflection_angle
        local_drag_deflection_angle = drag_deflection_angle_previous + 0.5*(drag_bending_moment_distribution[n]/E/moment_of_inertia_x+drag_bending_moment_distribution[n-1]/E/moment_of_inertia_x)*(spanwise_locations[n]-spanwise_locations[n-1])
        drag_deflection_angles.append(local_drag_deflection_angle)
        drag_deflection_angle_previous = local_drag_deflection_angle
    return deflection_angles, drag_deflection_angles

deflection_angles, drag_deflection_angles = local_deflection_angle(spanwise_locations, bending_moment_distribution, E, moment_of_inertia_y, drag_bending_moment_distribution, moment_of_inertia_x)

"""
This part calculates the deflection of the wing due to the bending/loads
"""
def local_deflection(spanwise_locations, deflection_angles, drag_deflection_angles):
    deflections = [0]
    deflection_previous = 0
    drag_deflections = [0]
    drag_deflection_previous = 0
    for n in range(1,number_of_wing_segments):
        local_deflection = deflection_previous + (deflection_angles[n]+deflection_angles[n-1])/2*(spanwise_locations[n]-spanwise_locations[n-1])
        deflections.append(local_deflection)
        deflection_previous = local_deflection
        local_drag_deflection = drag_deflection_previous + (drag_deflection_angles[n]+drag_deflection_angles[n-1])/2*(spanwise_locations[n]-spanwise_locations[n-1])
        drag_deflections.append(local_drag_deflection)
        drag_deflection_previous = local_drag_deflection
    return deflections, drag_deflections

deflections, drag_deflections = local_deflection(spanwise_locations, deflection_angles, drag_deflection_angles)
#plt.plot(spanwise_locations, deflections)
#plt.plot(spanwise_locations, drag_deflections)