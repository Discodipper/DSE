# -*- coding: utf-8 -*-
"""
Created on Mon Jan 06 13:33:76 2020

@author: Mart
"""

import matplotlib.pyplot as plt
import scipy as sp
import copy

wing_box_chord_root = 0.5 #m
wing_box_chord_tip = 0.4 #m
wing_box_height_root = 0.2 #m
wing_box_height_tip = 0.15 #m
number_of_wing_segments = 100 #-
span = 30 #m
lift = 42000 #kN
load_factor = 4 #-
wing_weight = 10000 #N
wing_surface_area = 40 #m^2
chord_root = 1 #m
chord_tip = 0.8 #m


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

wing_box_segmentation_outcomes = wing_box_segmentation(wing_box_chord_root, wing_box_chord_tip, wing_box_height_root, wing_box_height_tip, number_of_wing_segments, span)
#plt.plot(wing_box_segmentation_outcomes[1], wing_box_segmentation_outcomes[2])
#plt.show()

def wing_segmentation(chord_root, chord_tip, number_of_wing_segments, span):
    chord_lengths = []
    spanwise_locations = []
    for i in range(0,number_of_wing_segments):
        local_chord = chord_root - i*(chord_root-chord_tip)/number_of_wing_segments
        chord_lengths.append(local_chord)
        local_spanwise_location = i*span/2/number_of_wing_segments
        spanwise_locations.append(local_spanwise_location)
    return chord_lengths,spanwise_locations

def determine_Kq(lift, load_factor, wing_weight, wing_surface_area):
    Kq = (lift - load_factor*wing_weight)/wing_surface_area
    return Kq

def load_distribution_over_wing(Kq, chord_lengths):
    local_load_distribution = []
    for i in chord_lengths:
        local_load = Kq*i
        local_load_distribution.append(local_load)
    return local_load_distribution

def shear_force_wing(local_load_distribution, spanwise_locations):
    local_shear_distribution = [0]
    shear_force_previous = 0
    spanwise_locations.reverse()
    local_load_distribution.reverse()
    n = 1
    for i in spanwise_locations[1:]:
        local_shear = shear_force_previous - (local_load_distribution[n]+local_load_distribution[n-1]/2*(i-spanwise_locations[n-1]))
        n = n + 1
        shear_force_previous = local_shear
        local_shear_distribution.append(local_shear)
    spanwise_locations.reverse()
    local_load_distribution.reverse()
    local_shear_distribution.reverse()
    return local_shear_distribution
    
def bending_moment_distribution(local_shear_distribution, spanwise_locations):
    bending_moment_distribution = [0]
    bending_moment_previous = 0
    spanwise_locations.reverse()