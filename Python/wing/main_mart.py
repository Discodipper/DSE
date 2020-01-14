# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 11:52:22 2020

@author: Mart

This script describes the main part of defining material and properties of the wing box
"""


import matplotlib.pyplot as plt
import scipy as sp
import copy
import sys, os
sys.path.clear()
sys.path.append(os.path.realpath('..\\wing'))
sys.path.append(os.path.realpath('..\\Flight_Performance\\Shape\\weight_and_cg_estimation'))

from shape_parameters import sweep_angle_deg, wing_area, wing_span, chord_root, chord_tip
from wing_box_stress_calculations import bending_stress_wing_box, thickness_required_bending
from wing_loads_and_deflection_calculations import wing_box_segmentation, wing_segmentation

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

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

chord_lengths, spanwise_locations = wing_segmentation(chord_root, chord_tip, number_of_wing_segments, span)
plt.plot(spanwise_locations, chord_lengths)

#plate_thickness = thickness_required_bending(bending_moment_y, bending_moment_x, yield_stress_plate, safety_factor, wing_box_chord_height, wing_box_chord_width)
