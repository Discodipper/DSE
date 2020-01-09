# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 11:15:24 2020

@author: thoma
"""

import scipy as sp
import sys, os
import matplotlib.pyplot as plt
sys.path.clear()
sys.path.append(os.path.realpath('..\\wing'))
sys.path.append(os.path.realpath('..\\Flight_Performance\\Shape\\weight_and_cg_estimation'))

from shape_parameters import chord_tip, chord_root, wing_span_half, sweep_angle_rad
file_name = 'MRevE-v2.dat'


def retrieve_aerofoil_parameters(file_name):
    notepad = open(file_name, 'r')
    aerofoil_data = str.split(notepad.read())
    aerofoil_cross_section= []

# =============================================================================
#     #linear interpolation
#     aerofoil_data_upper = aerofoil_data[:len(aerofoil_data)//2 -1]
#     aerofoil_data_lower = aerofoil_data[len(aerofoil_data)//2 - 1:]
#     cambers = []
#     aerofoil_data = [float(step) for step in aerofoil_data]
#     for i in range (0, len(aerofoil_data), 2):
#         print(aerofoil_data[i-2])
#         print(aerofoil_data[i])
#         print(aerofoil_data[i] > 0.1, aerofoil_data[i-2] < 0.1)
#         if (aerofoil_data[i] > 0.1 and aerofoil_data[i-2] < 0.1):
#             camber = aerofoil_data[i-1] + (0.1 - aerofoil_data[i-2]) / (aerofoil_data[i] - aerofoil_data[i-2]) * (aerofoil_data[i-1]-aerofoil_data[i+1]) 
#             print(camber)
#             cambers = sp.append(cambers, camber)
#     
# =============================================================================
    
    #reading from the plot
    for i in range (0, len(aerofoil_data), 2):
        aerofoil_point = {
            'chord':float(aerofoil_data[i]),
            'camber':float(aerofoil_data[i + 1])
            }
        aerofoil_cross_section= sp.append(aerofoil_cross_section, aerofoil_point)
    #5plt.plot([step['chord'] for step in aerofoil_cross_section], [step['camber'] for step in aerofoil_cross_section])
    #plt.axis('equal')
    
    #as read from the plot
    top_z_position_at_10_percent_chord = 0.142
    low_z_position_at_10_percent_chord = -0.0742
    
    top_z_position_at_70_percent_chord = 0.117
    low_z_position_at_70_percent_chord = 0.02313
    wing_box_z_dimensional_properties = {
        'top_z_position_at_10_percent_chord': top_z_position_at_10_percent_chord,
        'low_z_position_at_10_percent_chord': low_z_position_at_10_percent_chord,
        'z_distance_at_10_percent_chord': top_z_position_at_10_percent_chord - low_z_position_at_10_percent_chord,
        'top_z_position_at_70_percent_chord': top_z_position_at_70_percent_chord,
        'low_z_position_at_70_percent_chord': low_z_position_at_70_percent_chord,
        'z_distance_at_70_percent_chord': top_z_position_at_70_percent_chord - low_z_position_at_70_percent_chord,
        'upper_camber_slope': top_z_position_at_70_percent_chord - top_z_position_at_10_percent_chord,
        'lower_camber_slope': low_z_position_at_70_percent_chord - low_z_position_at_10_percent_chord
        }
    return z_positions
    

def wing_box_width_calculator(y_position, chord_for_y_position):
    #track trailing edge and leading edge positions 
    LE_x_position_for_y_position = y_position * - sp.sin(sweep_angle_rad)
    TE_x_position_for_y_position = LE_x_position_for_y_position - chord_for_y_position

    #consider high-lift devices and front spacing
    HLD_LE_spacing_as_chord_percentage = 0.1
    HLD_TE_spacing_as_chord_percentage = 0.3 
    LE_wing_box_x_position = LE_x_position_for_y_position - HLD_LE_spacing_as_chord_percentage  * chord_for_y_position
    TE_wing_box_y_position = TE_x_position_for_y_position + HLD_TE_spacing_as_chord_percentage  * chord_for_y_position
    wing_box_x_distance = LE_wing_box_x_position - TE_wing_box_y_position 
    return LE_wing_box_x_position, TE_wing_box_y_position, wing_box_x_distance
    
# =============================================================================
#     #verify correctness of shape parameters   
#    LE_x_positions = sp.append(LE_x_positions, LE_x_position_for_y_position)
#    TE_x_positions =  sp.append(TE_x_positions, TE_x_position_for_y_position)
#     plt.plot(y_positions, LE_x_positions, label='LE x-position')
#     plt.plot(y_positions, TE_x_positions, label='TE x-position')    
#     plt.axis('equal')
#     plt.legend()
# =============================================================================


def wing_box_height_calculator():
    #values taken from MRevE-v2_HC
    height_at_10_percent_chord = 0.14585443 + (0.16861102 - 0.14585443) * (0.09344809 - 0.1) * (0.09344809 - 0.12877692)


def determine_wing_box_parameters(n_sections):
    #####################################################################
    #this is subject to change when the double fuselages is implemented
    chord_step_cross_section = (chord_tip- chord_root)/ n_sections
    #####################################################################
    y_step_cross_section = wing_span_half/n_sections
    TE_x_positions = []
    LE_x_positions = []
    y_positions = []
    for n in range(0, n_sections + 1) : 
        y_position = y_step_cross_section * n
        y_positions = sp.append(y_positions, y_position)
        chord_for_y_position = chord_root + chord_step_cross_section * n
        
        LE_wing_box_x_position, TE_wing_box_x_position, wing_box_x_distance = wing_box_width_calculator(y_position, chord_for_y_position)
        #wing_box_height_calculator(y_position, chord_for_y_position)
    
z_positions = retrieve_aerofoil_parameters(file_name)
determine_wing_box_parameters(1000)