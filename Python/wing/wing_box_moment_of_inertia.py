# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 11:06:19 2019

@author: thoma
"""
import pathlib
import scipy as sp
import copy

pth = pathlib.Path('tether_main')
corner_stringer_area = 0.005     #m^2
standard_stringer_area = 0.003   #m^2
n_stringers_top = 4
n_stringers_bottom = 4
chord_length = 0.5  #m
chord_height = 0.2  #m

def determine_stringer_spacing(length, number_of_strings):
    stringer_spacing = length / number_of_strings
    return stringer_spacing


def boom_position_calculator(n_stringers, stringer_spacing, z_position, booms):
    for i in range(n_stringers):
        x_stringer = (i+1) * stringer_spacing
        z_stringer = chord_height
        single_stringer = sp.array([standard_stringer_area, x_stringer, z_stringer])
        booms = sp.vstack((booms,single_stringer))
    return booms

def get_stringer_position_matrix(stringer_spacing_top, stringer_spacing_bottom):
    #array to track boom parameters  [A, x, z]
    booms = sp.array([
        [corner_stringer_area, 0, 0],
        [corner_stringer_area, 0, chord_height],
        [corner_stringer_area, chord_length, 0],
        [corner_stringer_area, chord_length, chord_height]
        ])
    n_stringers_top = int(chord_length // stringer_spacing_top)
    n_stringers_bottom = int(chord_length // stringer_spacing_bottom)
    booms = boom_position_calculator(n_stringers_top, stringer_spacing_top, chord_height, booms)
    booms = boom_position_calculator(n_stringers_bottom, stringer_spacing_bottom, 0, booms)
    return booms

def neutral_position_calculator(booms):
    area_distance_x = 0 
    area_distance_z = 0
    total_area = sp.sum(booms[:,0])
    for i in range(len(booms)):
        area_distance_x += booms[i, 0] * booms[i, 1]
        area_distance_z += booms[i, 0] * booms[i, 2]
    neutral_x_position = area_distance_x / total_area
    neutral_z_position = area_distance_z / total_area
    return neutral_x_position, neutral_z_position 

def correct_for_neutral_axis(booms, neutral_x_position, neutral_z_position):
    booms_corrected_for_neutral_axis= copy.copy(booms)
    for i in range(len(booms_corrected_for_neutral_axis)):
        booms_corrected_for_neutral_axis[i,1] -= neutral_x_position
        booms_corrected_for_neutral_axis[i,2] -= neutral_z_position
    return booms_corrected_for_neutral_axis

def moment_of_inertia_boom_method():
    
    stringer_spacing_top = determine_stringer_spacing(chord_length, n_stringers_top)
    stringer_spacing_bottom = determine_stringer_spacing(chord_length, n_stringers_bottom)
    moment_of_inertia_xx = 0 
    moment_of_inertia_zz = 0 
    moment_of_inertia_xz = 0 
    
    booms = get_stringer_position_matrix(stringer_spacing_top, stringer_spacing_bottom)
    neutral_x_position, neutral_z_position = neutral_position_calculator(booms)
    booms_corrected_for_neutral_axis = correct_for_neutral_axis(booms, neutral_x_position, neutral_z_position)
    
    for i in range(len(booms_corrected_for_neutral_axis)):
        moment_of_inertia_xx += booms_corrected_for_neutral_axis[i,0]  * (booms_corrected_for_neutral_axis[i,1] ** 2)
        moment_of_inertia_zz += booms_corrected_for_neutral_axis[i,0] * (booms_corrected_for_neutral_axis[i,2] ** 2)
    return moment_of_inertia_xx, moment_of_inertia_zz, moment_of_inertia_xz

moment_of_inertia_xx, moment_of_inertia_zz, moment_of_inertia_xz = moment_of_inertia_boom_method()