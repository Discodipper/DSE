# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 12:06:16 2020

@author: Mart

This script is built to use the loads on the wing to calculate the bending stress in the wing box
"""

def bending_stress_wing_box(bending_moment_x, bending_moment_y, moment_of_inertia_x, moment_of_inertia_y, x_coordinate_from_neutral_axis, y_coordinate_from_neutral_axis):
    bending_stress = (moment_of_inertia_y*bending_moment_x*y_coordinate_from_neutral_axis+bending_moment_y*moment_of_inertia_x*x_coordinate_from_neutral_axis)/moment_of_inertia_x*moment_of_inertia_y
    return bending_stress

def moments_of_inertia(height, chord):
    moment_of_inertia_y = 1/12*chord^3*height
    moment_of_inertia_x = 1/12*height*chord^3
    moment_of_inertia_z = chord*height^3*(16/3-3.36*height/chord*(1-height^4/12/chord^4))
    return moment_of_inertia_x, moment_of_inertia_y, moment_of_inertia_z

def torsional_stress(torsion, radius, moment_of_inertia_z ):
    torsional_stress = torsion*radius/moment_of_inertia_z
    return torsional_stress

def thickness_required_bending(bending_moment_y, bending_moment_x, yield_stress_plate, safety_factor, wing_box_chord_height, wing_box_chord_width):
    thickness_required_bending = 0.001
    running = True
    while running == True:
        bending_stress = bending_moment_x*0.5*wing_box_chord_height/2/(thickness_required_bending**(3)/12*wing_box_chord_width+wing_box_chord_width*thickness_required_bending*(0.5*wing_box_chord_height-thickness_required_bending/2)**(2))+bending_moment_y*0.5*wing_box_chord_width/2/(wing_box_chord_width**(3)/12*thickness_required_bending)
        thickness_required_bending = thickness_required_bending + 0.001
        if bending_stress < yield_stress_plate*safety_factor: 
            running = False
    return thickness_required_bending                  

#thickness = thickness_required_bending(3000000.0, 500000.0, 2500000.0, 3.0, 0.5, 3.0)
#print(thickness)