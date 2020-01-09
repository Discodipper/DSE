# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 12:06:16 2020

@author: Mart

This script is built to use the loads on the wing to calculate the bending stress in the wing box
"""

def bending_stress_wing_box(bending_moment_x, bending_moment_y, moment_of_inertia_x, moment_of_inertia_y, x_coordinate_from_neutral_axis, y_coordinate_from_neutral_axis):
    bending_stress = (moment_of_inertia_y*bending_moment_x*y_coordinate_from_neutral_axis+bending_moment_y*moment_of_inertia_x*x_coordinate_from_neutral_axis)/moment_of_inertia_x*moment_of_inertia_y
    return bending_stress

def torsional_stress(torsion, radius, polar_moment_of_inertia, )