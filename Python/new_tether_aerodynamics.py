# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:05:23 2020

@author: ErikJan
"""

from ISA_calculator import isa
from numpy import sin, cos, pi, arange

# general inputs
C_D_tether = 1
angle_ground_glider = 20 # degree
Y_tether_angle = 60 # degree
distance_gliders = 2*467 # m
V_glider = 87.7 # m/s
surface_area_glider = 60 # m^2
altitude_glider = 3000 # m

def Y_single_tether_drag(altitude_glider,Y_tether_diameter,V_glider,surface_area_glider):
    mesh_length = Y_tether_diameter
    Y_single_tether_length = distance_gliders/2/sin(Y_tether_angle/2/180*pi)
    D = 0
    for i in arange(0,Y_single_tether_length,mesh_length):
        L = i
        D = D+.5*C_D_tether*isa(altitude_glider)[2]*(L/Y_single_tether_length*V_glider)**2*mesh_length*Y_tether_diameter  
    C_D_gebeund = 2*D/isa(altitude_glider)[2]/V_glider**2/surface_area_glider
    return(D,C_D_gebeund)

def Y_single_tether_drag_iterative(Y_tether_diameter, distance_gliders):
    mesh_length = Y_tether_diameter
    Y_single_tether_length = distance_gliders/2/sin(Y_tether_angle/2/180*pi)
    D = 0
    for i in arange(0,Y_single_tether_length,mesh_length):
        L = i
        D = D+.5*C_D_tether*isa(altitude_glider)[2]*(L/Y_single_tether_length*V_glider)**2*mesh_length*Y_tether_diameter  
    C_D_gebeund = 2*D/isa(altitude_glider)[2]/V_glider**2/surface_area_glider
    return(D,C_D_gebeund)

def stationary_tether_drag(altitude_glider,stationary_tether_diameter,C_length_sag_cable):
    mesh_length = stationary_tether_diameter
    stationary_tether_length = altitude_glider/sin(angle_ground_glider/180*pi)
    D = 0
    for i in arange(0,stationary_tether_length,mesh_length):
        D = D+.5*C_D_tether*isa(i*sin(35/180*pi))[2]*(isa(i*sin(35/180*pi))[3])**2*mesh_length*stationary_tether_diameter
    
    D_X_W = C_length_sag_cable/stationary_tether_length*D*sin(angle_ground_glider/180*pi)*sin(angle_ground_glider/180*pi)
    D_Z_W = C_length_sag_cable/stationary_tether_length*D*sin(angle_ground_glider/180*pi)*cos(angle_ground_glider/180*pi)
    
    return(D_X_W,D_Z_W)

def moving_cable_drag(altitude_glider, tether_diameter, V_glider,C_length_sag_cable):
    mesh_length = tether_diameter
    tether_length = altitude_glider/sin(angle_ground_glider/180*pi)
    D = 0
    for i in arange(0,tether_length,mesh_length):
        D = D+.5*C_D_tether*isa(i*sin(35/180*pi))[2]*(i/tether_length*V_glider)**2*mesh_length*tether_diameter
    
    D = C_length_sag_cable/tether_length*D
    
    return(D)