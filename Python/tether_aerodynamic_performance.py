# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 11:57:40 2019

@author: ErikJan
"""


# imports
from ISA_calculator import isa
from numpy import tan, pi, arange, sin, cos

# drag of cables underneath 'Y'-connection calculators
def drag_stationary_part_tether_reeled_in(horizontal_distance_tether_2): #in Xw direction
    # test inputs
    diameter_tether_2 = .008 #m # from Caro and TA
    phi = 0 #degree # from Pranav and Xander
    
    # inputs
    Cd = 1 #median of Schmel
    mesh_length_2 = .100 #mesh length in x direction
    operation_angle = 25 #degree # output from Xander and Pranav
    
    Drag_tether_2 = []
    
    for x in arange(0,horizontal_distance_tether_2+mesh_length_2,mesh_length_2):
        h = x*tan(operation_angle/180*pi) # from Caro and TA
        Drag_tether_2_segment = .5*Cd*isa(h)[2]*isa(h)[3]**2*diameter_tether_2*mesh_length_2*sin(phi/180*pi) #worst case scenario
        Drag_tether_2.append(Drag_tether_2_segment)
        
    return(sum(Drag_tether_2))

def drag_stationary_part_tether_reeled_out(horizontal_distance_tether_1,horizontal_distance_tether_2): #in Xw direction
    # test inputs
    diameter_tether_1 = .016 #m # from Caro and TA
    diameter_tether_2 = .008 #m # from Caro and TA
    phi = 25 #degree # from Pranav and Xander
    
    # inputs
    Cd = 1 #median of Schmel
    mesh_length_1 = .100
    mesh_length_2 = .100
    
    operation_angle = 25 #degree # output from Xander and Pranav
    
    Drag_tether_1 = []
    Drag_tether_2 = []
    
    for x in arange(0,horizontal_distance_tether_1+mesh_length_1,mesh_length_1):
        h = x*tan(operation_angle/180*pi) # from Caro and TA
        Drag_tether_1_segment = .5*Cd*isa(h)[2]*isa(h)[3]**2*diameter_tether_1*mesh_length_1*sin(phi/180*pi)
        Drag_tether_1.append(Drag_tether_1_segment)
    
    for x in arange(horizontal_distance_tether_1,horizontal_distance_tether_2+mesh_length_2,mesh_length_2):
        h = x*tan(operation_angle/180*pi) # from Caro and TA
        Drag_tether_2_segment = .5*Cd*isa(h)[2]*isa(h)[3]**2*diameter_tether_2*mesh_length_2*sin(phi/180*pi)
        Drag_tether_2.append(Drag_tether_2_segment)
        
    return(sum(Drag_tether_1)+sum(Drag_tether_2))

def Y_single_tether_drag(Y_tether_length,Y_tether_diameter,h,angle_between_Y_tethers):
    # inputs
    Vkr = 30 # from Pranav and Xander
    Vkt = 30 # from Pranav and Xander
    Cd = 1 # median of Schmel
    phi = 25 # degree # from Pranav and Xander
    phi = phi/180*pi
    theta = 35 # degree # from Pranav and Xander
    theta = theta/180*pi
    gamma = angle_between_Y_tethers/360*pi
    Y = Y_tether_length
    
    # va[zero[x,y,z],...,twohundredseventy[x,y,x]]
    va = []
    for xhi in arange(0,2*pi,.5*pi):
        vax = sin(theta)*cos(phi)*isa(h)[3]-Vkr
        vay = cos(theta)*cos(phi)*isa(h)[3]-cos(xhi)*Vkt
        vaz = -sin(phi)*isa(h)[3]-sin(xhi)*Vkt
        va.append([vax,vay,vaz])
    
    if phi/.25/pi <= 1:
        if theta/.25/pi <= 1:
            if gamma/.25/pi <= 1:
                Y_tether_area = [[Y,Y,0],[Y,Y,0],[Y,Y,0],[Y,Y,0]]
            else:
                Y_tether_area = [[0,Y,Y],[Y,0,Y],[0,Y,Y],[Y,0,Y]]
        else:
            if gamma/.25/pi <= 1:
                Y_tether_area = [[0,Y,Y],[0,Y,Y],[0,Y,Y],[0,Y,Y]]
            else:
                Y_tether_area = [[Y,Y,0],[Y,0,Y],[Y,Y,0],[Y,0,Y]]
    else:
        if theta/.25/pi <= 1:
            if gamma/.25/pi <= 1:
                Y_tether_area = [[Y,Y,0],[Y,Y,0],[Y,Y,0],[Y,Y,0]]
            else:
                Y_tether_area = [[Y,0,Y],[0,Y,Y],[Y,0,Y],[0,Y,Y]]
        else:
            if gamma/.25/pi <= 1:
                Y_tether_area = [[Y,0,Y],[Y,0,Y],[Y,0,Y],[Y,0,Y]]
            else:
                Y_tether_area = [[Y,Y,0],[0,Y,Y],[Y,Y,0],[0,Y,Y]]
    
    # calculate drag in x, y and z direction for xhi is zero, ..., twohundredseventy
    drag_xyz_zero = []
    drag_xyz_ninety = []
    drag_xyz_hundredeighty = []
    drag_xyz_twohundredseventy = []
    
    for i in range(0,3):
        zero = .5*Cd*isa(h)[2]*(.5*va[0][i])**2*Y_tether_diameter*Y_tether_area[0][i]
        drag_xyz_zero.append(zero)
        ninety = .5*Cd*isa(h)[2]*(.5*va[1][i])**2*Y_tether_diameter*Y_tether_area[1][i]
        drag_xyz_ninety.append(ninety)
        hundredeighty = .5*Cd*isa(h)[2]*(.5*va[2][i])**2*Y_tether_diameter*Y_tether_area[2][i]
        drag_xyz_hundredeighty.append(hundredeighty)
        twohundredseventy = .5*Cd*isa(h)[2]*(.5*va[3][i])**2*Y_tether_diameter*Y_tether_area[3][i]
        drag_xyz_twohundredseventy.append(twohundredseventy)
    
    return(drag_xyz_zero,drag_xyz_ninety,drag_xyz_hundredeighty,drag_xyz_twohundredseventy)


