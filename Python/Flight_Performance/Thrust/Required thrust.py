# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 16:47:26 2020

@author: ErikJan
"""

#imports
import sys, os
#!!!!! DO NOT REMOVE UNLESS YOU KNOW HOW TO OPTIMISE!!!!
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath('...'))))

from ISA_calculator import isa
from numpy import arange
from math import sqrt, pi

g = 9.80665


# def take_off_thrust(surface_area_main_wing, coefficients_CL, coefficients_CL_CD, glider_mass, h, alpha_max):
#     'surface_area__main_wing = single value in m^2, coefficients_CL = list of coefficients of length 7 of 6th degree polynomial in order from alpha^0 to alpha^6, coefficients_CL_CD = list of coefficients of length 7 of 6th degree polynomial in order from CL^0 to CL^6, glider_mass = single number in kg, h = operating altitude in m, alpha_max = integer value of maximum angle of attack in degree'
#     # assume that high lift devices are deployed until operating altitude is reached.
        
#     # CL and CD w.r.t. angle of attack
#     CL_alpha = []
#     CD_CL_alpha = []
#     for alpha in arange(0,alpha_max,.1):
#         CL = coefficients_CL[0]+coefficients_CL[1]*alpha+coefficients_CL[2]*alpha**2+coefficients_CL[3]*alpha**3+coefficients_CL[4]*alpha**4+coefficients_CL[5]*alpha**5+coefficients_CL[6]*alpha**6
#         CL_alpha.append(CL)
#         CD = coefficients_CL_CD[0]+coefficients_CL_CD[1]*CL+coefficients_CL_CD[2]*CL**2+coefficients_CL_CD[3]*CL**3+coefficients_CL_CD[4]*CL**4+coefficients_CL_CD[5]*CL**5+coefficients_CL_CD[6]*CL**6
#         CD_CL_alpha.append(CD)
    
#     # stall speed
#     CL_max =  max(CL_alpha)
#     V_stall = sqrt(2*g*glider_mass/isa(h)[2]/surface_area_main_wing/CL_max)
    
#     # lift and drag calculations at ground level assuming V is constant and bigger than stall speed at operating altitude
#     L_0_alpha_V = []
#     #L_h_alpha_V = []
#     D_0_alpha_V = []
#     #D_h_alpha_V = []
#     for i in range(len(CL_alpha)):
#         for V in range(int(V_stall),100):
#             L_0 = .5*isa(0)[2]*V**2*surface_area_main_wing*CL_alpha[i]
#             #L_h = .5*isa(h)[2]*V**2*surface_area_main_wing*CL_alpha[i] 
#             D_0 = .5*isa(0)[2]*V**2*surface_area_main_wing*CD_CL_alpha[i] 
#             #D_h = .5*isa(h)[2]*V**2*surface_area_main_wing*CD_alpha[i]
#             L_0_alpha_V.append(L_0)
#             #L_h_alpha_V.append(L_h)
#             D_0_alpha_V.append(D_0)
#             #D_h_alpha_V.append(D_h)
    
#     # minimal possible drag and thrust assuming thrust is 6 times drag of main wing
#     D_min = min(D_0_alpha_V)
#     required_thrust = 6*D_min

#     return(required_thrust)

# engine selection
def engine(required_thrust,V_take_off):
    propeller_efficiency = .64
    C_T = 8/9#.065 #.1064,.1004,.0924,.0846,.0696
    n = 6000 # RPM from literature test data
    diameter_single_propeller = sqrt(sqrt(.5*required_thrust/isa(0)[2]/(n*2*pi)**2/C_T))
    required_engine_power = 1/propeller_efficiency*required_thrust*V_take_off
    total_engine_mass = required_engine_power/13000
    engine_diameter = .45 #m
    engine_length = .12 #m
    return(required_engine_power,total_engine_mass, engine_diameter,engine_length,diameter_single_propeller)
    
#!!!!! DO NOT REMOVE UNLESS YOU KNOW HOW TO OPTIMISE!!!!