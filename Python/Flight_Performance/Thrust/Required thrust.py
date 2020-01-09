# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 16:47:26 2020

@author: ErikJan
"""

#imports
from ISA_calculator import isa
import sys, os
sys.path.append(os.path.realpath(

def take_off_thrust(surface_area__main_wing, coefficients_CL, coefficients_CL_CD, glider_mass):
    'surface_area__main_wing = single value in m^2, coefficients_CL = list of coefficients of length n of nth degree polynomial in order from alpha^0 to alpha^n, coefficients_CL_CD = list of coefficients of length n of nth degree polynomial in order from CL^0 to CL^n, glider_mass = single number in kg'
    