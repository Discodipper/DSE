#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 13:35:31 2019

@author: pranav
"""



"""ISA calculations"""
def isa(alt):
    a = -0.0065
    T_0 = 288.15
    g = 9.80665
    p_0 = 101325
    rho_0 = 1.225
    
    
    T = T_0 + (a*alt)
    p = (p_0) * (T/T_0)**((-g)/(a*287))
    rho = (rho_0) * (T/T_0)**(((-g)/(a*287))-1)
    V_w = 7*10**-8*alt**2 + 0.0011*alt + 7.5388 #m/s
    
    return(T, p, rho, V_w)