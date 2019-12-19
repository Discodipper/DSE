#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:37:38 2019

@author: pranav
"""

from math import *
import matplotlib.pyplot as plt
import numpy as np
import time
start = time.time()

altitude = np.arange(2000, 16500, 500) #M
volume = np.arange(2000000, 5000000, 10000) #m^3
g = 9.80665 
rho_h = 0.0837765634 #hydrogen in kg/m3
voltage = 130 #V
m_cable_lst = []
m_solarpanel_lst = []
A_solarpanel_lst = []
eff_lst = []
P_out_lst = []
alt_lst = []
V_lst = []
cable_thickness_lst = []
I_lst = []
Tot_area_needed = []
delta_T_lst = []
t_c = []
W_cables = []


rho_e = 1.72*10**(-8)
I = 9.91 #A
V_in = 31 #V  
R_cable = (0.2*V_in)/I

nu = 0.8
P_req = 8700000 #W
P_req2 = 8000000 #W


radiation_lst = [0,1.243505181776943, 2.439532955387986, 3.5893544197450122, 4.6942214931527735, 5.755366973179216, 6.774004597073443, 7.751329102742852, 8.688516290301822, 9.586723084204786, 10.447087595976798, 11.270729187555142, 12.058748535256372, 12.812227694383278, 13.532230164486982, 14.219800955300181, 14.875966653357649, 15.501735489321305, 16.098097406027478, 16.767416608916115, 17.407330073079443, 17.998751436019802, 18.54535538816441, 19.050538155007143, 19.517438598945684, 19.948957722039346, 20.347776690864386, 20.716373495461028, 21.057038345879562]
for alt in altitude:
    #----------------------- efficiency -----------------------------------#
    i = (alt/500)-4    
    radiation = radiation_lst[int(i)] #[-]
    solar_eff = (radiation /100)+1 #[-]
    operating_angle = 75*pi/180 #rad
#        """ISA calculations"""
#        a = -0.0065
#        T_0 = 288.15
#        g = 9.80665
#        p_0 = 101325
#        rho_0 = 1.225
#        T = T_0 + (a*alt)
#        p = (p_0) * (T/T_0)**((-g)/(a*287))
#        rho = (rho_0) * (T/T_0)**(((-g)/(a*287))-1)
#        rho_h = 0.083777 #kg/m^3
    if (0<  alt <= 11000):
        p_0 = 101325 #Pa
        T_0 = 288.15 #K
        h_0 = 0 #m
        R_air = 287
        a_lapse = -0.0065
        T = T_0 + a_lapse * (alt-h_0) #Temperature for gradient layer
        p = p_0 * (T/T_0)**(-g/(a_lapse*R_air))
        rho = p/(R_air*T)
                
    elif (11000 <alt <= 20000):
        p_0 = 22700 #Pa
        h_0 = 11000
        T = 216.8 #K
        rho_0 = 0.364805
        p = p_0 * e**((-g/(R_air*T))*(alt-h_0))
        rho = rho_0 * (p/p_0)
    

    
    temp_coef = -0.004
    base_eff = 0.18
    delta_T = T - 298.15
    delta_T_lst.append(delta_T)
    frac_diff = delta_T * temp_coef
    SF = 1
    efficiency = ((base_eff + frac_diff) * solar_eff) * (5/SF)
    Area_solp = 1.65 #Area of one solar panel
    P_out = (305/base_eff)* efficiency  #W output of one solar panel
    P_out_sqrm = P_out / Area_solp #Power output per squared meter
    Tot_area_needed.append(P_req / P_out_sqrm) #m^2
    
    L_cable = alt/(np.sin(operating_angle))
    t_c = 2 * sqrt(I*rho_e*(L_cable/(pi*V_in*(1-nu))))
    radius_cable = t_c/2 #np.sqrt(((resistivity * L_cable)/R_cable)/pi)
    copper_density = 8890 #kg/m^3
    m_cable = pi * ((radius_cable)**2) *L_cable * copper_density
    W_cables = m_cable * g
    
    
    
    L = V * g * (rho - rho_h)
    W_solp = L-W_cables 
    
    
    
    
    
    
#        nu = 0.8
#        t = 2*np.sqrt(I*resistivity *L_cable/((1-nu)*pi*voltage))
        
        
        
        
#
#        
#        
#        

#        
#        
#        W_solarpanel = L - W_cable
#        
#        
#        if W_solarpanel > 0:
#            rho_solarpanel = 0.135/(0.0027*0.008)
#            m_solarpanel = W_solarpanel/g #kg
#            A_solarpanel = (m_solarpanel/rho_solarpanel)/0.008 #m^2
#            number_solarpanel = A_solarpanel/0.0027
#            current_solarpanel = (17.4/1000)*(A_solarpanel*10000)
#            "This is the surface area of the solar panel sub system that can be lifted by the balloon."
#            
#        
#            
#            P_out_lst.append(P_out)
#            
#            cable_thickness_lst.append(t)
#            m_cable_lst.append(m_cable)
#            m_solarpanel_lst.append(m_solarpanel)
#            A_solarpanel_lst.append(A_solarpanel)
#            eff_lst.append(efficiency)
#            alt_lst.append(alt)
#            V_lst.append(V)
#            I_lst.append(I)
#    
#    
    
    
    
        




end = time.time()
print("Elapsed time: ", end-start)