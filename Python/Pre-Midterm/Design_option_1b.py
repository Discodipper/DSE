#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 09:35:43 2019

@author: pranav
"""

from math import *
import matplotlib.pyplot as plt
import numpy as np
import time
start = time.time()

altitude = np.arange(2000, 16500, 500) #M
balloon_count = np.arange(1, 101, 1)
Delft = 8.76*10**6 #W
P_needed = 10* 10**6
cable_eff = Delft/P_needed
g = 9.80665
SP_density = 18.3/(1.65*0.992) #kg/m^2
power_density = 305/(1.65*0.992) #W/m^2
operating_angle = 75*pi/180 #rad
rho_e = 2.72 * 10**(-8) #ohm * m
aluminium_density = 2700 #kg/m^3 
rho_h = 0.083777

m_cable_lst = []
A_solarpanel_lst = []
balloon_lst = []
P_out_lst = []
alt_lst = []
volume_lst = []
cable_length_lst = []
I_lst = []
W_sp_lst = []
tot_eff_lst = []
t_cable_lst = []
P_density_WRTground_lst = []
spacing_lst = []

"""This is the % of radiation increase for 500m increase."""
radiation_lst = [4.6942214931527735, 5.755366973179216, 6.774004597073443, 7.751329102742852, 8.688516290301822, 9.586723084204786, 10.447087595976798, 11.270729187555142, 12.058748535256372, 12.812227694383278, 13.532230164486982, 14.219800955300181, 14.875966653357649, 15.501735489321305, 16.098097406027478, 16.767416608916115, 17.407330073079443, 17.998751436019802, 18.54535538816441, 19.050538155007143, 19.517438598945684, 19.948957722039346, 20.347776690864386, 20.716373495461028, 21.057038345879562, 23, 23, 23, 23]

for alt in altitude:
    for balloon_number in balloon_count:
        i = (alt/500)-4
        radiation = radiation_lst[int(i)] #[-]
        solar_eff = (radiation /100)+1 #[-]
        
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
        frac_diff = delta_T * temp_coef
        SF = 1
        tot_eff = ((base_eff + frac_diff) * solar_eff) * (5/SF)
        power_density_ = (305/(1.65*0.992*base_eff)) * tot_eff #W/m^2 power density
        A_sp = P_needed/(power_density * balloon_number) #m^2 per balloon
        
        W_sp = A_sp * SP_density * g
        
        
        
        """"Cable stuff below"""
        L_cable = alt/(np.sin(operating_angle))
        V_in = 400 #V
        I_in = P_needed/(V_in * balloon_number) #A per balloon
        
        V_out = V_in * cable_eff #V
        V_drop = V_in - V_out #V
        
        R_cable = V_drop/I_in #ohm per baloon
        
        A_cable = (rho_e * L_cable)/R_cable #m^2
        t_cable = np.sqrt(A_cable*4/pi) #m per balloon
        
        W_cable = g * A_cable * L_cable *aluminium_density #N per balloon
        m_cable = W_cable/g
        
        Lift = W_sp + W_cable #N per baloon
        volume_h = Lift/(g * (rho - rho_h)) #m^3 per balloon
        
        
        """Ground spacing"""
        max_angle = pi/2 #rad
        angle_mean = (max_angle + operating_angle)/2
        angle_diff = angle_mean - operating_angle
        d_u = L_cable/(sin(angle_mean - angle_diff) * (1/tan(angle_mean - angle_diff)) * (1/(tan(2*angle_diff))))
        P_density_WRTground = Delft/(balloon_number * d_u**2)
        
        if P_density_WRTground >= 0.7165042753074458:
            volume_lst.append(volume_h)
            m_cable_lst.append(m_cable)
            A_solarpanel_lst.append(A_sp)
            alt_lst.append(alt)
            cable_length_lst.append(L_cable)
            balloon_lst.append(balloon_number)
            W_sp_lst.append(W_sp)
            tot_eff_lst.append(tot_eff)
            t_cable_lst.append(t_cable)
            P_density_WRTground_lst.append(P_density_WRTground)
            spacing_lst.append(d_u)
    

print("Operational altitude = ", alt_lst[1], "m")
print("H2 volume needed per balloon = ", volume_lst[1], "m^3")
print("Area of solar panels per balloon = ", A_solarpanel_lst[1], "m")
print("Power of single unit = ", P_needed/balloon_lst[1], "W")
print("Number of balloons = ", balloon_lst[1])
print("Thickness of the cable = ", t_cable_lst[1], "m")
print("Spacing between balloons = ", spacing_lst[1], "m")
print("Power density with respect to ground = ", P_density_WRTground_lst[1], "W/m^2")    


end = time.time()
print("Elapsed time: ", end-start)        