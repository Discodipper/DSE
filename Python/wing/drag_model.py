# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:59:40 2020

@author: wesle
"""

import numpy as np
import matplotlib.pyplot as plt



v_speed = 88.4
apparent_windspeed = v_speed
v_speed_lst = []
cd0_total_lst = []
wingdrag_surface_area = np.arange(20, 140, 40)
wingdrag_surface_area_lst = []
total_list=[]
drag_polar_lst =[]
lift_over_drag_coefficient_lst = []
lift_coefficient_lst = []
cd0_wing_lst =[]
drag_polar_wing_lst =[]
drag_polar_class0_lst = []
drag_polar_class1_lst = []

#for parameter in arange
for surface_ref in wingdrag_surface_area: #wing reference area
    surface_ratio = 3.8
    
    equivalent_skin_friction_coefficient = 0.0045
    zeroliftdrag = equivalent_skin_friction_coefficient*surface_ratio
    
    
    cdc_wing = 0.007    
    cdc_fus = 0.08
    cdc_nac = 0.06
    cdc_tail = 0.008
    cd_misc = 0.15
    ac_wing = 60#wing reference area
    ac_fus = 2*np.pi*0.5**2#maximum frontal area of fuselage
    ac_nac = 0#maximum frontal area of nacelles
    ac_tail_h = 14.8 #reference area H tail
    ac_tail_v = 3.18 #reference area V tail
    zeroliftdrag_clean = (1 / surface_ref * (cdc_wing * ac_wing + cdc_fus * ac_fus + cdc_nac * ac_nac + cdc_tail * (ac_tail_h + ac_tail_v))) * (1 + cd_misc)
    
    
    
    # ----------- Component drag build-up method -----------
    #This method estimates the subsonic parasite drag of each aircraft component
    #calculate surface areas
    L1 = 3.0 #m
    L2 = 3.0 #m
    L3 = 12.649 - L1 - L2 #m 
    d_fus = 0.5 #m
    s_w_exp = surface_ref #exposed area wing
    s_ht_exp = 14.8 #exposed area horizontal tail
    s_vt_exp = 0.7 #exposed area vertical tail
    laminar_flow_percentage_fus = 0.25 #laminar boundary layer of this amount for fuselage
    laminar_flow_percentage_wing = 0.5 #laminar boundary layer of this amount for wing and tail
    nu = 1.27E-5 #at 0 degrees celcius, about 2km altitude
    rho = 1.007 #at 2km
    speed_of_sound = 331.30 #at 0 degrees celcius, about 2km altitude
    mach_number = apparent_windspeed / speed_of_sound
    smoothness_of_material = 0.052E-5 #1.015E-5 for camouflage paint on alu, 0.634 for smooth paint, 0.405 for production sheet metal, 0.152 for polished sheet metal, 0.052 for smooth molded composite
    
    x_over_c_maximum_thickness = 0.2715
    t_over_c_average = 0.18
    length_wing = 2.385 #m #mean aerodynamic chord
    sweep_at_maximum_thickness = 0 #degrees
    
    x_over_c_maximum_thickness_tail_h = 0.2903 #NACA 0012
    t_over_c_average_tail_h = 0.1 #NACA 0012
    length_tail_h = 1.7 #m
    sweep_at_maximum_thickness_tail_h = 0 #degrees
    
    length_fus = L1+L2+L3
    maximum_frontal_area_fus = np.pi*(d_fus/2)**2
    
    #wing
    subsonic_reynolds = np.minimum(rho*apparent_windspeed*length_wing/nu, 38.21*(length_wing/smoothness_of_material)**1.053)
    laminar_friction_coefficient_part = 1.328/np.sqrt(subsonic_reynolds)
    turbulent_friction_coefficient_part = 0.455/((np.log10(subsonic_reynolds))**2.58 * (1+0.144*mach_number**2)**0.65)
    
    flat_plate_skin_friction_coefficient_wing = laminar_flow_percentage_wing * laminar_friction_coefficient_part + (1-laminar_flow_percentage_wing)*turbulent_friction_coefficient_part #estimates component friction drag
    component_form_factor_wing = (1+0.6/(x_over_c_maximum_thickness)*t_over_c_average+100*t_over_c_average**4)*(1.34*mach_number**0.18*(np.cos(np.deg2rad(sweep_at_maximum_thickness)))**0.28) #estimates pressure drag due to viscous separation
    interference_factor_wing = 1 #estimates effect of each components on drag of other components
    s_wet_wing = 1.07*2*s_w_exp #wing wetted area
    cd0_wing = 1/surface_ref*(flat_plate_skin_friction_coefficient_wing*component_form_factor_wing*interference_factor_wing*s_wet_wing)
    
    
    #horizontal tail
    subsonic_reynolds = np.minimum(rho*apparent_windspeed*length_tail_h/nu, 38.21*(length_tail_h/smoothness_of_material)**1.053)
    laminar_friction_coefficient_part = 1.328/np.sqrt(subsonic_reynolds)
    turbulent_friction_coefficient_part = 0.455/((np.log10(subsonic_reynolds))**2.58 * (1+0.144*mach_number**2)**0.65)
    
    flat_plate_skin_friction_coefficient_tail_h = laminar_flow_percentage_wing * laminar_friction_coefficient_part + (1-laminar_flow_percentage_wing)*turbulent_friction_coefficient_part #estimates component friction drag
    component_form_factor_tail_h = (1+0.6/(x_over_c_maximum_thickness_tail_h)*t_over_c_average_tail_h+100*t_over_c_average_tail_h**4)*(1.34*mach_number**0.18*(np.cos(np.rad2deg(sweep_at_maximum_thickness_tail_h)))**0.28) #estimates pressure drag due to viscous separation
    interference_factor_tail_h = 1.08
    s_wet_tail_h = 1.05*2*s_ht_exp #horizontal tail wetted area
    cd0_tail_h = 1/surface_ref*(flat_plate_skin_friction_coefficient_tail_h*component_form_factor_tail_h*interference_factor_tail_h*s_wet_tail_h)
    
    
    #vertical tail
    subsonic_reynolds = np.minimum(rho*apparent_windspeed*length_tail_h/nu, 38.21*(length_tail_h/smoothness_of_material)**1.053)
    laminar_friction_coefficient_part = 1.328/np.sqrt(subsonic_reynolds)
    turbulent_friction_coefficient_part = 0.455/((np.log10(subsonic_reynolds))**2.58 * (1+0.144*mach_number**2)**0.65)
    
    flat_plate_skin_friction_coefficient_tail_v = laminar_flow_percentage_wing * laminar_friction_coefficient_part + (1-laminar_flow_percentage_wing)*turbulent_friction_coefficient_part #estimates component friction drag
    component_form_factor_tail_v = (1+0.6/(x_over_c_maximum_thickness_tail_h)*t_over_c_average_tail_h+100*t_over_c_average_tail_h**4)*(1.34*mach_number**0.18*(np.cos(np.rad2deg(sweep_at_maximum_thickness_tail_h)))**0.28) #estimates pressure drag due to viscous separation
    interference_factor_tail_v = 1.08
    s_wet_tail_v = 1.05*2*s_vt_exp #horizontal tail wetted area
    cd0_tail_v = 1/surface_ref*(flat_plate_skin_friction_coefficient_tail_v*component_form_factor_tail_v*interference_factor_tail_v*s_wet_tail_v)
    
    #fuselage
    subsonic_reynolds = np.minimum(rho*apparent_windspeed*length_fus/nu, 38.21*(length_fus/smoothness_of_material)**1.053)
    laminar_friction_coefficient_part = 1.328/np.sqrt(subsonic_reynolds)
    turbulent_friction_coefficient_part = 0.455/((np.log10(subsonic_reynolds))**2.58 * (1+0.144*mach_number**2)**0.65)
    
    flat_plate_skin_friction_coefficient_fus = laminar_flow_percentage_fus * laminar_friction_coefficient_part + (1-laminar_flow_percentage_fus) *turbulent_friction_coefficient_part #estimates component friction drag
    f_fus=length_fus/np.sqrt(4/np.pi*maximum_frontal_area_fus)
    component_form_factor_fus = 1+60/f_fus**3 + f_fus/400
    interference_factor_fus = 1.0
    s_wet_fus = np.pi * d_fus / 4.0 * (1.0/ (3.0*L1**(2.0)) * ((4.0*L1**(2.0) + d_fus**(2.0)/4.0)**(1.5) - (d_fus**(3.0)/8.0)) - d_fus + 4.0*L2 + 2.0 * np.sqrt(L3**(2.0) + d_fus**(2.0)/4.0))#fuselage wetted area
    cd0_fus = (1/surface_ref*(flat_plate_skin_friction_coefficient_fus*component_form_factor_fus*interference_factor_fus*s_wet_fus))
    
    #other
    surface_base = 2*0.150**2*np.pi
    fuselage_base_drag_coefficient = (0.139 + 0.419*(mach_number - 0.161)**2)/(surface_base)/surface_ref
    # fuselage_base_drag_coefficient = 0
    excrescence_and_leakage_drag = 0.075 #between 5 and 10 percent
    
    
    #final formula
    cd0_total = (cd0_wing + cd0_tail_h + cd0_tail_v + 2* cd0_fus + fuselage_base_drag_coefficient)*(1+excrescence_and_leakage_drag)
    v_speed_lst.append(apparent_windspeed)
    cd0_total_lst.append(cd0_total)
    wingdrag_surface_area_lst.append(surface_ref)
    lift_coefficients = np.arange(0, 2.3, 0.1)
    for lift_coefficient in lift_coefficients:
        lift_coefficient_lst.append(lift_coefficient)
        aspect_ratio = 12
        oswald_factor = 1.78*(1 - 0.045*aspect_ratio**0.68)-0.64
        induced_drag = lift_coefficient**2/(np.pi*aspect_ratio*oswald_factor)
        drag_polar = cd0_total + induced_drag
        drag_polar_lst.append(drag_polar)
        lift_over_drag_coefficient = lift_coefficient/drag_polar
        lift_over_drag_coefficient_lst.append(lift_over_drag_coefficient)
        drag_polar_wing = cd0_wing + induced_drag
        drag_polar_wing_lst.append(drag_polar_wing)
        drag_polar_class0 = zeroliftdrag + induced_drag
        drag_polar_class0_lst.append(drag_polar_class0)
        drag_polar_class1 = zeroliftdrag_clean + induced_drag
        drag_polar_class1_lst.append(drag_polar_class1)
    
total_list.append(cd0_total_lst)
total_list.append(v_speed_lst)
total_list.append(wingdrag_surface_area_lst)
total_list.append(drag_polar_lst)
total_list.append(lift_over_drag_coefficient_lst)
total_list_array = np.array(total_list)
print("v_speed_lst =", v_speed_lst)
print("Surface area =", wingdrag_surface_area_lst)    
print("cd0_total_lst =", cd0_total_lst)
print("drag_polar_lst ="), drag_polar_lst


from cycler import cycler
color_c = cycler('color', ['b','r','k','g'])
style_c = cycler('linestyle', ['-','--',':', '-.'])
markr_c = cycler('marker', ['', '.', 'o'])
c_cms = color_c * markr_c * style_c
c_csm = color_c * style_c * markr_c



plt.rc('axes',prop_cycle=c_csm)
drag_polar_class0_60 = plt.plot(drag_polar_class0_lst[23:46],lift_coefficient_lst[23:46], label = "Class 0: wing + fuselage + tail", color='k')
drag_polar_class1_60 = plt.plot(drag_polar_class1_lst[23:46],lift_coefficient_lst[23:46], label = "Class I: wing + fuselage + tail", color= 'r')
# drag_polar_20 = plt.plot(drag_polar_lst[:23],lift_coefficient_lst[:23], label = "$ S = 20 m^2$", linestyle = '--')
drag_polar_60 = plt.plot(drag_polar_lst[23:46],lift_coefficient_lst[23:46], label = "Class II: wing + fuselage + tail", color ='g')
# drag_polar_100 = plt.plot(drag_polar_lst[46:],lift_coefficient_lst[46:], label = "$S = 100 m^2$", linestyle = '--')
# drag_polar_wing_20 = plt.plot(drag_polar_wing_lst[:23],lift_coefficient_lst[:23], label = "$ S = 20 m^2$", linestyle = '--')
drag_polar_wing_60 = plt.plot(drag_polar_wing_lst[23:46],lift_coefficient_lst[23:46], label = "Class II: wing", color = 'b')
# drag_polar_wing_100 = plt.plot(drag_polar_wing_lst[46:],lift_coefficient_lst[46:], label = "$S = 100 m^2$", linestyle = '--')





plt.ylabel("$C_{L}$")
plt.xlabel("$C_D$")
plt.title("Drag polar aircraft and wing (analytical model)")


plt.legend()
plt.grid()
plt.show()
   # s_wet_tail_v = 1.05*2*s_vt_exp #vertical tail wetted area
    # #calculate coefficients
    # flat_plate_skin_friction_coefficient = #estimates component friction drag
    # #For the wing, tail, strut and pylon:
    # component_form_factor_wing = (1+0.6/(x_over_c_maximum_thickness)*t_over_c_average+100*t_over_c_average**4)*(1.34*mach_number**0.18*(np.cos(np.rad2deg(sweep_at_maximum_thickness)))**0.28) #estimates pressure drag due to viscous separation
    # #For the fuselage and smooth canopy:
    
    # #For the nacelle and smooth external store:
    # component_form_factor_nac = 1+0.35/f
    
    # fuselage_base_drag_coefficient = (0.139+0.419(mach_number-0.161)**2)
    # excrescence_and_leakage_drag = 0.075 #between 5 and 10 percent
    
    
    
    
    
    
    # smoothnesss_of_material = 0.052E-5 #1.015E-5 for camouflage paint on alu, 0.634 for smooth paint, 0.405 for production sheet metal, 0.152 for polished sheet metal, 0.052 for smooth molded composite
    # subsonic_reynolds = np.min(rho*V*l/nu, 38.21*(l/smoothness_of_material)**1.053)
    # laminar_friction_coefficient_part = 1.328*np.sqrt(subsonic_reynolds)
    # turbulent_friction_coefficient_part = 0.455/((np.log10(subsonic_reynolds))**2.58 * (1+0.144*mach_number**2)**0.65)
    
    # cd_misc_component = 
    # component_drag_buildup_method_cd0 = 1/surface_ref*(flat_plate_skin_friction_coefficient*component_form_factor*interference_factor*surface_wet) + cd_misc_component #estimates component friction drag
