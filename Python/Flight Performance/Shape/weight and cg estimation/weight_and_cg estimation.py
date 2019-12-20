import math as m
import numpy as np
# wing_aspect_ratio = 12 # [-]
# wing_area = 40 # [m^2]
# rho = 1.007 #density at 2 km [kg/m^3]
# V_wind = 10 # [m/s]
# V_glider = 30 # [m/s]
# V = V_wind + V_glider
# wing_span = m.sqrt(wing_aspect_ratio*wing_area)
# lift_coefficient = 1.5
# lift = 0.5 * rho * V**2 * wing_area * lift_coefficient
length_fuselage_dylan = 20
wing_area_dylan = 150 
wing_area = 60
g = 9.80665
ratio_from_dylan = m.sqrt(wing_area/wing_area_dylan)
length_fuselage = ratio_from_dylan * length_fuselage_dylan

#------------------------------------------- COGs --------------------------------------------------------
center_of_gravity_engines = 0*length_fuselage  # [m]
center_of_gravity_tail = 0.95*length_fuselage # [m]
center_of_gravity_fuselage = 0.5*length_fuselage # [m]
center_of_gravity_wheels_nose = 0.1* length_fuselage # [m]
center_of_gravity_wheels_main = 0.55*length_fuselage # [m]
center_of_gravity_wing = 0.42 * length_fuselage # [m]

#------------------------------------------- WEIGHTS ----------------------------------------------------
W_glider = ratio_from_dylan * 6885.2 *g # [N] 
W_tail = 0.08*W_glider# [N]
W_wing = 0.50*W_glider # [N]
W_wheels_nose = 0.02*W_glider # [N] # Weight of entire nose wheels 
W_wheels_main = 0.04*W_glider # [N] Weight of entire main wheels 
W_fuselage = 0.20 * W_glider # [N] Weight of both fuselages combined 
W_engines = 0.16*W_glider # [N] Weight of both engines combined 

# import stuff here
def Center_of_Gravity(W_tail, W_wheels_main, W_wheels_nose, W_fuselage, W_engines, W_wing, \
                      center_of_gravity_tail, center_of_gravity_wheels_main, center_of_gravity_wheels_nose, center_of_gravity_fuselage, center_of_gravity_engines, center_of_gravity_wing):
    W_fuselage_group = W_tail + W_fuselage + W_wheels_nose + W_engines
    W_wing_group = W_wing  + W_wheels_main 
    center_of_gravity_fuselage_group = (W_tail*center_of_gravity_tail + W_fuselage*center_of_gravity_fuselage + \
                                        W_wheels_nose*center_of_gravity_wheels_nose + W_tail*center_of_gravity_tail + W_engines * center_of_gravity_engines) / W_fuselage_group # in [m]
    center_of_gravity_wing_group = (W_wheels_main * center_of_gravity_wheels_main + W_wing*center_of_gravity_wing) / W_wing_group # in [m]
    center_of_gravity =  (W_wing_group * center_of_gravity_wing_group + W_fuselage_group * center_of_gravity_fuselage_group ) / W_glider  # in [m] from the nose
    if center_of_gravity > center_of_gravity_wheels_main:
        print("Error: main wheel cg in front of total cg")
    center_of_gravity_ratio = center_of_gravity/length_fuselage # cg as a percentage of length fuselage taken from nose to aft
    return(center_of_gravity, center_of_gravity_ratio)