import math as m
g = 9.80665
sweep_angle_deg = 20 # degrees at LE 
sweep_angle_rad = (1*m.pi)/180 # sweep angle at LE in radians
taper_ratio = 0.4 # to make lift distribution elliptical 
aspect_ratio = 12 # [-]
wing_area = 60 # m^2
wing_span = m.sqrt(wing_area * aspect_ratio) # [m]
chord_root = wing_area/(wing_span*0.7) # [m]
chord_tip = taper_ratio * chord_root # [m]
chord_mean_aerodynamic = (2/3)*chord_root * ((1+taper_ratio + \
taper_ratio**2)/(1+taper_ratio)) # [m]
mean_aerodynamic_y = (wing_span/6) * ((1+2*taper_ratio)/(1+taper_ratio)) # [m] # taken from middle 
mean_aerodynamic_x = mean_aerodynamic_y * m.tan(sweep_angle_rad) # [m] TAKEN FROM WHERE?????
quarter_chord_line_angle_rad = m.atan(m.tan(sweep_angle_rad) + (chord_root/(2*wing_span)) * (taper_ratio-1))
aspect_ratio_limit = 17.7 * (2-taper_ratio) * m.exp(-0.043 * quarter_chord_line_angle_rad)

import ISA_calculator as isa_calc 
gamma = 1.4
R = 287.00
alt = 2000
speed_apparent = 99 # Apparent speed [m/s] SUBJECTED TO CHANGE!!! import from pranav's code 0
T, p, rho, V_w = isa_calc.isa(alt)


def mach_number(speed_apparent,T):
    gamma = 1.4
    R = 287.00
    speed_of_sound = m.sqrt(gamma*R*T)
    mach = speed_apparent/speed_of_sound # [-]
    return (mach)

beta = m.sqrt(1-mach**2) # compressibility factor 
half_chord_line_angle_rad = m.atan(m.tan(sweep_angle_rad) - 0.5 * ((2*chord_root)/wing_span)*(1-taper_ratio))
eta = 0.95 # airfoil efficiency factor 0.95 for good approximation
def CL_alpha(aspect_ratio,beta,eta,half_chord_line_angle_rad):
        
    lift_coefficient_alpha_rad = (2*m.pi*aspect_ratio) / (2+m.sqrt(4+ ((aspect_ratio*beta/eta)**2) * (1+(m.tan(half_chord_line_angle_rad)**2)/beta**2))) # [1/rad]
    lift_coefficient_alpha_deg = lift_coefficient_alpha_rad *m.pi /180 # lift coefficient for entire wing 
    return (lift_coefficient_alpha_rad,lift_coefficient_alpha_deg)

max_camber = 6              # [%]
W_glider_dylan = 6885.2 *g
W_glider = 42703.869701564894
dynamic_viscosity = 1.729 * 10**(-5)
reynolds_number = ((rho * speed_apparent * chord_mean_aerodynamic)/dynamic_viscosity)
q = 0.5*rho*speed_apparent**2
oswald = 0.95
CD0 = 0.01

lift_coefficient = 18*(1/q)*(W_glider/wing_area)
CD = CD0 + (lift_coefficient**2/(m.pi*aspect_ratio*oswald))
lift_over_drag_ratio = lift_coefficient/CD






