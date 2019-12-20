import math as m
sweep_angle_deg = 1 # degrees at LE 
sweep_angle_rad = (1*m.pi)/180 # sweep angle at LE in radians
taper_ratio = 0.4 # to make lift distribution elliptical 
aspect_ratio = 12 # [-]
wing_area = 40 # m^2
wing_span = m.sqrt(wing_area * aspect_ratio) # [m]
chord_root = wing_area/(wing_span*0.7) # [m]
chord_tip = taper_ratio * chord_root # [m]
chord_mean_aerodynamic = (2/3)*chord_root * ((1+taper_ratio + \
taper_ratio**2)/(1+taper_ratio)) # [m]
mean_aerodynamic_y = (wing_span/6) * ((1+2*taper_ratio)/(1+taper_ratio)) # [m] # taken from middle 
mean_aerodynamic_x = mean_aerodynamic_y * m.tan(sweep_angle_rad) # [m]
tan_quarter_chord_line_angle_rad = m.tan(sweep_angle_rad) + (chord_root/(2*wing_span)) * (taper_ratio-1)
quarter_chord_line_angle_rad = m.atan(tan_quarter_chord_line_angle_rad)
aspect_ratio_limit = 17.7 * (2-taper_ratio) * m.exp(-0.043 * quarter_chord_line_angle_rad)













