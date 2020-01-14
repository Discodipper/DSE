import math as m
import numpy as np
combinations = []
angle_of_attack = 3.0 #degrees
V = 100 #m/s
C_L_slope_tail_deg_lst = [0.026788416666667,0.0444725,0.05625483333333334,0.06455525,0.070693583333,0.075410833333] #AR1,2,3,4,5,6    
C_L_slope_tail_rad_lst = [i*180/m.pi for i in C_L_slope_tail_deg_lst]

rho = 1.006 #Consider designing for density at higher altitudes (less lift is generated)
q = 0.5*rho*V**2
angle_of_incidence = 0
sweep_angle_rad = 0 #at LE
x = 0.25
tail_area = 13 # [m^2, laat in report zien dat tail surface area nauwelijks effect heeft op cl en cd]
from shape_parameters import quarter_chord_line_angle_rad
zero_lift_aoa_deg = 7 #IN NEGATIVE DIRECTION [degrees]
z_distance_tail_to_wing = 1 + m.sin(zero_lift_aoa_deg*m.pi/180)# [meter] 
C_L_slope_deg_lst =  []
V_tail = V #SUBJECT TO CHANGE, DISCUSS WITH GROUP
rho_tail = rho
C_m_alpha_lst= []
C_m_alpha_lst2 = []
#length_ac_wing_to_tail
C_D_zero = [0.056,0.019,0.011]
C_L_alpha = [0.09,0.09867,0.10109090909,0.09125,0.1019167,0.10425,0.0919167,0.1033,0.1067] # [1/DEGREES] used for interpolation S20AR6, S20AR12 , S20AR18, S60AR6, S60AR12, S60AR18, S100AR6, S100AR12, S100AR18
C_L_zero_aoa_lst = [0.652,0.748,0.785,0.662,0.762,0.804,0.666,0.766,0.808]
wing_area = 60
for C_L_slope_tail_deg in C_L_slope_tail_deg_lst:
    #for wing_area in range(20,100,5):
    for aspect_ratio in range(6,19,2):
        wing_span = m.sqrt(wing_area*aspect_ratio) #wing span in m
        wing_chord_average = wing_area/wing_span #WHICH CHORD DO WE NEED IN DERIVATION?
        
#--------------------------- WING CL_ALPHA ----------------------------------                        
        if 21<=wing_area<=59 and 7<=aspect_ratio<=11:
            C_L_zero_aoa = ((aspect_ratio-6)/6)*(0.755-0.657)+0.652
            C_L_slope_deg = (0.755-0.657) * ((aspect_ratio - 6)/6) + 0.657
            #C_L_slope_deg = (C_L_alpha[0] + ((wing_area-20)/40) * (C_L_alpha[3]-C_L_alpha[0])) * ((((((aspect_ratio-6)/6)/(C_L_alpha[1]-C_L_alpha[0])/((aspect_ratio-6)/6)/(C_L_alpha[4]-C_L_alpha[3]))-1)*((aspect_ratio-6)/6))+1)
        if 21<=wing_area<=59 and 13<=aspect_ratio<=17:
            C_L_zero_aoa = ((aspect_ratio-12)/6)*(0.7945-0.755)+0.755
            C_L_slope_deg = (0.7945-0.755) * ((aspect_ratio - 12)/6) + 0.755
            #C_L_slope_deg = C_L_alpha[0] + ((wing_area-20)/40) * (C_L_alpha[3]-C_L_alpha[0]) * ((((((aspect_ratio-12)/6)/(C_L_alpha[2]-C_L_alpha[0])/((aspect_ratio-12)/6)/(C_L_alpha[5]-C_L_alpha[3]))-1)*((aspect_ratio-12)/6))+1)
        if 61<=wing_area<=99 and 7<=aspect_ratio<=11:
            C_L_zero_aoa = ((aspect_ratio-6)/6) * (0.764-0.664) + 0.664
            C_L_slope_deg = (0.764-0.664) * ((aspect_ratio - 6)/6) + 0.664
            #C_L_slope_deg = C_L_alpha[0] + ((wing_area-60)/40) * (C_L_alpha[6]-C_L_alpha[3]) * ((((((aspect_ratio-6)/6)/(C_L_alpha[4]-C_L_alpha[3])/((aspect_ratio-6)/6)/(C_L_alpha[7]-C_L_alpha[6]))-1)*((aspect_ratio-6)/6))+1)         
        if 61<=wing_area<=99 and 13<=aspect_ratio<=17:       
            C_L_zero_aoa = ((aspect_ratio-12)/6)*(0.806-0.764)+0.764
            #C_L_slope_deg = C_L_alpha[3] + ((wing_area-60)/40) * (C_L_alpha[6]-C_L_alpha[3])  * ((((    ((aspect_ratio-12)/6)/(C_L_alpha[5]-C_L_alpha[3]) * ((aspect_ratio-12)/6)/(C_L_alpha[8]-C_L_alpha[6]))-1)*((aspect_ratio-12)/6))+1)
            C_L_slope_deg = (0.806-0.764) * ((aspect_ratio - 12)/6) + 0.764
        if wing_area == 20 and aspect_ratio == 6:
            C_L_zero_aoa = C_L_zero_aoa_lst[0]
            C_L_slope_deg = C_L_alpha[0] 
        if wing_area == 20 and aspect_ratio == 12:
            C_L_zero_aoa = C_L_zero_aoa_lst[1]
            C_L_slope_deg = C_L_alpha[1]
        if wing_area == 20 and aspect_ratio == 18:
            C_L_zero_aoa = C_L_zero_aoa_lst[2]
            C_L_slope_deg = C_L_alpha[2]
        if wing_area == 60 and aspect_ratio == 6:
            C_L_zero_aoa = C_L_zero_aoa_lst[3]
            C_L_slope_deg = C_L_alpha[3]
        if wing_area == 60 and aspect_ratio == 12:
            C_L_zero_aoa = C_L_zero_aoa_lst[4]
            C_L_slope_deg = C_L_alpha[4]
        if wing_area == 60 and aspect_ratio == 18:
            C_L_zero_aoa = C_L_zero_aoa_lst[5]
            C_L_slope_deg = C_L_alpha[5]
        if wing_area == 100 and aspect_ratio == 6:
            C_L_zero_aoa = C_L_zero_aoa_lst[6]
            C_L_slope_deg = C_L_alpha[6]
        if wing_area == 100 and aspect_ratio == 12:
            C_L_zero_aoa = C_L_zero_aoa_lst[7]
            C_L_slope_deg = C_L_alpha[7]
        if wing_area == 100 and aspect_ratio == 18:
            C_L_zero_aoa = C_L_zero_aoa_lst[8]
            C_L_slope_deg = C_L_alpha[8]
        C_L_slope_deg_lst.append(C_L_slope_deg)
        #C_L_zero_aoa = 
        C_L = C_L_zero_aoa + C_L_slope_deg*angle_of_attack
        
        C_L_slope_rad = C_L_slope_deg * 180/m.pi
        
        wing_span =m.sqrt( wing_area*aspect_ratio)
        wing_span_half = wing_span/2
        taper_ratio = 0.4
        chord_root = m.sqrt(wing_area/aspect_ratio) * (0.5/(0.35+0.15*x))
        chord_mean_aerodynamic = (2/3)*chord_root * ((1+taper_ratio + \
         import math as m
# #C_D_zero_lst = [0.056,0.012,0.011]
# aspect_ratio = 18
# oswald = 1.78*(1-0.045*aspect_ratio**0.68) - 0.64
# C_D_zero = 0.012
#     C_D = C_D_zero + (C_L**2)/(m.pi * aspect_ratio * oswald)
#                                              taper_ratio**2)/(1+taper_ratio)) # [m]
        quarter_chord_line_angle_rad = m.atan(m.tan(sweep_angle_rad) + (chord_root/(2*wing_span)) * (taper_ratio-1))
        
                
#--------------------------------- TAIL C_L_ALPHA ----------------------------------------  
        for length_acwing_to_actail in range(10,30,2):
            for length_ac_to_cg in np.arange(0.1,5.1,0.1):
                
                m_tv = z_distance_tail_to_wing/ (wing_span_half)

                    
                r = length_acwing_to_actail/(wing_span_half)
                K_epsilon_sweep = (0.1124+0.1265*quarter_chord_line_angle_rad + 0.1766 *quarter_chord_line_angle_rad**2) / (r**2) + (0.1024/r) + 2 
                K_epsilon_zerosweep = (0.1124/(r**2)) + (0.1024/r) + 2 
                downwash_derivative_alpha = (K_epsilon_sweep/K_epsilon_zerosweep) *( (r/(r**2+m_tv**2)) * (0.4876/m.sqrt(r**2 + 0.6319 + m_tv**2)) + (1+(r**2/(r**2 + 0.7915 + 5.0734*m_tv**2))**(0.3113)) * (1-m.sqrt(m_tv**2/(1+m_tv**2)))) *(C_L_slope_rad/(m.pi*aspect_ratio))  # [rad/rad] sead slides
                angle_downwash = (2*C_L/(m.pi * aspect_ratio))*180/m.pi #DEGREES
                angle_of_attack_tail = angle_of_attack - angle_downwash + angle_of_incidence # degrees
                tail_volume = tail_area * length_acwing_to_actail / (wing_area*chord_mean_aerodynamic)
                C_m_alpha = (C_L_slope_deg * length_ac_to_cg / chord_mean_aerodynamic ) - tail_volume * C_L_slope_tail_deg * (1-downwash_derivative_alpha) #WHICH CHORD SHOULD I TAKE?
                L = C_L*q*wing_area + C_L_slope_tail_deg * angle_of_attack_tail *tail_area *0.5*rho_tail * V_tail**2   
                #D = 0.5*C_D_wing*rho*wing_area*V**2 + 
                C_m_alpha_lst.append(C_m_alpha)
                
                
                if C_m_alpha<0 and L>250000:
                    combinations.append([wing_area,aspect_ratio,tail_volume])#,length_acwing_to_actail,length_ac_to_cg, C_L_slope_tail_deg_lst[C_L_slope_tail_deg_lst.index(C_L_slope_tail_deg)]]) #kijken welke mogelijkheden er zijn
                #Thrust_moment = 3*D*x*wing_span_half
                #Rudder_moment = Rudder_force * ac_vertical_tail-
C_m_alpha_lst2.append(C_m_alpha_lst[600000:638000])
                      