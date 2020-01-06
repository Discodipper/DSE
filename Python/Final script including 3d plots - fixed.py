# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 10:27:25 2019

@author: wesley & pranav
"""


from math import *
import matplotlib as mpl
import numpy as np
import time
from mpl_toolkits import mplot3d

#parameters to include a legend into the plot
mpl.rcParams['legend.fontsize'] = 10

#end of legend to plot


start = time.time()

#emptyweight = np.arange(250, 500, 10)#300 #kg
altitude = np.arange(2000, 2600, 100) #10000 #m
surfacearea = np.arange(15, 41, 1) #m^2
reelspeed = np.arange(10, 21, 1)
WingLoading = np.arange(400, 610, 10)
OEW_lst = []
alt_lst = []
CL_lst = []
S_lst = []
V_reel_lst = []
P_lst = []
wingloading_lst = []
V_a_lst = []
gliders = []
thrust_lst = []
VTO_lst = []
dneeded_lst = []
tension_lst = []
spacing_lst = []
W_cable_lst = []
P_density_WRTground_lst = []
for alt in altitude:
    for S in surfacearea:
        for V_reel in reelspeed:
            for wingloading in WingLoading:
                L = wingloading * S #N
                OEW = (S + 0.0013)/0.0375 #kg
                operating_angle = 22*pi/180 #rad
                max_angle = 30*pi/180
                
                
                """ISA calculations"""
                a = -0.0065
                T_0 = 288.15
                g = 9.80665
                p_0 = 101325
                rho_0 = 1.225
                
                
                T = T_0 + (a*alt)
                p = (p_0) * (T/T_0)**((-g)/(a*287))
                rho = (rho_0) * (T/T_0)**(((-g)/(a*287))-1)            
                
                
                
                V_w = 7*10**-8*alt**2 + 0.0011*alt + 7.5388 #m/s #between 1-10km altitude, R^2 = 0.997
                """ https://kitepower.nl/resources/ch14_awe15_fechner.pdf"""
                V_k = V_reel * np.cos(operating_angle) #m/s
                V_a = V_w + V_k #m/s
                
                d_cable = 0.008 #m
                cable_density = 970 #kg/m^3
                L_cable = alt/(np.sin(operating_angle)) #m
                W_cable = cable_density*L_cable*pi*(d_cable/2)**2 * g #N
                W_glider = OEW*g #N
                T_y = L - W_cable - W_glider #N
                tension = T_y/np.sin(operating_angle) #N
                if tension > 0:
                    SF = 3
                    d_needed = 2*((tension*SF)/(3.9*(10**9)*np.pi))**0.5 #m
                    gen_eff = 0.8
                    generator_power = gen_eff*V_reel*tension/1000 #kW
                    glider_number =  10*10**3/(generator_power)
                    
                    CL = L/(0.5*rho*V_a**2*S)
                    
                    

                    angle_mean = (max_angle + operating_angle)/2
                    angle_diff = angle_mean - operating_angle
                    d_u = L_cable/(sin(angle_mean - angle_diff) * (1/tan(angle_mean - angle_diff)) * (1/(tan(2*angle_diff))))
                    P_density_WRTground = generator_power/(d_u**2)
                    
                    if CL <= 2.2 and d_needed <= d_cable and generator_power >= 300:
                        OEW_lst.append(OEW)
                        alt_lst.append(alt)
                        CL_lst.append(CL)
                        V_reel_lst.append(V_reel)
                        S_lst.append(S)
                        P_lst.append(generator_power)
                        wingloading_lst.append(wingloading)
                        V_a_lst.append(V_a)
                        gliders.append(glider_number)
                        dneeded_lst.append(d_needed)
                        tension_lst.append(tension)
                        
                        V_TO = np.sqrt((2*OEW*g)/(rho_0 * S * (CL+0.3))) #m/s
                        climb_angle = 10*pi/180 #rad
                        CD = 0.07
                        drag = 0.5*CD*rho_0*S*V_TO**2
                        thrust = drag + OEW*g*sin(climb_angle)
                        
                        thrust_lst.append(thrust)
                        VTO_lst.append(V_TO)
                        
                        spacing_lst.append(d_u)
                        P_density_WRTground_lst.append(P_density_WRTground)
                        W_cable_lst.append(W_cable)
design_number =  263  #insert here                       
print("OEW = ", OEW_lst[design_number], "kg")
print("Operational altitude = ", alt_lst[design_number], "m")
print("CL needed = ", CL_lst[design_number])
print("Reel speed = ", V_reel_lst[design_number], "m/s")
print("Wing surface area = ", S_lst[design_number], "m^2")
print("Power of single unit = ", P_lst[design_number], "kW")
print("Number of gliders = ", gliders[design_number])
print("Thrust needed for takeoff = ", thrust_lst[design_number], "N")
print("Spacing between gliders = ", spacing_lst[design_number], "m")
print("Power density with respect to ground = ", P_density_WRTground_lst[design_number]*1000, "W/m^2")                    
                    

fig1 = mpl.pyplot.figure(1)
ax = mpl.pyplot.axes(projection="3d")
scatter1 = ax.scatter3D(OEW_lst, wingloading_lst, P_lst, c=P_lst, cmap="inferno")
ax.set_title("Structural plot design 2")
ax.set_xlabel("OEW [kg]")
ax.set_ylabel("Wing loading [$N/m^2$]")
ax.set_zlabel("Power [kW]")
legend1 = ax.legend(*scatter1.legend_elements(),
                    loc="lower left", title="Classes")
ax.add_artist(legend1)

fig2 = mpl.pyplot.figure(2)
ax = mpl.pyplot.axes(projection="3d")
scatter2 = ax.scatter3D(alt_lst, CL_lst, P_lst, c=P_lst, cmap="inferno")
ax.set_title("Operational plot design 2")
ax.set_xlabel("Altitude [m]")
ax.set_ylabel("$C_L$ [-]")
ax.set_zlabel("Power [kW]")
legend2 = ax.legend(*scatter2.legend_elements(),
                    loc="lower left", title="Classes")
ax.add_artist(legend2)
             
fig3 = mpl.pyplot.figure(3)
ax = mpl.pyplot.axes(projection="3d")
scatter3 = ax.scatter3D(OEW_lst, S_lst, P_lst, c=P_lst, cmap="inferno")
ax.set_title("Sizing plot design 2")
ax.set_xlabel("OEW [kg]")
ax.set_ylabel("Surface area [$m^2$]")
ax.set_zlabel("Power [kW]")
legend2 = ax.legend(*scatter3.legend_elements(),
                    loc="lower left", title="Classes")
ax.add_artist(legend2)
mpl.pyplot.show()




end = time.time()
print("Elapsed time: ", end-start)