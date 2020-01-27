import math as m
drum_diameter = 0.53 #m
drum_radius = 0.5*drum_diameter 
angle_force_leverarm = 90 #degrees
drum_length = 3.45 # meters
density_rod = 1500 #kg/m3 https://www.wiley.com/college/callister/1118061608/case_studies/ch01.pdf
#--------------------------------- Reel out -----------------------------------------
power = 635000 #[W]
reel_out_speed = 3.44 # [m/s]
tether_force = power/reel_out_speed # [N]

cable_length = 800 # [m]
cable_weight = 300#kg
reel_out_time = cable_length/reel_out_speed # seconds
reel_out_time_min = reel_out_time /60
energy_generated = power * reel_out_time # [J]
energy_generated_kwh = energy_generated/(1000*3600)
rot_speed = reel_out_speed/drum_radius # [[rad/s]
rot_speed_rpm = rot_speed * 60/(2*m.pi)
torque_out = drum_radius * tether_force * m.sin(angle_force_leverarm*m.pi/180)
#--------------------------------- Reel in ----------------------------------------

reel_in_speed_mmin = 25*60# m/min
reel_in_speed = reel_in_speed_mmin/60
reel_in_length = cable_length #m
reel_in_time_min = reel_in_length/reel_in_speed_mmin #min
reel_in_time = reel_in_time_min * 60 #seconds
reel_in_rot_speed = reel_in_speed/drum_radius # rad/s

mass_rod = 0.5 * cable_weight + m.pi*(drum_radius**2)* drum_length * density_rod
I = 0.5 * mass_rod *drum_radius**2
winch_energy_needed = 0.5*I*reel_in_rot_speed**2 # [J]
winch_energy_needed_kwh = winch_energy_needed/(1000*3600) # [kWh]
power_needed = winch_energy_needed/reel_in_time

#--------------------how many systems do we need
net_energy = energy_generated - winch_energy_needed
net_energy_kwh = net_energy / (1000*3600)
average_net_power = net_energy/(reel_out_time+reel_in_time)
tu_needed_energy_kwh = 87600000
amount_of_cycles_needed_year = tu_needed_energy_kwh/net_energy_kwh
number_of_cycles_year = (60*24*365)/(reel_in_time_min+reel_out_time_min)
number_of_systems_needed = amount_of_cycles_needed_year / number_of_cycles_year
tu_needed_power = 10000000
amount_needed = tu_needed_power/average_net_power 
generator_efficiency = 0.97
amount_needed_after_generator = amount_needed/generator_efficiency
winch_efficiency = net_energy/energy_generated

area_runway_onesystem = m.pi*(300+10)**2 - m.pi*300**2
area_onesystem = area_runway_onesystem + 6.93 + 4.14
tdfyguh = (m.pi*310**2)*19
# power_density = (average_net_power*19) / (area_onesystem + (6.93 + 4.14)*18)
power_density = (average_net_power*19) / (tdfyguh)

rot_frequency = 0.0166666666666667*rot_speed_rpm

kwh_day = ((average_net_power * 3600* 24) / (1000*3600)) * amount_needed 
kwh_year = kwh_day * 365
cost_generator = kwh_year/13.5326