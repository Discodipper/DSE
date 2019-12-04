from math import *
# ISA CALCULATIONS
g = 9.80665
h = []
F_b = []
F_bc = []
#Outp_year = []
hc = []
Num_solp = []
# SOLAR PANEL CALCULATIONS
Num_balloon = 100
Area_solp = 0.0027 #[m^2]
t_solp = 0.008 #[m]
W_solp = 0.135*g # [N]
rho_solp = (W_solp/g)/(Area_solp*t_solp)
T_b = 298.15 # [K]
Outp_0 = 1.155 # Output van 1 solar panel [J/s]  
eff_0 = 0.316 #  [frac]
tot_outp_0 = Outp_0/eff_0 # [W]
temp_coef = -0.004 # [frac/deg C] 
Req_eng = 70000 #[MWh] TU Delft yearly energy requirement

alt = range(2000,32500,500)
rho_lst = [];
for i in alt:
    #[T, rho, p] = airdensity_calculator(alt);
    
    if 0 <= i <=11000:
        p_0 = 101325 #Pa
        T_0 = 288.15 #K
        h_0 = 0. #@m
        R_air = 287.
        a_lapse = -0.0065
        T = T_0 + a_lapse * (i-h_0) #Temperature for gradient layer
        p = p_0 * (T/T_0)**(-g/(a_lapse*R_air))
        rho = p/(R_air*T)
        
        
    elif 11000 < i <=20000:
        p_0 = 22700. #Pa
        h_0 = 11000.
        T = 216.8 #K
        rho_0 = 0.364805
        p = p_0 * exp((-g/(R_air*T))*(i-h_0))
        rho = rho_0 * (p/p_0)
        
        
        
    #elif 20000 <= i <=32000:
    #    p_0 = 5474.9
    #    h_0 = 20000.
     #   T_0 = 216.8
      #  R = 287.
       # a_lapse = -0.001
        #T = T_0 + a_lapse * (i-h_0) #Temperature for gradient layer
         #p = p_0 * (T/T_0)**(-g/(a_lapse*R_air))
          #rho = p/(R_air*T)
    
    rho_lst.append(rho)
    
    
perc = []

for i in range(len(rho_lst)): 
    perc.append(((rho_lst[i-1] - rho_lst[i]) / (rho_lst[0] - rho_lst[len(rho_lst)-1])) * 100)
perc = perc[1:]
print(perc)
o_lst = []

#h = int(input("What's your altitude in m per 0/500/1000/1500/etc.: "))
for h in range(2000,20500,500):
    d = perc[alt.index(h)-1]
    x = sum(perc[:(perc.index(d))+1]) 
    o = 0.23*x

    o_lst.append(o)
print(o_lst[1:])
#i_lst = []
#for i in o_lst[1:]:
#    i_lst.append(i)
#print(i_lst)



