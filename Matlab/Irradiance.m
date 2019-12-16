rho_lst = [];
for alt = 0:500:32000
    %[T, rho, p] = airdensity_calculator(alt);
    
    if (0<alt)&& (alt<=11000)
        p_0 = 101325; %Pa
        T_0 = 288.15; %K
        h_0 = 0; %m
        R_air = 287;
        a_lapse = -0.0065;
        T = T_0 + a_lapse * (alt-h_0); %Temperature for gradient layer
        p = p_0 * (T/T_0)^(-g/(a_lapse*R_air));
        rho = p/(R_air*T);
        
        
    elseif (11000<alt) && (alt<=20000)
        p_0 = 22700; %Pa
        h_0 = 11000;
        T = 216.8; %K
        rho_0 = 0.364805;
        p = p_0 * exp((-g/(R_air*T))*(alt-h_0));
        rho = rho_0 * (p/p_0);
        
        
        
    elseif (20000<alt) && (alt<=32000)
        p_0 = 5474.9;
        h_0 = 20000;
        T_0 = 216.8;
        R = 287;
        a_lapse = -0.001;
        T = T_0 + a_lapse * (alt-h_0); %Temperature for gradient layer
        p = p_0 * (T/T_0)^(-g/(a_lapse*R_air));
        rho = p/(R_air*T);
    end
    rho_lst = [rho_lst,rho];
    
    
end

for i in rho_lst
    


    