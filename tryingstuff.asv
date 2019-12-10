g = 9.80665;
rho_i = 0.0837765634;
i = 1;
j = 1;
for alt = 2000:500:20000
    for Vol = 3000:100:20000
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

        end
        
        L(i) = Vol(i)*g* (rho - rho_i);
    
    
        j = j+1;
    i=i+1;
    
end

