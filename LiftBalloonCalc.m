% ISA CALCULATIONS
h = [];
i = 1;
F_b = [];
F_bc = [];
V = [];
Vc = [];
hc = [];

for alt = 2000:500:20000
    hc(i) = alt;
    j = 1;
    for Vol = 0:500:500000
        Vc(j) = Vol;
        g = 9.80665; %gravity constant m/s^2
        h = [h, alt];
        V = [V, Vol];
        %h = 8000; %m
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
        
        %MASS ESTIMATIONS
        % Super-pressure balloon
        
        % Zero-pressure balloon
        
        % Balloon calculations:
        rho_i = 0.0837765634; % density inside balloon
        % V = 17000; %m^3
        Force = Vol * g * (rho - rho_i);
        F_b = [F_b,Force];
        F_bc(j,i) = Force;
        
        j = j +1;
        
    end
   i = i +1; 
end



%plot3(h,V,F_b)
contour3(hc,Vc,F_bc,25)
