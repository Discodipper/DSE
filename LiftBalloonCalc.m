% ISA CALCULATIONS
g = 9.80665; %gravity constant m/s^2
h = [];
i = 1;
F_b = [];
F_bc = [];
%Outp_year = [];
hc = [];
Num_solp = [];
% SOLAR PANEL CALCULATIONS
Num_balloon = 50;
Area_solp = 1.63; % [m^2]
t_solp = 0.048; %[m]
W_solp = 18.6*g; % [N]
rho_solp = (W_solp/g)/(Area_solp*t_solp); %[kg]
T_b = 298.15; % [K]
Outp_0 = 370; % Output van 1 solar panel [J/s]  
eff_0 = 0.227; %  [frac]
tot_outp_0 = Outp_0/eff_0; % [W]
temp_coef = -0.0029; % [frac/deg C] 
Req_eng = 70000; % [MWh] TU Delft yearly energy requirement
for alt = 2000:500:20000
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
    
    Temp_diff = T - T_b;
    frac_diff = Temp_diff * temp_coef;
    eff = eff_0 + frac_diff;
    Outp(i) = tot_outp_0 * eff * 5; % Total output * efficiency * 5 for being above clouds LOL
    Outp_year(i) = (Outp(i)*3600*24*365)/(3600*1000000); %Output of 1 solar panel per year in [MWh]
    Num_solp(i) = ceil(Req_eng/Outp_year(i)); % Minimum number of solar panels
    Tot_W_solp(i) = Num_solp(i) * W_solp;
    %MASS ESTIMATIONS
    % Super-pressure balloon
    
    % Zero-pressure balloon
    
    % Balloon calculations:
    rho_i = 0.0837765634; % density inside balloon
    Vol(i) = Tot_W_solp(i) / (g * (rho - rho_i)); %Volume entire system [m^3]
    Tot_area_solp(i) = (Tot_W_solp(i) / rho_solp) / t_solp; %Total area solp of system [m^2]
    
    %Horizontal spacing calculations
    
    Vol_balloon(i) = Vol(i) / Num_balloon; %Volume per balloon [m^3]
    Area_solp(i) = Tot_area_solp(i) / Num_balloon; %area of solp per balloon
    Spacing(i)  = (Area_solp(i) + 200*200) * Num_balloon; %Spacing needed for entire system 
 
    i = i+1;
    


end




%plot3(h,V,F_b)
%contour3(hc,Vc,F_bc,25)
