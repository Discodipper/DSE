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
Num_balloon = 100;
A_sp = 0.0027; % [m^2]
t_solp = 0.008; %[m]
W_solp = 0.135*g; % [N]
rho_solp = (W_solp/g)/(A_sp*t_solp); %[kg]
T_b = 298.15; % [K]
Outp_0 = 1.155; % Output van 1 solar panel [J/s]  
eff_0 = 0.316; %  [frac]
tot_outp_0 = Outp_0/eff_0; % [W]
temp_coef = -0.004; % [frac/deg C] 
Req_eng = 70000; % [MWh] TU Delft yearly energy requirement

%Weights
W_antenna = 5278; %[N]
W_balloon = 0; %[N]
W_cables = 300; %[N]
W_battery = 10000; %[N]


for alt = 2000:500:20000
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
     
    end
    
    Temp_diff = T - T_b;
    frac_diff = Temp_diff * temp_coef;
    eff = eff_0 + frac_diff;
    Outp(i) = tot_outp_0 * eff * 5; % Total output * efficiency * 5 for being above clouds LOL
    Outp_year(i) = (Outp(i)*3600*24*365)/(3600*1000000); %Output of 1 solar panel per year in [MWh]
    Num_solp(i) = ceil(Req_eng/Outp_year(i)); % Minimum number of solar panels
    Tot_W_solp(i) = Num_solp(i) * W_solp;
    
    % Balloon calculations when stationary altitude:
    rho_i = 0.0837765634; % density inside balloon
    Vol(i) = Tot_W_solp(i) / (g * (rho - rho_i)); % Volume entire system when stationary[m^3]
    Tot_area_solp(i) = (Tot_W_solp(i) / rho_solp) / t_solp; %Total area solp of system [m^2]
    
    %Horizontal spacing calculations
    
    Vol_balloon(i) = Vol(i) / Num_balloon; %Volume per balloon [m^3]
    Area_solp(i) = Tot_area_solp(i) / Num_balloon; %area of solp per balloon
    %Spacing(i)  = (Area_solp(i) + 200*200) * Num_balloon; %Spacing needed for entire system 
    
    
    
    %D = C_D * 0.5 * V^2 * Area_balloon * rho
    %W(i) = (Tot_W_solp(i)/Num_balloon) + W_balloon + W_cables + W_antenna + W_battery; %Total weight of one balloon [N]
    %L(i) = Vol(i)*g*(rho - rho_i) - W(i); %- D; 
    %a(i) = (L(i)*g)/W(i);
    
    i = i+1;
    
    

end
nexttile
plot(2000:500:20000,Vol_balloon)
nexttile
plot(2000:500:20000,Area_solp)

%plot3(h,V,F_b)
%contour3(hc,Vc,F_bc,25)

