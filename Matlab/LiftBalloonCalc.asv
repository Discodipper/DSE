% ISA CALCULATIONS
g = 9.80665; %gravity constant m/s^2
radiation = [0,1.243505181776943, 2.439532955387986, 3.5893544197450122, 4.6942214931527735, 5.755366973179216, 6.774004597073443, 7.751329102742852, 8.688516290301822, 9.586723084204786, 10.447087595976798, 11.270729187555142, 12.058748535256372, 12.812227694383278, 13.532230164486982, 14.219800955300181, 14.875966653357649, 15.501735489321305, 16.098097406027478, 16.767416608916115, 17.407330073079443, 17.998751436019802, 18.54535538816441, 19.050538155007143, 19.517438598945684, 19.948957722039346, 20.347776690864386, 20.716373495461028, 21.057038345879562, 21.37188790198919, 21.66287842496402, 21.93181793216053, 22.180377430908745, 22.41010130101615, 22.62241689049432, 22.818643384129427, 23.000000000000004];
i = 1;
Num_solp = [];
Num_balloon = 75;
A_sp = 0.0027; % [m^2] % Area of 1 solar panel
t_solp = 0.008; %[m]
W_solp = 0.135*g; % [N]
rho_solp = (W_solp/g)/(A_sp*t_solp); %[kg/m3]
T_b = 298.15; % [K]
Outp_0 = 305; % Output van 1 solar panel [J/s]  
eff_0 = 0.316; %  [frac]
tot_outp_0 = Outp_0/eff_0; % [W]
temp_coef = -0.004; % [frac/deg C] 
Req_eng = 67908; % [MWh] TU Delft yearly energy requirement
Vol_balloon = [];
Area_solp = [];
Tot_W_solp = [];
Tot_area_solp = [];
Vol = [];
Outp = [];

% ---------------------------- WEIGHTS --------------------------------
W_antenna = 0; %[N]
W_film = 0; %[N]
W_cables = 0*g*Num_balloon; %[N] 
W_battery = 0; %[N]

% ---------------------------- MICROWAVES ----------------------------------------
D = 2000; % [m] seperation between antennas
c = 299792458; % [m/s] speed of microwaves
f = 2.45 * 10^9; % [1/s] frequency of waves
labda = c/f; % [m] wavelength
eff_dc_rf = 0.85; % [-]
eff_atm = 0.722358; % [-]
eff_wireless = eff_dc_rf * eff_atm; % [-]
A_r =5000; % [m^2] 
A_t = [];
% ---------------------------- ALTITUDE-----------------------------

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
    x(i) = ((radiation(i)/100)+1); % Increase of output due to solar radiation 
    Outp(i) = tot_outp_0 * eff * 5 * x(i); % Total output * efficiency * 5 for being above clouds LOL
    Outp_year(i) = (Outp(i)*3600*12*365)/(3600*1000000); %Output of 1 solar panel per year in [MWh]
    Num_solp(i) = ceil(Req_eng/Outp_year(i)); % Minimum number of solar panels
    Tot_Outp(i) = Outp(i) * Num_solp(i); %Power output of entire system 

    Tot_W_solp(i) = Num_solp(i) * W_solp; 
    Tot_W(i) = Tot_W_solp(i) + W_antenna + W_film + W_cables + W_battery;
    
    
    % Balloon calculations when stationary altitude:
    rho_i = 0.0837765634; % density inside balloon
    Vol(i) = Tot_W(i) / (g * (rho - rho_i)); % Volume entire system when stationary[m^3]
    Tot_area_solp(i) = (Tot_W_solp(i) / rho_solp) / t_solp; %Total area solp of system [m^2]
    
    %Horizontal spacing calculations
    
    Vol_balloon(i) = Vol(i) / Num_balloon; %Volume per balloon [m^3]
    Area_solp(i) = Tot_area_solp(i) / Num_balloon; % area of solp per balloon
    Spacing(i)  = (Area_solp(i) + 500 *500) *0.75* Num_balloon; %Spacing needed for entire system 
    
    % ---------------------------- MICROWAVES ----------------------------------------
    Num_solp_balloon(i) = Num_solp(i)/Num_balloon; 
    P_t_tot(i) = Outp(i) * Num_solp_balloon(i) *eff_dc_rf; % [W] Power transmitted per balloon
    P_r_tot(i) = P_t_tot(i) * eff_atm; % [W] Power recieved of one balloon
    A_t(i) = (P_r_tot(i) * labda^2 * alt^2)/(A_r*P_t_tot(i)); 
    Pow_outp_density_ground(i) = Tot_Outp(i) / Spacing(i); 
    Pow_outp_density_air(i) = (Tot_Outp(i)/Num_balloon) / (Area_solp(i));  

    i = i+1;
end
nexttile
plot(Pow_outp_density_air,2000:500:20000)
nexttile
plot(Area_solp,2000:500:20000)
hold on
plot(A_t,2000:500:20000)
hold on




%plot3(Area_solp,2000:500:20000)

