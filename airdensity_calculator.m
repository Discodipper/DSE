function [temperature_at_flight_altitude, rho_at_flight_altitude, pressure_at_flight_altitude] = airdensity_calculator(flight_altitude)
% Summary of this function goes here
%   Detailed explanation goes here

R_air = 8.3144598; %N·m/(mol·K)
g_0 = 9.80665; %m/s2
molar_mass_of_air_earth = 0.0289644; %kg/mol

% Determine the temperature lapse rate from the flight altitude
if (flight_altitude <11000) 
    temperature_lapse_rate = -0.0065;
    T_0 = 288.15;
    h_0 = 0;
    p_0 = 101325;
    rho_0 = 1.2250;
elseif (flight_altitude >= 11000 && flight_altitude < 20000)
    temperature_lapse_rate = 0;
    T_0 = 216.65;
    h_0 = 11000;
    p_0 = 22632.10;
    rho_0 = 0.36391;
elseif (flight_altitude >= 20000 && flight_altitude < 32000)
    temperature_lapse_rate = 0.001;
    T_0 = 216.65;
    h_0 = 20000;
    p_0 = 5474.89;
    rho_0 = 0.08803;
elseif (flight_altitude >= 32000 && flight_altitude < 47000)
    temperature_lapse_rate = 0.0028;
    T_0 = 228.65;
    h_0 = 32000;
    p_0 = 868.02;
    rho_0 = 0.01322;
elseif (flight_altitude >= 47000 && flight_altitude < 51000)
    temperature_lapse_rate = 0;
    T_0 = 270.65;
    h_0 = 47000;
    p_0 = 110.91;
    rho_0 = 0.00143;
elseif (flight_altitude >= 51000 && flight_altitude < 71000)
    temperature_lapse_rate = -0.0028;
    T_0 = 270.65;
    h_0 = 51000;
    p_0 = 66.94;
    rho_0 = 0.00086;
elseif (flight_altitude >= 71000 && flight_altitude < 71000)
    temperature_lapse_rate = -0.002;
    T_0 = 214.65;
    h_0 = 71000;
    p_0 = 3.96;
    rho_0 = 0.000064;
end

% to track the height from the relevant height level
dh = flight_altitude - h_0;
temperature_at_flight_altitude = T_0 + temperature_lapse_rate * dh;
pressure_at_flight_altitude = p_0 * (T_0 / temperature_at_flight_altitude)^(molar_mass_of_air_earth* g_0 / (R_air * temperature_lapse_rate)); %exp((-g_0/(R_air*temperature_at_flight_altitude))*dh);
rho_at_flight_altitude = rho_0 * (pressure_at_flight_altitude/p_0);

end

