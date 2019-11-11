function [temperature_at_flight_height, rho_at_flight_height] = airdensity_calculator(flight_altitude)
% Summary of this function goes here
%   Detailed explanation goes here

% Determine the temperature lapse rate from the flight altitude
if (flight_altitude <11000) 
    temperature_lapse_rate = -0.0065;
    T_0 = 288.15;
    h_0 = 0;
elseif (flight_altitude >= 11000 && flight_altitude < 20000)
    temperature_lapse_rate = 0;
    T_0 = 216.65;
    h_0 = 11000;
elseif (flight_altitude >= 20000 && flight_altitude < 32000)
    temperature_lapse_rate = 0.001;
    T_0 = 216.65;
    h_0 = 20000;
elseif (flight_altitude >= 32000 && flight_altitude < 47000)
    temperature_lapse_rate = 0.0028;
    T_0 = 228.65;
    h_0 = 32000;
elseif (flight_altitude >= 47000 && flight_altitude < 51000)
    temperature_lapse_rate = 0;
    T_0 = 270.65;
    h_0 = 47000;
elseif (flight_altitude >= 51000 && flight_altitude < 71000)
    temperature_lapse_rate = -0.0028;
    T_0 = 270.65;
    h_0 = 51000;
elseif (flight_altitude >= 71000 && flight_altitude < 71000)
    temperature_lapse_rate = -0.002;
    T_0 = 214.65;
    h_0 = 71000;
end

% to track the height from the relevant height level
dh = flight_height - h_0;
temperature_at_flight_height = T_0 + temperature_lapse_rate * dh;
rho_at_flight_height = rho_0;
end

