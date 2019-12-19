function [power_req] = kWh_to_W_converter_in_days(energy, time)
%Function to transfer kWh to W 
%   With the given kWh required for a certain time in days)
power_req= energy / (time * 24);
end

