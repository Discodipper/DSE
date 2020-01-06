function [cableweight] = ElectricalCableCalculations[t, L, rho]

%Input of diameter, resistivity, length
%I = input(prompt) %A
%t = 1 %mm
%rho = 1 %kg/m^3
%resistivity = 1 %Ohm*m
%L = 1 %m


cableweight = pi*((t/1000)/2)^2*L*rho

end

function[cablethickness] = ElectricalCableCalculations[eff, L]

t = 