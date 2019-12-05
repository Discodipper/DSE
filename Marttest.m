delta = 1;
hlist = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000];
Turbinelist = [1 2 3 4];
Airfoillist = [1 2 3 4];
n_turbineslist = [2 4 6 8];
eff_cablelist = [0.5 0.6 0.7 0.8 0.9]; 
g = 9.80665; %unit 
mturblist = [11 22 33 44]; %kg 
Esystemlist = [11 22 33 44] ;
rholist = [11 22 33 44];
Vlist = [11 22 33 44];
Cllist = [11 22 33 44];
ClCdlist = [11 22 33 44];
Cdturb = 1 ;
rhocable = 8000; %kg/m^3
Ilist = [11 22 33 44];
rho_el = 1.72*10^-8; %Ohm*m
Ulist = [11 22 33 44];
p0 = 101325;
t0 = 288.15;
h0 = 0;
R_air = 287;
A_lapse = -0.0065;

h = 2000 ;
n_turbines = 1 ;
eff_cable = 0.8 ;
Edemand = 10000000 ; %kW
Esysteminit = 5000 ; %kW
mturb = 15 ; %kg
Vinit = 10 ; %m/s
Cl = 1 ;
ClCd = 25 ;
rhocable = 8960 ; %kg/m^3
I = 100 ; %A
U = 1200 ; %V
Cd = Cl/ClCd;
A = 5; %m/s

%S = pils(h, n_turbines, eff_cable, Edemand, Esysteminit, g, mturb, rho, Cl, Cd, ClCd, Cdturb, rhocable, rho_el, U, A, p0, t0, h0, R_air, A_lapse);
%disp(S);
start = 1;
hlist2 = [];
Slist = [];
for h = hlist
    hlist2(start) = h;
    S = pils(h, n_turbines, eff_cable, Edemand, Esysteminit, g, mturb, rho, Cl, Cd, ClCd, Cdturb, rhocable, rho_el, U, A, p0, t0, h0, R_air, A_lapse, Vinit);
    Slist(start) = S;
    start = start+1;
end

%plot(hlist2, S)

function [S] = pils(h, n_turbines, eff_cable, Edemand, Esysteminit, g, mturb, rho, Cl, Cd, ClCd, Cdturb, rhocable, rho_el, U, A, p0, t0, h0, R_air, A_lapse, Vinit)
    
    
    T = t0+A_lapse*(h-h0);
    P = p0*(T/t0)^(-g/(A_lapse*R_air));
    rho = P/(R_air*T); 
    V = (0.002*h+5.0582);
    Esystem = n_turbines*Esysteminit*rho/1.225*(V^3)/(Vinit^3);
    %disp(Esystem)
    I = Esystem/U;
    Lcable = h*sqrt(2);
    t = 2*sqrt(I*rho_el*Lcable/(1-eff_cable)/pi/U);
    mcable = rhocable*Lcable*pi*(t/2)^2;
    Dturb = 0.5*rho*V^2*Cdturb*n_turbines*A;
    S = (Dturb+(mturb*n_turbines+mcable)*g)/(0.5*rho*V^2*(Cl-Cd));
    s_total = S*(Edemand/Esystem/n_turbines);
    dev = s_total/S;
    disp(Esystem)
end
