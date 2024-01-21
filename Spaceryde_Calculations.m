%%%%%%%%%%%%%%%%%%%%%%%
% SpaceRyde Design Optimization Challenge
% 
% Name : Sepehr Salimi
%
% Acknowledgements: 
% N/A
%
% Description:
% MATLAB code to use numerical optimization to maximize the delta-v of the
% rocket subject to various constraints and parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; % reset
tic % timer

%% -------------------------------------------------
% Exhaust velocity

% Values
g = 9.80665; % m/s^2. Assumed gravity value is at sea level throughout flight as a conservative value.
Isp = 250; % Seconds. Specific impulse is time for which the thrust equals the weight of propellant consumed, when thrust and flow rate are constant.
k_air = 1.3928; R = 287; % Values for air at sea level to calculate speed of sound. First is specific heat, second is specific gas constant in J/kg/k.
temp_sea = 273.15 + 25; % assuming room temperature (25 deg celsuis)

% Calculations
global ve; % makes exhaust velocity into global, to be fed into Optimize.m
ve = Isp*g; % exhaust velocity in m/s
v_sound_sea = sqrt(k_air*R*temp_sea); % speed of sound at sea level

%% -------------------------------------------------
% Combustion product adiabatic flame temperature (assuming constant volume)

% Values
N_kerosene = 1; % mole values
N_o2 = 12.5;
N_co2 = 12;
N_h2o = 13;

temp_ad_guess = 2300; % based on kerosene adiabatic temperature with air [8]
t_evaluated = 0.5 * (temp_sea + temp_ad_guess); % will be using 1300 Kelvin for Cp values to obtain an estimate for adiabatic temperature. Already implemented.

Cp_kerosene_298 = 238.141; % kJ/(kmol * K). from table [6]
Cp_o2_298 = 29.315; % rest are from Appendix A and B tables [8]
Cp_co2_1300 = 56.987; 
Cp_h2o_1300 = 45.027;

h_kerosene_298 = -276.813 * 1000; % kJ/kmol. % from table [6]
h_o2_298 = 0; % rest are from Appendix A and B tables [8]
h_co2_298 = -393546; 
h_h2o_298 = -24845;

Ru = 8.315; % J/(mol * K). ideal gas constant

% Calculations

syms temp_ad

Hr_tot = N_kerosene * h_kerosene_298 + N_o2 * h_o2_298; % total enthalpy of reactants
Hp_tot = N_co2 * (h_co2_298 + Cp_co2_1300 * (temp_ad - temp_sea) ) + N_h2o * (h_h2o_298 + Cp_h2o_1300 * (temp_ad - temp_sea)); % total enthalpy of products

N_react = N_kerosene + N_o2; % sum of N reactants
N_prod = N_co2 + N_h2o; % sum of N products

v_deltaP = Ru * (N_react * temp_sea - N_prod * temp_ad); % delta internal energy = 0, or delta enthalpy = change in pressure * volume

temp_ad = vpasolve(v_deltaP == Hr_tot + Hp_tot, temp_ad); % adiabatic flame temperature of products in Kelvin
temp_ad = double(temp_ad); % makes into double instead of syms

%% ----------------------------------------
% Combustion product pressure

% Values
P = 2.5 * 10^6; % Pa. Pressure per reactant
Pr = 2 * P; % 2 tanks, so total reactant pressure is 2*P

% Calculations
Pc = Pr * N_prod * temp_ad / (N_react * temp_sea); % pressure of products after combustion, using ideal gas law

%% ------------------------------------------------
% Exhaust temperature and mach #

% Values
k = 1.28; % Assumed value of specific heat of combustion (CO2's value) for simplicity [10]

% Calculations
syms temp_exhaust v_sound_exhaust M_exhaust

E1 = temp_exhaust == temp_ad * ( 1 + M_exhaust^2 * (k-1)/2 )^(-1);  % exhaust temperature though converging diverging nozzle
E2 = v_sound_exhaust == sqrt(k * R * temp_exhaust); % speed of sound at the end of the nozzle
E3 = M_exhaust == ve / v_sound_exhaust; % Mach # of engine

result = solve(E1,E2,E3); % solving system of equations

temp_exhaust = vpa(result.temp_exhaust); temp_exhaust = double (temp_exhaust); % making exhaust temperature into double format
M_exhaust = vpa(result.M_exhaust); M_exhaust = double (M_exhaust); % making exhaust mach # into double format
v_sound_exhaust = vpa(result.v_sound_exhaust); v_sound_exhaust = double (v_sound_exhaust); % making exhaust speed of sound into double format


%% ------------------------------------------------
% Combustion and exhaust density

% Values
M_CO2 = 44.01 ; % g/mol. molar mass of CO2 product found through summation of elements in periodic table
M_H2O = 18.01528 ; % g/mol. molar mass of H2O product found through summation of elements in periodic table

% Calculations
Mc_avg = ( (12/25) * M_CO2 + (13/25) * M_H2O ) / 1000 ; % Average density of gas mixture in kg/mol

global rho_combustion; % make rho_combustion global to be fed into Optimize.m
rho_combustion = (Pc * Mc_avg) / (Ru * temp_ad); % assumed ideal gas and using average molar mass of product mixture

rho_exhaust = rho_combustion * (1 + M_exhaust^2 * (k-1) / 2) ^ (k-1); % exit density of products, using density and speed at end of converging-diverging nozzle

%% -------------------------------------------------
% Exhaust Pressure

% Variables
Patm = 101325; % atm pressure at sea level

% Calculations
Pe = Pc * (1 / (1 + M_exhaust^2 * (k-1)/2) ) ^ (k/(k-1)); % exhaust pressure, using combustion pressure and exhaust speed at end of converging-diverging nozzle. 

if Pe >= Patm
    fprintf('Since Pe > Patm, no shockwave occurs within the engine!\n Proceeding with Optimization! \n\n')
else 
    msg = "RESULTS INVALID, SHOCKWAVE OCCURS SINCE Pe < Patm!";
    f = msgbox(msg)
    error(msg)
end

%% ------------------------------------------
% Exhaust area and diameter

% Values
T = 10000; % Assume constant thrust = 10000 N

% Calculations

Ar = 1 / M_exhaust * ((1 + M_exhaust^2 * (k-1)/2) / (1+(k-1)/2) )^ ((k+1) /(2*k-2)); % Area ratio of Aexhaust / Athroat

syms At; % Letting At become a sysmatic variable to be solved in the next line
At = vpasolve( T/(Pc*At) == Ar * (Pe/Pc - Patm/Pc) + k * sqrt( 2/(k-1) * (2/(k+1)) ^ ( (k+1)/(k-1) ) * ( 1 - (Pe/Pc) ^ ( (k-1)/k ) ) ), At); % Area of throat calculated through convering diverging nozzle
At = double (At); % Making At to double format

global Ae; % making Ae global value to be fed into Optimize.m
Ae = At * Ar; % exhaust area
de = sqrt(4*Ae/pi); % exhaust diameter

%% -----------------------------------------------------------
% Running Optimize.m

ObjFcn = @OptimizeParameters;
nvars = 4; % Number of variables
LB = [0.5 0 0 0]; % Lower Bound, 0.5m is to prevent the solution being in centimeters for diamter and height as it's not practical.
UB = []; % Upper Bound

A=[4 -1 0 0;0 0 3 3]; % Inequalities
b=[0;10000];

rng default; % resets random generator

options = optimoptions('ga', 'PopulationSize', 10, 'Generations', 25,...
    'MaxGenerations',35,'MaxStallGenerations', 100,'TolFun',0.5e-2,...
    'PlotFcn', @gaplotbestf,'Display','iter');

[x, fval, exitFlag, Output] = ga(ObjFcn, nvars, A, b, [], [], LB, UB, [], options);


% Printing

fprintf('The number of generations was : %d\n', Output.generations); % # of generations
fprintf('The number of function evaluations was : %d\n', Output.funccount); % # of functions evaluated

%% ----------------------------------------------------------
% Time to burn all propelent, maximum velocity, max height and its time

% First, the rocket burns propellent and reaches a height in which it no
% longer has any fuel left. This is its ballistic stage. The solution will
% be broken into two phases.

% Calculations

m_propellent_rate = rho_combustion * Ae * ve; % kg/s. rate which propellent is used up
t2 = x(3)/m_propellent_rate; % seconds. The time in which propellent is all used up and rocket enters ballistic stage.
v_max = ve * log(1 + (x(3)/x(4)) ) - g*t2; % tells maximum velocity

m_propellent_rate = rho_combustion * Ae * ve; % kg/s. rate which propellent is used up

t2 = x(3)/m_propellent_rate; % seconds. The time in which propellent is all used up and rocket enters ballistic stage.

height2 = v_max * t2 / 2; % height at which propellent is all burned, using average velocity
height3 = height2 + v_max^2 / (2*g); % max height

syms t3 % solving for t3

t3 = vpasolve(height3 - height2 == v_max * t3 - 0.5 * g * t3^2, t3); 
t3 = double(t3); % making syms into double
t3 = t3 + t2; % total time till max height

% Print

fprintf('The maximum velocity is : %g m/s\n\n', v_max); % maximum velocity
fprintf('The diameter of the rocket was found to be %g meters, and the total height to be %g meters \n',x(1),x(2));
fprintf('The mass of the propellents initially will weigh %g kg, and the mass of the rocket with empty tanks weighs %g kg \n\n', x(3),x(4));
fprintf('The amount of time that the fuel is burned for is : %g seconds\n', t2); % maximum velocity
fprintf('The maximum height it reaches is : %g meters\n', height3); % maximum height
fprintf('And the time it takes for it to reach this is : %g seconds\n\n', t3); % maximum height's time

save Results

toc