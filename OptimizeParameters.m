function y = OptimizeParameters(x)
% x(1) is diameter of rocket
% x(2) is total height of rocket
% x(3) is the mass of the propellent before any is burned
% x(4) is the mass of the struct of the rocket (with tanks empty)

% Values
rho_oxid = 1141; % kg/m^3. O2 density
rho_fuel = 810; % kg/m^3. C12H26, Kerosene density
rho_al = 2810; % kg/m^3. 7075 Aluminum Alloy density, assuming entire rocket is made of this.
g = 9.80665; % m/s^2. Assumed gravity value is at sea level throughout flight as a conservative value
m_engine = 15; % kg. Assumed only 1 engine of 15 kg
m_payload = 10; % kg.
T = 10000; % N. Assumed constant thrust
thickness = 0.002; % m. The thickness of the body tube and fairing. 
m_ratio = 2.5; % 2.5 O2 : 1 kerosene ratio.

% Global Values (RUN Spaceryde_Calculations FIRST!)
global rho_combustion; % kg/m^3. Exhaust combustion. Calculated in Spaceryde_Calculations.m 
global ve; % m/s. Exhaust speed. Calculated in Spaceryde_Calculations.m 
global Ae; % m^2. Exhaust area. Calculated in Spaceryde_Calculations.


% Calculations
% global h_fuel; h_fuel = 4 * x(1) * rho_oxid / (15 * rho_fuel); % the fuel is a cylinder
% global h_oxid; h_oxid = x(1); % since O2 tank is a sphere, the height is the diameter

h_fuel = x(1); % since kerosene tank is a sphere, the height is the diameter
h_oxid = 2.5 * x(1) * rho_fuel / (3 * rho_oxid); % the oxid tank is a cylinder

m_fairing = rho_al * thickness * pi * x(1)/2 * (x(1)/2 + sqrt((x(2)-h_oxid-h_fuel)^2 + x(1)^2/4)); % density * thickness * cone shape surface area [5]
m_body = pi * x(1) * (h_oxid + h_fuel) * thickness * rho_al; % Assuming height of body = h_oxid and h_fuel, since most is used for fuel.
m_fueltankempty = rho_al * thickness * pi * x(1) * h_fuel; % mass of cylindrical fuel tank, empty
m_oxidtankempty = rho_al * thickness * pi * x(1)^2; % mass of O2 spherical tank, empty

m_oxid = pi/6 * x(1)^3 * rho_oxid; % volume * density for spherical shape
m_fuel = pi/4 * x(1)^2 * h_fuel * rho_fuel; % surface area * height * density for cylindrical shape

x(3) = m_oxid + m_fuel; % mass of propellent before any burned
x(4) = m_fairing + m_payload + m_engine + m_body + m_fueltankempty + m_oxidtankempty; % mass of rocket structure (empty tanks)

m_propellent_rate = rho_combustion * Ae * ve; % kg/s. rate which propellent is used up
t2 = x(3)/m_propellent_rate; % seconds. The time in which propellent is all used up and rocket enters ballistic stage.

% Objective function
y = -1 * (ve * log(1 + (x(3)/x(4)) ) - g*t2); % maximize objective function, delta-v. Assume no drag.

save OptmizeResults
end

