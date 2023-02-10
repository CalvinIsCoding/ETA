unction [Outputs] = Center_of_Mass(M_b, V_b, rho, V_s, M_s, M1, M2, M3, M5, M6, M7, M8, M9, M10, M11)
% Last File Edit: 2 - 9 - 23
% Last Editor: Mitchell Heim 
%
% This purpose of this function is so calculate location of the center of
% mass on the aircraft under various weight conditions. 
%
% 2 weight conditions exist:
%   1) Empty Weight 
%           Weight of 02 and H2 = Zero
%
%   2) Full Weight 
%           Weight of 02 and H2 = Full Mass
%
% The result of this function is a plot of how the location of the COM
% changes based on percentage of fuel weight remaining 
%
% Inputs: 
%   M_b: Weight of large H2 tank [lb]
%   V_b: Volume of large H2 tank [ft^3]
%   rho: Density of fuel tank material [lb/ft^3]
%   V_s: Volume of smaller 02 tank [ft^3]
%   M_s: Weight of smaller 02 tank [lb]
%   M1: Weight of Fuel Cell Stack 
%   M2: Weight of O2 Tank
%   M3: Weight of Battery 
%   M4: NOT AN INPUT 
%   M5: Weight of Motor Assembly 
%   M6: Weight of Propeller 
%   M7: Weight of Nose Section
%   M9: Weight of Cabin Section 
%   M9: Weight of Tail Section 
%   M10: Weight of Wings 
%   M11: Weight of Tail

%
% Outputs 
%   x_CM_empty: Location of COM with empty fuel weight [ft]
%   x_CM_full: Location of COM with full fuel weight [ft]
%   x_CM_interm: Vector locatio of COM based on various fuel weight [ft]
%



%% Code Section 

% -------------------------------------------------------------------------
% Code used for setting up the loop for moving COM plots

totalSteps = 1001;
stepPercentage = linspace(1, 0, totalSteps); % 1x1001 vector from 100% to 0%
% 100% is full fuel weights, 0% is empty fuel weight

 %% This section should be placed in the master file except for in code calculations 
% -------------------------------------------------------------------------
% Mass of hygrogen tanks
M_b=52.895; % [lb] Mass of large hydrogen tank
V_b=9.239; % [ft^3]Volume of large hydrogen tank (NOTE: taken from sheets,
% not same as one given in discord)
rho=M_b/V_b; % [lb/ft^3] density of material assuming uniform denisty
V_s=pi*(9.84/12)^2*(26.26/12); % [ft^3] v=pi*r^2*h:volume of small tank
M_s=rho*V_s; % [lb] mass of small tank assuming uniform density

% -------------------------------------------------------------------------
% Masses of all components that do not change *FIXED COMPONENT WEIGHT* 
M1=800.052; % [lb] Fuel Cell Stack
M2=22.04; % [lb] Weight of O2 tank 
M3=116.812; % [lb] Weight of Battery
M4=2*M_s; % [lb] Weight of both H2 tanks beacuse they are in the same position
M5=198.36; % [lb] Wieght of motor assembly
M6=78; % [lb] Weight of propeller
M7=51.3861; % [lb] Weight of nose section
M8=97.3858; % [lb] Weight of cabin section
M9=133.4139; % [lb] Weight of tail section
M10=623.71; % [lb] Weight of wings
M11=38.5+61; % [lb] Weight of Tail


%% Above is stand alone section to be moved 
% ------------------------------------------------------------------------- 
% Adjustable / Constantly changing weights (Fuel Weight)
M2_2=22.04+11.02; %[lb] Weight of O2 tank and O2
M4_2=2*(M_s+32.619); % [lb] Weight of H2 tanks and H2

% Changing Weights - These are vectors that indicate the fuel weight based
% on percentage of fuel in the aircraft 
M2_2per = (M2_2-M2)*stepPercentage;
M4_2per = (M4_2-M4)*stepPercentage;


% -------------------------------------------------------------------------
% Distances from quarter chord with quarter chord being x=0
x1=8.43/12; % [ft] Distance from fuel cell stack to 1/4 chord
x2=-75.26/12; % [ft] Distance from O2 tank to 1/4 chord
x3=-73.91/12; % [ft] Distance from battery to 1/4 chord
x4=-54.22/12; % [ft] Distance from H2 tanks to 1/4 chord
x5=-71.54/12; % [ft] Distance from motor to 1/4 chord
x6=-110.28/12; % [ft] Distance from nose to 1/4 chord
x7=(-22.17-80.37)/12; %[ft] Distance from nose center of mass to 1/4 chord
x8=(42.58-80.37)/12; % [ft] Distance from cabin center of mass to 1/4 chord
x9=(142.72-80.37)/12; % [ft] Distance from tail center of mass to 1/4 chord
x10=0; % [ft] Distance from wing to 1/4 chord
x11=214.03/12; % [ft] Distance from Tail to 1/4 chord

% -------------------------------------------------------------------------
% Calculate the x location of the center of mass

% NOTE: These values are distances from the quarter chord
%       Negative Value: Towards Back of Aircraft
%       Positive Value: Towards Front of Aircraft 

x_CM_empty = (M1*x1+M2*x2+M3*x3+M4*x4+M5*x5+M6*x6+M7*x7+M8*x8+M9*x9+M10*x10+M11*x11)...
    /(M1+M2+M3+M4+M5+M6+M7+M8+M9+M10+M11) % [ft] Empty Fuel Weight Location
x_CM_full = (M1*x1+M2_2*x2+M3*x3+M4_2*x4+M5*x5+M6*x6+M7*x7+M8*x8+M9*x9+M10*x10+M11*x11)...
    /(M1+M2+M3+M4+M5+M6+M7+M8+M9+M10+M11) % [ft] Full Fuel Weight Location 

for i=1:totalSteps
    
    %x_Cm_interm is a vector of the location of the center of mass based on
    % percentage of fuel in the tanks
    x_CM_interm(i) = (M1*x1+  M2_2per(i)  *x2+M3*x3+  M4_2per(i)  ...
        *x4+M5*x5+M6*x6+M7*x7+M8*x8+M9*x9+M10*x10+M11*x11)...
    /(M1+M2+M3+M4+M5+M6+M7+M8+M9+M10+M11); % [ft] Full Fuel Weight Location

end

% -------------------------------------------------------------------------
%Define the function Outputs 
Outputs = [x_CM_empty, x_CM_full, x_CM_interm];

% Plot the results 

%Define an x vector
xPos = linspace(-2, 2, totalSteps);

figure();
plot(stepPercentage, x_CM_interm)
grid on;
xlabel('Percentage of Fuel Weight')
ylabel('Position of COM in Feet from Cbar')


%% End of Code Section 


end
