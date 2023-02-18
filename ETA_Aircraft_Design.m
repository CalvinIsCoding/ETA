%% ETA Aircraft Design
% 
%   Course: AEM 4333 - Senior Design
%   Term: Spring 2023
%   Revision #: 1 
%   Last Modified: February 18, 2023
% 
%   Group Name: Ecological Training Aircraft (ETA) 
%   Members: Michael Smith, Erne Habegger-McCabe, Ben Schley, Calvin
%   Karalus, Mitchell Heim, Clare Graney-Dolan, Ebenezer Pierre
%
%   Usage: The purpose of this script is to perform all necessary
%          calculations in order to obtain all geometric quantities for our
%          general aviation aircraft design.

%*************************************************************************%
% Define Known Constants for our Airplane

AR = 10;  % Wing aspect ratio (dimensionless) 
Lambda = 0.55; % Taper ratio of the wing (dimensionless)
Rho = 1.226 * (0.062428); % Atmospheric density at altitude (lbs/ft^3)
twist = -4; % Geometric twist of the wing (in degrees)
V = 277.8 * (0.911344); % Aircraft velocity (ft/s)
WLoad = 1197.01 * (0.020885472); % Wing loading (lbs/ft^2)

N = 9; % (number of wing segments - 1)
S = 11.6 * (10.7639); % ft^2
i_w = 5; % wing setting angle (in degrees) (also called angle of incidence)
CLalpha = 5.4113; % lift-curve slope for NACA 644-421 (1/rad)
alpha_0 = -2.5; % zero-lift angle of attack (deg)
b = sqrt(AR*S); % wing span (ft)
MAC = S/b; % Mean Aerodynamic Chord (ft)
Croot = (1.5*(1+Lambda)*MAC)/(1+Lambda+Lambda^2); % root chord (ft)
theta = pi/(2*N):pi/(2*N):pi/2;
alpha = i_w+twist:-twist/(N-1):i_w;
flaps = 'false'; % no flaps

% Mass of hygrogen tanks
M_b=52.895; % [lb] Mass of large hydrogen tank
V_b=9.239; % [ft^3]Volume of large hydrogen tank (NOTE: taken from sheets,
% not same as one given in discord)
rho=M_b/V_b; % [lb/ft^3] density of material assuming uniform denisty
V_s=pi*(9.84/12)^2*(26.26/12); % [ft^3] v=pi*r^2*h:volume of small tank
M_s=rho*V_s; % [lb] mass of small tank assuming uniform density

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
%*************************************************************************%

% Call applicable functions to calculate geometric parameters
% Wing Design 
[CLWF, CDindW] = FiniteWing(AR, Lambda, Rho, twist, V, WLoad)

% Tail Design %% NEED MORE INFO!!!
% Don't know some of these variables

% [Sh,Ch_tip,Ch_root,bh,ih,ARh,lambdah,Lambdah,Gammah,CLh,Cma,CLt]= ...
%     HorizontalTailSizing(W,Df,Vcruise,Mac,AR,lambda,iw,S,Lambda, Gamma,CLaw,at,Vh,aoaW)

% Center of Mass
[Outputs] = Center_of_Mass...
    (M_b, V_b, rho, V_s, M_s, M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11)

% Flap Sizing (Lifting Line Theory)
[CL1,y_s,CL_wing] = liftingLineTheory...
    (N,S,AR,Lambda,twist,i_w,CLalpha,alpha_0,b,MAC,Croot,theta,alpha,flaps)
