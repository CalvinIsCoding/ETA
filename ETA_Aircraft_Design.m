%% ETA Aircraft Design
% 
%   Course: AEM 4333 - Senior Design
%   Term: Spring 2023
%   Revision #: 4 
%   Last Modified: February 28, 2023
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
Df = 50 /(12); % Largest aft fuselage diameter (ft)
Sweep = 0; % Sweep angle of wing (deg)
Gamma = 5.5; % Dihedral angle of wing (deg) 
Vh = 0.7; % Horizontal Tail Volume Coefficient:
          % Choose from Table 6.4 (pg 303)
aoaW = 2; % angle of attack of the wing
q = 0.5*Rho*V^2; % Dynamic pressure of the wing ()

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

% Vertical Tail Sizing parameters
Vv = 0.04;
ARv = 1.84;
taperv = .55;
iv = 0;
lamdaV = 18;

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

% Center of Mass
[Outputs] = Center_of_Mass...
    (M_b, V_b, rho, V_s, M_s, M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11);

W = Outputs(3);
AR = linspace(0,15,10);
Lambda = linspace(0,0.99,10);
S = linspace(0,500,10);
i_w = linspace(-5,10,10);

% % Optimized wing parameters
% [best] = optimizeWingDistribustion(W, q, S, AR, Lambda, i_w)
% 
% % best = [S AR i_w lamda twist CL_wing]
% W = best(1);
% q = best(2);
% S = best(3);
% AR = best(4);
% Lambda = best(5);
% i_w = best(6);

% Wing Design 
[CLWF, CDindW] = FiniteWing(AR, Lambda, Rho, twist, V, WLoad);

% Horizontal Tail Design
[Sh,Ch_tip,Ch_root,bh,ih,ARh,lambdah,Lambdah,Gammah,CLh,Cma,CLt]= ...
    HorizontalTailSizing(W,Df,V,MAC,AR,Lambda,i_w,S,Sweep,Gamma,CLalpha,twist,Vh,aoaW);

% Vertical Tail Design
[Sv,bv,Cv,Cvr,Cvt] = VerticalTailSizing(Vv,b,ARv,taperv,iv,lambdaV,MAC,Df,S);

% Flap Sizing (Lifting Line Theory)
[CL1,y_s,CL_wing] = liftingLineTheory...
    (N,S,AR,Lambda,twist,i_w,CLalpha,alpha_0,b,MAC,Croot,theta,alpha,flaps);
