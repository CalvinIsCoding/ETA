function [CL1,y_s,CL_TO] = liftingLineTheory(N,S,AR,lambda,alpha_twist,i_w, CLalpha,alpha0,b,MAC,Croot,theta,alpha)
%% Lifting Line Theory Code
% This code calculates the coefficient of lift over an entire wing using
% the Prandtl thin airfoil theory.
%
% Source: Aircraft Design: A Systems Engineering Approach 
% Author: Mohammad H. Sadraey
% Edited by: Michael Smith
%
% N = 9; % (number of segments - 1)
% S = 11.6; % m^2
% AR = 10; % Aspect ratio
% lambda = 0.55; % Taper ratio
% alpha_twist = -4; % Twist angle (deg)
% i_w = 5; % wing setting angle (deg) (also called angle of incidence)
% CLalpha = 5.4113; % lift-curve slope for NACA 644-421 (1/rad)
% alpha_0 = -2.5; % zero-lift angle of attack (deg)
% b = sqrt(AR*S); % wing span (m)
% MAC = S/b; % Mean Aerodynamic Chord (m)
% Croot = (1.5*(1+lambda)*MAC)/(1+lambda+lambda^2); % root chord (m)
% theta = pi/(2*N):pi/(2*N):pi/2;
% alpha = i_w+alpha_twist:-alpha_twist/(N-1):i_w;

%% Neglecting Flaps
clc
clear




% segment’s angle of attack
z = (b/2)*cos(theta);
c = Croot * (1 - (1-lambda)*cos(theta)); % Mean Aerodynamics
% Chord at each segment (m)
mu = c * CLalpha / (4 * b);
LHS = mu .* (alpha-alpha_0)/57.3; % Left Hand Side


% Solving N equations to find coefficients A(i):
for i=1:N
    for j=1:N
        B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) * (2*j-1)) / sin(theta(i)));
    end
end

A=B\transpose(LHS);

for i = 1:N
    sum1(i) = 0;
    sum2(i) = 0;
    for j = 1 : N
        sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
        sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
    end
end

CL = 4*b*sum2 ./ c;
CL1=[0 CL(1) CL(2) CL(3) CL(4) CL(5) CL(6) CL(7) CL(8) CL(9)];
y_s=[b/2 z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9)];


% Plot results
figure(1)
plot(y_s,CL1,'-o')
grid
title('Lift distribution with no flaps');
subtitle('using NACA 64_4-421 airfoil');
xlabel('Semi-span location (m)');
ylabel ('Lift coefficient');
CL_wing = pi * AR * A(1);

%% With Flaps

clc
clear

N = 9; % (number of segments-1)
S = 11.6; % m^2
AR = 10; % Aspect ratio
lambda = 0.55; % Taper ratio
alpha_twist = -4; % Twist angle (deg)
i_w = 5; % wing setting angle (deg)
a_2d = 5.4113; % lift curve slope (1/rad) (for NACA 64_4-421)
cf_c = 0.2; % flap-to-wing-chord ratio
bf_b = 0.6; % flap-to-wing span ratio
delta_f = 26; % flap deflection (deg)

a_0 = -1.15*(cf_c)*delta_f; % flap up zero-lift angle of attack (deg)
a_0_fd = a_0*2; % flap down zero-lift angle of attack (deg) 
b = sqrt(AR*S); % wing span
MAC = S/b; % Mean Aerodynamic Chord
Croot = (1.5*(1+lambda)*MAC)/(1+lambda+lambda^2); % root chord

theta = pi/(2*N) : pi/(2*N) : pi/2;
alpha = i_w+alpha_twist : -alpha_twist/(N-1) : i_w;

% segment’s angle of attack
for i=1:N
    if (i/N)>(1-bf_b)
        alpha_0(i)=a_0_fd; %flap down zero lift AOA
    else
        alpha_0(i)=a_0; %flap up zero lift AOA
    end
end

z = (b/2)*cos(theta);
c = Croot * (1 - (1-lambda)*cos(theta)); % MAC at each segment
mu = c * a_2d / (4 * b);
LHS = mu .* (alpha-alpha_0)/57.3; % Left Hand Side

% Solving N equations to find coefficients A(i):
for i=1:N
    for j=1:N
        B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) *(2*j-1)) / sin(theta(i)));
    end
end

A=B\transpose(LHS);

for i = 1:N
    sum1(i) = 0;
    sum2(i) = 0;
    for j = 1 : N
        sum1(i) = sum1(i) + (2*j-1) * A(j) *sin((2*j-1)*theta(i));
        sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
    end
end

CL = 4*b*sum2 ./ c;
CL1=[0 CL(1) CL(2) CL(3) CL(4) CL(5) CL(6) CL(7) CL(8) CL(9)];
y_s=[b/2 z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9)];

% Plot results
figure(2)
plot(y_s,CL1,'-o')
grid
title('Lift distribution with \delta_f = 26 degrees');
subtitle('using NACA 64_4-421 airfoil');
xlabel('Semi-span location (m)');
ylabel ('Lift coefficient');

CL_TO = pi * AR * A(1);

