function [Sh,Ch_tip,Ch_root,bh,ih,ARh,lambdah,Lambdah,Gammah,CLh,Cma,CLt]= HorizontalTailSizing(W,Df,Vcruise,Mac,AR,lambda,iw,S,Lambda, Gamma, CLaw,at, Vh,aoaW)
% Computes all necessary tail size parameters, reference example 6.2 on pg
% 331 in Mohammad H Sadraey's textbook
%
% USAGE:[Sh,Ch_tip,Ch_root,bh,ih,ARh,lambdah,Lambdah,Gammah]= 
%       HorizontalTailSizing(M,D,Vc,Mac,AR,lambda,iw,S,alpha_twist,Lambda, Gamma, Claw)
%
% INPUT:
%
% W             (1,1)              Weight at takeoff
% Df            (1,1)              Largest aft fusealge diameter
% Vcruise       (1,1)              Cruise Velocity
% Mac           (1,1)              Mean aerodynamic chord
% AR            (1,1)              Aspect Ratio of Wing
% lambda        (1,1)              Taper Ratio
% iw            (1,1)              Incidence angle of wing
% S             (1,1)              Wing area
% Lambda        (1,1)              Sweep Angle
% Gamma         (1,1)              Dihedral angle of the wing
% CLaw          (1,1)              Coeficient of lift slope vs alpha for
%                                  the wing
% at            (1,1)              Twist of the wing
% Vh            (1,1)              Horizontal Tail Volume Coefficient:
%                                  Choose from Table 6.4 (pg 303) in
%                                  Mohommad H Sadraey)
% aoaW          (1,1)              Angle of attack of the wing
%
% OUTPUT: 
%
% Sh            (1,1)              Platform area of horizontal tail
% Ch_tip        (1,1)              Horizontal tail tip chord
% Ch_root       (1,1)              Horizontal tail root chord
% bh            (1,1)              Horizontal tail span
% ih            (1,1)              Incidence angle of horizontal tail
% ARh           (1,1)              Aspect Ratio of horizontal tail
% lambdah       (1,1)              Horizontal tail taper ratio
% Lambdah       (1,1)              Sweep angle of horizontal tail
% Gammah        (1,1)              Dihedral angle of horizontal tail
% CLh           (1,1)              Horizontal tail lift coeficcinet at
%                                  cruise
% Cma           (1,1)              Moment coefficent at cruise angle of
%                                  attack for horizontal tail
% CLt           (1,1)              Created lift coefficeint of the tail

%% DEFINE CONSTANT VARIABLES
h0 = .25;
ath = 0.000001; % twist angle of horizontal tail

% Calcualte Tail Platform area
lopt = 1.2*sqrt(4*Mac*S*Vh/(Df*pi)); % Eqn 6.47: Optimal tail arm to minimize aircraft drag [m]
Sh = Mac*S*Vh/lopt; % Eqn 6.24: Tail platform area [m^2]

CL = 2*W/(1.225*Vcruise^2*S); % Eqn 6.27: Cruise Lift Coefficient
Cmaf = -0.025; % Airfoil Sectional pitching moment coefficient From table 5.2 Sadraey\

Cmowf = Cmaf*((AR*cosd(Lambda)^2)/(AR+2*cosd(Lambda)))+0.01*at; % Eqn 6.26: Pitching moment coefficient wings/fuselage

Lf = lopt/0.6; % Fuselage length [m] (reference Table 5.2 of Sadraey)

Xcg = 0.25*Mac-0.114; % Location of center of mass from leading edge in terms of MAC [m]
h = Xcg/Mac; % Percent of MAC

CLh = (Cmowf + CL*(h-h0))/Vh; % Cruise Tail lift coefficient

ARh = 2/3*AR; % Eqn 6.59: Initial horizontal tail aspect ratio
lambdah = lambda; % initially same as wing 
Lambdah = Lambda; % Same as wing Sweep angle [deg]
Gammah=Gamma; % same as wing dihedral [deg]
Clah = .1111111; %Find from book - double check this

CLah = Clah/(1+Clah/(pi*ARh)); % Eqn 6.57: tail lift curve slope [1/rad]
ah = CLh/CLah; % Eqn 6.51: Tail aoa at cruise [deg]

% To calculate the tail created lift coefficient, the lifting line theory 
% is employed. The following calculates the tail lift coefficient with an 
% angle of attack of ah deg.

N = 9; % number of segments
bh = sqrt(ARh*Sh); % tail span [m^2]
MACh = Sh/bh; % mean aerodynamic chord
Ch_root = (1.5*(1+lambdah)*MACh)/(1+lambdah+lambdah^2); % root chord
theta = pi/(2*N):pi/(2*N):pi/2;
alpha = ah+ath:-ath/(N-1):ah; % segmnets angle of attack

c = Ch_root*(1-(1-lambdah)*cos(theta)); % mean aerodynmaic chord at each segment
mu = c*CLah/(4*bh);
LHS = mu.*(alpha/57.3); % left hand side

% Solving N equations to find coefficients A(i)
for i=1:N
    for j=1:N
        B(i,j)= sin((2*j-1)*theta(i))*(1+(mu(i)*(2*j-1))/sin(theta(i)));
    end
end
A = B\transpose(LHS);
for i=1:N
    sum1(i) = 0;
    sum2(i) = 0;
    for j=1:N
        sum1(i) = sum1(i)+(2*j-1)*A(j)*sin((2*j-1)*theta(i));
        sum2(i) = sum2(i)+A(j)*sin((2*j-1)*theta(i));
    end
end
CLt = pi*AR*A(1); % created lift coefficeint of the tail

epi0 = 2*CL/(pi*AR); % Eqn 6.55 for downwash
depida = 2*CLaw/(pi*AR); % Eqn 6.56
epi = epi0+depida*aoaW; % Eqn 6.54
ih = ah-1+epi; % Eqn 6.53: Incidence angle/tail setting angle for horizontal tail

% calculate the rest of the horisontal tail paramenters
MACh=2/3*Ch_root; % Eqn 6.65: mean aerodynamic chord of horizontal tail
ARh=bh/MACh; % Eqn: 6.63: Aspect ratio of horizontal tail
Cma = CLaw*(h-.25)-CLah*.98*Sh/S*(lopt/C-h)*(1-depida); % Eqn 6.67
Ch_tip=gammah*Ch_root; % Eqn 6.64: Chord at the tip of the horizontal tail


