function [Sh,Ch_tip,Ch_root,bh,ih,ARh,lambdah,Lambdah,Gammah]= HorizontalTailSizing(M,D,Vc,Mac,AR,lambda,iw,S,alpha_twist,Lambda, Gamma, Claw)
% Computes all necessary tail size parameters 
%
% USAGE:[Sh,Ch_tip,Ch_root,bh,ih,ARh,lambdah,Lambdah,Gammah]= 
%       HorizontalTailSizing(M,D,Vc,Mac,AR,lambda,iw,S,alpha_twist,Lambda, Gamma, Claw)
%
% INPUT:
% M             (1,1)              Mass at takeoff
% D             (1,1)              Maximum fusealge diameter
% Vc            (1,1)              Cruise Velocity
% Mac           (1,1)              Mean aerodynamic chord
% AR            (1,1)              Aspect Ratio
% lambda        (1,1)              Taper Ratio
% iw            (1,1)              Incidence angle of wing
% S             (1,1)              Wing area
% alpha_twist   (1,1)              Angle of twist of the wing
% Lambda        (1,1)              
% Gamma         (1,1)              Dihedral angle of the wing
% Claw          (1,1)              Coeficient of lift slope vs alpha for
%                                  the wing
% OUTPUT: 
% Sh            (1,1)              Platform area of horizontal tail
% Ch_tip        (1,1)              Horizontal tail tip chord
% Ch_root       (1,1)              Horizontal tail root chord
% bh            (1,1)              Horizontal tail span
% ih            (1,1)              Incidence angle of horizontal tail
% ARh           (1,1)              Aspect Ratio of horizontal tail
% lambdah       (1,1)              Horizontal tail taper ratio
% Lambdah       (1,1)              Sweep angle of horizontal tail
% Gammah        (1,1)              Dihedral angle of horizontal tail
%% DEFINE CONSTANT VARIABLES
Vh = 0.7; %horizontal tail Volume coefficient
C = 1.077; %Wing MAC (m)
S = 11.6; %Wing Area (m^2)
Df = 1.524; %Largest Aft fuselage diameter (m)
CLaw = 5.4113;
W = 12946.1; %Weight
Vcruise = 46.3; %m/s

lopt = 1.4*sqrt(4*C*S*Vh/(Df*pi)); %Optimal tail arm in meters
Sh = C*S*Vh/lopt; %taill area m^2

CL = 2*W/(1.225*Vcruise^2*S); %Cruise Lift Coefficient
AR = 10; %Wing AR
Cmaf = -0.025; %Airfoil Sectional pitching moment coefficient From table 5.2 Sadraey
Lamda = 4.5; %sweep angle
at = -1e-6; %Twist of the wing
Cmowf = Cmaf*((AR*cosd(Lamda)^2)/(AR+2*cosd(Lamda)))+0.01*at; %Pitching moment coefficient of the wings and uselage but not really the second one

Lf = lopt/0.6;
Xapex = -0.25*C+0.32*Lf+.144; %I don't know what this is
Xcg = 0.25*C-0.114; %meters from leading edge
h = Xcg/C; % %MAC
h0 = .25;
CLh = (Cmowf + CL*(h-h0))/Vh; %Cruise Tail lift coefficient

ARh = 5.59;
taperh = .55; %initially same as wing
Lamdah = 5; %Same as wing Sweep
Gamma = 0.00001; %Same as wing Dihedral
Clah = .1111111; %Find from book - double check this

CLah = Clah/(1+Clah/(pi*ARh)); %tail lift curve slope
%ah = CLh/CLah; %Tail aoa at cruise
ah = -3.89;

%Last part
N = 9;
bh = sqrt(ARh*Sh);
Ch = Sh/bh;
Cr = (1.5*(1+taperh)*Ch)/(1+taperh+taperh^2);
theta = pi/(2*N):pi/(2*N):pi/2;
ath = 0.000001;
alpha = ah+ath:-ath/(N-1):ah;

z = (bh/2)*cos(theta);
c = Cr*(1-(1-Lamdah)*cos(theta));
mu = c*CLah/(4*bh);
LHS = mu.*(alpha/57.3);

for i=1:N
    for j=1:N
        B(i,j)= sin((2*j-1)*theta(i)*(1+(mu(i)*(2*j-1))/sin(theta(i))));
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
CLt = pi*AR*A(1);

epi0 = 2*CL/(pi*AR);
depida = 2*CLaw/(pi*AR);
aw = 2; %Check for cruise aoa
epi = epi0+depida*aw;
ih = ah-1+epi;

Cma = CLaw*(h-.25)-CLah*.98*Sh/S*(lopt/C-h)*(1-depida);

%% Vertical Tail

Vv = 0.04;
b = 10.7703;
Sv = b*S*Vv/lopt;
ARv = 1.84;
taperv = .55;
iv = 0;
lamdaV = 18;

bv = sqrt(Sv*ARv);
Cv = bv/ARv;
Cvr = (3/2)*Cv*((1+taperv)/(1+taperv+taperv^2));
Cvt = taperv*Cvr;

%% Elevator Sizing

%Rotation time during take off: 1-3 seconds
%Take-off pitch angular acceleration: 8-10 deg/s/s

be = .9*bh; %Elevator span (m)
deltapmax = 20; %max positive deflection (down)
deltanmax = -25; %Max nedative deflection (up)

Vstall = 28.2944; %Stall speed in m/s - 55 knots
Vr = 30.8666; %Rotation Speed 60 knots
angR = 9; %Take off pitch angular acceleration (deg/s^2)

qr = .5*1.225*Vr^2;
CLTO = 1.1589; %Coefficient of lift for take off TODO get exact value
Lwf = qr*CLTO*S; %Take off Lift
Lhto = qr*Sh*CLh; %Take off tail lift
Lto = Lwf+Lhto;

CDTO = .05; %Coefficient of Drag at take off TODO Get number
Dto = qr*CDTO*S;

Cmacwf = -0.04; %Coefficient of moment about ac for the wings number
Macwf = qr*Cmacwf*S*C;

T = 125/Vr; %Thrust at take off get this number
m = W/9.81; %Mass
a = (T-Dto)/m; %Linear Acceleration at take off

Xmg = 3.149; %Distance of main landing gear from nose cone (m)
Xcgn = 2.8;
Mw = W*(Xmg-Xcgn); %Weight moment

Zd = .889; %Drag Height I guess
Zmg = 0; %just the wheels
Zcg = .889; %Center of gravity height
Md = Dto*(Zd-Zmg); %Drag moment

Zt = 1.016; %thrust height
Mt = T*(Zt-Zmg); %Thrust moment

Xacwf = .889; %wing ac
Mlwf = Lwf*(Xmg-Xacwf); %Wing lift moment

Ma = m*a*(Zcg-Zmg); %Acceleration moment

Xach = lopt;
Iyy = (1/12*282.91*C^2+(3/5*23.308*(.25*(Df/2)^2+1.69^2)...
+23.308*Xmg^2)+(3/5*89.533*(.25*(Df/2)^2+4.386^2)+...
89.533*(Lf-Xmg)^2)+1/12*60.516*(3*(Df/2)^2+2.163^2))+m;
Lh = (Mlwf+Macwf+Ma+Mw+Md+Mt+Iyy*angR)/(Xach-Xmg);

CLhto = 2*Lh/(qr*Sh);

%te = ((ah/57.3)+(CLhto/CLah))/(deltanmax/57.3);
te = .644;
Ce = Ch/2;

%% Rudder Sizing

Vapr = 1.1*Vstall; %Approach landing velocity
Vw = 8.746; %Max Cross wind m/s - 17 kn
Vtot= sqrt(Vapr^2+Vw^2);

Sf = Lf*Df; %Side fuselage area
Ss = 1.02*(Sf+Sv); %Side surface area

Xca = (Sf*Lf/2+Sv*(Lf-Cv/2))/(Sf+Sv);

dc = Xca-Xcgn;

Fw = .5*1.225*(Vw^2)*Ss*.6;

sideslip = atan(Vw/Vapr);
%Cnb = .75*

