% Comppute aerodynamic proprteries of a trapezoidal wing with a geometric
% twist

% describe wing specifications
AR=6;
Lambda=0.55; 
rho=1.226;
twist=-4; % degrees
V=250; % km/hr
WLoad=800;
A=zeros(1,8);

% define initial conditions
k=8; % spanwise station

for j=1:k % j=1 is wing tip
    theta(j)=pi*j/(2*k); 
    costh(j)=cos(theta(j));
    sinth(j)=sin(theta(j));
    Y(j)=costh(j); % nondimensionalized spanwise distance relative to root chord
    C(j)=1-(1-Lambda)*costh(j); % ratio of wing station to root chord
end
    
% Compute coefficients D(j,N) in Eqn (6.40a)
for j=1:k
    D1=1/C(j);
    D2=pi/(AR*(1+Lambda)*sinth(j));
    for N=1:k
        I=2*N-1;
        D(j,N)=(D1+D2*I)*sin(I*theta(j));
    end
end

% Choose two values of Al1 and Al2 for absolute angle of attack at the root
% Calculate the absolute angle of attach at all wing stations, solve Eqns
% (6.40a) using cramers rule and then compute lift coefficient and
% sectional lift coefficients using Eqn (6.33) and (6.41a)

AL1=3; % degrees
AL2=6; % degrees

% Calculations for AL1
for j=1:k
    Alabs(j)=(AL1+twist*costh(j))*pi/180;
end
A=Cramer(D,Alabs,A,k);
CLW1=pi^2*A(1)/(1+Lambda);
for j=1:k
    sum=0;
    for N=1:k
        sum=sum+A(N)*sin((2*N-1)*theta(j));
        CL1(j)=2*pi/C(j)*sum;
    end
end

% Calculations for AL2
A=zeros(1,8);
for j=1:k
    Alabs(j)=(AL2+twist*costh(j))*pi/180;
end
A=Cramer(D,Alabs,A,k);
CLW2=pi^2*A(1)/(1+Lambda);
for j=1:k
    sum=0;
    for N=1:k
        sum=sum+A(N)*sin((2*N-1)*theta(j));
        CL2(j)=2*pi/C(j)*sum;
    end
end

% solve Eqns (6.50) for CLA and CLB
for j=1:k
    CLA(j)=(CL2(j)-CL1(j))/(CLW2-CLW1);
    CLB(j)=CL1(j)-CLA(j)*CLW1;
end

% compute the Wing Loading
CLWF=WLoad/(0.5*rho*(V*1000/3600)^2);

% assuming a linear relation between clw and the absolute angle of attack
% at the root, that angles under flights condition is 
ALF=AL1+(AL2-AL1)*(CLWF-CLW1)/(CLW2-CLW1);

% calcualte section lift coeff using either Eqn (6.46) or (6.41a)
% with Eqn (6.46)
for j=1:k
    CL_1(j)=CLB(j)+CLA(j)*CLWF;
end
% with Eqn (6.41a)
A=zeros(1,8);
for j=1:k
    Alabs(j)=(ALF+twist*costh(j))*pi/180;
end
A=Cramer(D,Alabs,A,k);
for j=1:k
    sum=0;
    for N=1:k
        sum=sum+A(N)*sin((2*N-1)*theta(j));
        CL_2(j)=2*pi/C(j)*sum;
    end
end
% We now compute sectional induced angles of attack using eq. (6..42a)
% sectional induced drag coefficients using eq (6.51), and sectional lift
% curve slope using eqs (6.48) and then proint results.
for j=1:k
    sum=0;
    for N=1:k
        I=2*N-1;
        sum=sum+I*A(N)*sin(I*theta(j))/sinth(j);
        Alind(j)=-pi/(AR*(1+Lambda))*sum;
        CDind(j)=-CL_1(j)*Alind(j); % change this line for the right CL
        M(j)=CL_1(j)/Alabs(j); % change this line for the righ CL
        Alind(j)=Alind(j)*180/pi;
    end
end
% compute induced drag coefficient, CDindW, weighted mean slope, mbar, and
% the additional absolute angle of attack, AlabsSW, of the wing using
% Eqn(6.35a), (6.53a), and (6.54) respectively and print results
sum=0;
for N=1:k
    sum=sum+(2*N-1)*A(N)^2;
    CDindW=pi^3/(AR*(1+Lambda)^2)*sum;
end
sum=0;
KM1=k-1;
for j=1:KM1
    sum=sum+(M(j)*C(j)+M(j+1)*C(j+1))*(Y(j)-Y(j+1));
    Mbar=sum/(1+Lambda);
    AlabsSW=CLWF/Mbar*180/pi;
end
