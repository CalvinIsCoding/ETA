function [best] = optimizeWingDistribustion(W, q, S, AR, lamda, i_w)
%Computes wing design parameters to obtain a lift distribusion that is
%close to elliptical, minimizing induced drag
%
%output: [best] = [S AR i_w lamda twist CL_wing]
% S - wing surface area
% AR - aspect ratio of the wing
% i_w - setting angle of the wing root
% lamda - taper ratio of the wing
% twist - geometric twist of the wing
% CL_wing - Coefficient of lift of the wings for the geometry calculated
%
% Returns [-1 -1 -1 -1 -1 -1] if the geometric parameters given did not produce
% enough lift for the weight and cruise conditions given
%
%input
% W - estimate of the aircraft weight
% q - dynamic pressure of the aircraft at cruise conditions
% S - 1xN matrix of possible wing areas
% AR - 1xM matrix of possible aspect ratios
% lamda - 1xP matrix of possible taper ratios for the wing
% i_w - 1xR matrix of possible wing setting angles

best = [-1 -1 -1 -1 -1 -1];
bestDiff = 1000000000;

for a=1:length(S)
    for d=1:length(AR)
        for e=1:length(i_w)
            twist = -i_w(e)-1:0.1:-i_w(e)+0.9;
            for j=1:length(lamda)
                for k=1:length(twist)
                    N = 15;
                    a_2d = 5.902;
                    alpha_0 = -2.8;
                    b = sqrt(AR(d)*S(a));
                    MAC = S(a)/b;
                    Croot = (1.5*(1+lamda(j))*MAC)/(1+lamda(j)+lamda(j)^2);
                    theta = pi/(2*N):pi/(2*N):pi/2;
                    alpha = i_w(e)+twist(k):-twist(k)/(N-1):i_w(e);
                    z = (b/2)*cos(theta);
                    c = Croot*(1-(1-lamda(j))*cos(theta));
                    mu = c*a_2d/(4*b);
                    LHS = mu.*(alpha-alpha_0)/57.3;
                    
                    for i=1:N
                        for l=1:N
                            B(i,l) = sin((2*l-1)*theta(i))*(1+(mu(i)*(2*l-1))/sin(theta(i)));
                        end
                    end
                    
                    A = B\transpose(LHS);
                    
                    for i=1:N
                        sum1(i) = 0;
                        sum2(i) = 0;
                        for l=1:N
                            sum1(i) = sum1(i)+(2*l-1)*A(l)*sin((2*l-1)*theta(i));
                            sum2(i) = sum2(i)+A(l)*sin((2*l-1)*theta(i));
                        end
                    end
                    CL = 4*b*sum2./c;
                    CL1 = [0 CL];
                    y_s = [b/2 z];
                    figure (1), clf
                    plot(y_s, CL1, '-o'); hold on;
                    y = sqrt(CL(N)^2-CL(N)^2*(y_s.^2/((b/2)^2)));
                    plot(y_s,y, '-*');
                    CL_wing = pi*AR(d)*A(1);
                    CLreq = W/(q*S(a));
                    if CLreq <= CL_wing
                        diffScore = sum((y-CL1).^2);
                        if diffScore <= bestDiff
                            if diffScore == bestDiff
                                if S(a) <= best(1)
                                    best = [S(a) AR(d) i_w(e) lamda(j) twist(k) CL_wing];
                                    bestDiff = diffScore;
                                end
                            else
                                best = [S(a) AR(d) i_w(e) lamda(j) twist(k) CL_wing];
                                bestDiff = diffScore;
                            end
                        end
                    end
                end
            end
        end
    end
end
end