% NormInvGaussianScript simulates a log asset price
% and thus the asset price path by a repeated call to
% the function NormInvGaussianGen, which is subordinated
% by the InvGaussianGen time process

close all; clear all; clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

global dX

dt=0.01
T=0:dt:1;

delta = 1 
alpha = 1
beta = 0
mu=0

if (-alpha>beta)
    disp('error: -alpha > beta');
elseif (beta > alpha)
    disp('error: beta > alpha');
end
    
dtvector=ones(1,length(T))*dt;
for j = 1:length(T)
    [ dX(j), dIG(j) ] = NormInvGaussianGen(alpha, beta, delta, dt );
end

IG=cumsum(dIG);
X=cumsum(dX);

figure
subplot (4,1,1)
plot (T, IG)
title ('Inverse-Gaussian Process')
ylabel ('IG')
xlabel ('Time')

subplot (4,1,2)
hist(dX,length(T)/4)
ylabel ('Frequency (dX)')
xlabel ('dX')

subplot (4,1,3)
plot (T, X)
title ('Normal Inverse-Gaussian log (Price) ')
ylabel ('X')
xlabel ('Time')
axis tight

eX=exp(X);
subplot (4,1,4)
plot (T,eX)
ylabel ('e^X')
title ('Normal Inverse-Gaussian Asset Price ')
xlabel ('Time')
axis tight

figure
stem3(T,IG,X)
hold on
plot3(T,IG,X,'gr')
axis tight
xlabel ('Real Time')
ylabel ('Inverse-Gaussian Time')
zlabel ('X')
title ('Normal Inverse Gaussian Process')

E = mean(dX)
%myE = ( (1/length(dX))*sum((dX)) )
V = var(dX)
%myV = ( (1/length(dX))*sum((dX-E).^2) )
S = ( (1/length(dX))*sum((dX-E).^3) ) / V^1.5
K = ( (1/length(dX))*sum((dX-E).^4) ) / V^2

% use defaults when Method of Moments calculates
% nonsensical numbers 
if (( (3*K-4*S^2) < 9 ) |...
        ((3*K-3*S^2-9) / ( S*(V-(5/4)*S^2-3)^2) <0) ) 
    disp('MoM Error') ;
    alphaMoM = 1;
    betaMoM = 0;
    deltaMoM = 0.01;
    muMoM = 0;
else
    alphaMoM = sqrt( (3*K-3*S^2-9) / ( S*(V-(5/4)*S^2-3)^2))
    betaMoM = S/ ( sqrt(V) * (K-(5/3)*S^2-3))
    deltaMoM = sqrt( (27*V*(K-(5/3)*S^2-3))...
        / (3*K-4*S^2-9) )
    muMoM = E - ( (3*S*sqrt(V)) / (3*K-4*S^2-9) )
end

alphaEst = alphaMoM;
betaEst = betaMoM; 
deltaEst = deltaMoM; 
muEst = muMoM; 

% MIGlikelihood function can accept three (mu assumed zero)
% or four parameters
param=[ alphaEst, betaEst, deltaEst, muEst];
%param=[ alphaEst, betaEst, deltaEst]

%Opt=optimset('Display','iter');
Opt=[];
[pnew,likelihood]=fminsearch('NIGlikelihood',param,Opt); 
likelihood=-likelihood;
fprintf(1, 'Likelihood = %6.4f \n', likelihood);

calcAlpha=pnew(1)
calcBeta=pnew(2)
calcDelta=pnew(3)
calcMu=pnew(4)







