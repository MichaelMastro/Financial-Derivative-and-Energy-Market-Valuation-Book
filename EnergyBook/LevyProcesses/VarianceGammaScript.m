% VarianceGammaScript simulates a log asset price
% and thus the asset price path by a repeated call
% to the function VarianceGammaGen

close all; clear all; clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

dt=0.01;
T=0:dt:1;
sigma=0.1
v=0.1
theta=0.1

dtvector=ones(1,length(T))*dt;
for j = 1:length(T)
    [ dX(j), dG(j) ] = VarianceGammaGen( sigma, v, theta, dt );
end

G=cumsum(dG);
X=cumsum(dX);

figure
subplot (4,1,1)
plot (T, G)
title ('Gamma Process')
ylabel ('G')
xlabel ('Time')

subplot (4,1,2)
hist(dX,length(T)/4)
ylabel ('Frequency (dX)')
xlabel ('dX')

subplot (4,1,3)
plot (T, X)
title ('Variance Gamma log (Price) ')
ylabel ('X')
xlabel ('Time')
axis tight

eX=exp(X);
subplot (4,1,4)
plot (T,eX)
ylabel ('e^X')
title ('Variance Gamma Asset Price ')
xlabel ('Time')
axis tight

figure
stem3(T,G,X)
hold on
plot3(T,G,X,'gr')
axis tight
xlabel ('Real Time')
ylabel ('Gamma Time')
zlabel ('X')
title ('Variance Gamma Process')

% Generate same VG process in CGM parameters
% convert v,theta, sigma to C,G,M 
Cp=1/v
Gp=1/(sqrt(0.25*theta^2*v^2+0.5*sigma^2*v) - 0.5*theta*v)
Mp=1/(sqrt(0.25*theta^2*v^2+0.5*sigma^2*v) + 0.5*theta*v)
%Up and down gamma process: Gu, Gd
au=Cp; ad=Cp; 
bu=Mp; bd=Gp;
%Difference of up and down gamma process: X=Gu-Gd
for j = 1:length(T)
    dGu(j)=GammaRand (au.*dt,bu); % up Gamma process
    dGd(j)=GammaRand (ad.*dt,bd); % up Gamma process
end
Gu=cumsum(dGu);
Gd=cumsum(dGd);
VGdiff=Gu-Gd;
figure
plot (T, VGdiff, T, X, '--')
legend('CGM','\sigmav\theta', 'location','NorthWest')
ylabel ('Log Price')
xlabel ('Time')
title('Comparison of VG Notation')


