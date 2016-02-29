function [ EuroCall AmerCall EuroPut AmerPut ] =...
    ExplicitPDE( S0, K, T, r, vol, Smax, M, q)
% ExplicitPDE employs PDE grid to calculate  price
% of European and American option. PlotFlag Option
% turns on plotting of option value as a function
% of time and underlying asset price
% Working in linear (as apposed to log) underlying
% asset values can lead to oscillations in solution
% Must resort to small time steps to suppress 
% oscillations

close all;
clc;
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);
    
if (nargin < 8), q = 0; end
if (nargin < 7), M = 21, end
if (nargin < 6), Smax = 100; end 
if (nargin < 5), vol = 0.2; end
if (nargin < 4), r = 0.03; end
if (nargin < 3), T = 5; end
if (nargin < 2), K = 50; end
if (nargin < 1), S0 = 50; end
    
dS=Smax/(M-1)
S=[0:dS:Smax]';
% need very small dt relative to dS (and volatility)
% else severe oscillations
% One simple step size calibration is
tstepApprox = dS^2/sqrt(vol)/1000
Napprox = round (T/tstepApprox+1)
t = linspace(0,T,Napprox);
N = length(t)
dt = t(2)-t(1)

c = zeros(M,N); C = zeros(M,N);
maxSK = (max(S-K,0))'; 
c(:,N) = maxSK; C(:,N) = maxSK; %time T expiration

p = zeros(M,N); P = zeros(M,N);
maxKS = (max(K-S,0))'; 
p(:,N) = maxKS; P(:,N) = maxKS; %time T expiration

j = [2:M-1]'; % Stock Price Interior Nodes
% Coefficients for linear explicit PDE
aCoeffStar = (1/(1+r*dt))*...
    (-0.5*(r-q)*j*dt + 0.5*vol^2*j.^2*dt);
bCoeffStar = (1/(1+r*dt))*(1 - vol^2*j.^2*dt);
cCoeffStar = (1/(1+r*dt))*...
    (0.5*(r-q)*j*dt + 0.5*vol^2*j.^2*dt);

for i=N-1:-1:1
% Time step for European (p) and American (P) put  
    c(j,i) = aCoeffStar.*c(j-1,i+1) + bCoeffStar.*c(j,i+1)...
                +cCoeffStar.*c(j+1,i+1);
    c(1,i) = 0; % OTM
    c(M,i) = Smax*exp(-q*((i-1)*dt))-K*exp(-r*(N-i)*dt); %ITM
    % Compare European Value vs. Early Exercise
    C(:,i) = max ( c(:,i), (S-K) ); % S=(jAll-1)*dS
% Repeat Time step for European (p) and American (P) put    
    p(j,i) = aCoeffStar.*p(j-1,i+1) + bCoeffStar.*p(j,i+1)...
                +cCoeffStar.*p(j+1,i+1);
    p(M,i) = 0; % OTM
    p(1,i) = K*exp(-r*(N-i)*dt); %ITM %Smin=0
    P(:,i) = max ( p(:,i), (K-S) );
end

% Option to plot numerical solution to PDE 
% as a function of time and underlying asset price
PlotFlag=1;
if (PlotFlag) 
    Smatrix=repmat(S,1,N);
    tmatrix=repmat(t,M,1);
    surf (Smatrix,tmatrix,c)
    shading interp
    xlabel ('Asset Price'); ylabel ('Time [Years]');
    zlabel ('Call Option Price');
    title ('American Call = European Call') 

    figure
    surf(Smatrix,tmatrix,P)
    hold on
    surf (Smatrix,tmatrix,p)
    surf (Smatrix,tmatrix,P-p)
    alpha(.7)
    shading interp
    hold off
    title ('Added Value of American vs. European Put') 
    xlabel ('Asset Price'); ylabel ('Time [Years]');
    zlabel ('Put Option Price');
end

% Interpolate as asset grid point may not match 
% initial Asset price S0
disp('European Put')
EuroPut = interp1 (S, p(:,1), S0); disp (EuroPut)
disp('American Put')
AmerPut = interp1 (S, P(:,1), S0); disp (AmerPut)
% Hull (2006) suggests control variate technique to  
% improve American Put price accuracy. Essential idea 
% is to compare PDE and analytical PDE. The difference
% is used to adjust American put value

% European call should equal american call value
% unless numerical error is present
disp('European Call')
EuroCall = interp1 (S, c(:,1), S0); disp (EuroCall)
disp('American Call')
AmerCall = interp1 (S, C(:,1), S0); disp (AmerCall)
end