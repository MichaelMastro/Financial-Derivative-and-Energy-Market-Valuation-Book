function [ EuroCall AmerCall EuroPut AmerPut ] =...
    ExplicitLogPDE ( S0, K, T, r, vol, M, q)
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
if (nargin < 6), Smax = 150; end 
if (nargin < 5), vol = 0.2; end
if (nargin < 4), r = 0.03; end
if (nargin < 3), T = 5; end
if (nargin < 2), K = 50; end
if (nargin < 1), S0 = 50; end
    
dX=log(Smax/S0)/(M/2-1)
%dt=0.1
%dX=vol*sqrt(3*dt)
X=linspace((-dX*(M-1)/2) , (dX*(M-1)/2) , M);
%Mapprox=length(X)
S=S0*exp(X');
disp ('Smin=');disp(S(1));
disp ('Smax=');disp(S(end));
dt=(dX)^2/(3*vol)
% S=S0*exp(dX); dX log spacing is more stable PDE
% dX=vol*sqrt(3*dt) -> dt=dX^2/(3*vol)
t = 0:dt:T;
N = length(t)

c = zeros(M,N); C = zeros(M,N);
maxSK = (max(S-K,0))'; 
c(:,N) = maxSK; C(:,N) = maxSK; %time T expiration

p = zeros(M,N); P = zeros(M,N);
maxKS = (max(K-S,0))'; 
p(:,N) = maxKS; P(:,N) = maxKS; %time T expiration

j = [2:M-1]; % log Stock Price Interior Nodes
% Coefficients for explicit log-price PDE
% j term is not present
alphaCoeffStar = (1/(1+r*dt))*...
    (-0.5*(r-q-0.5*vol^2)*dt/dX + 0.5*vol^2*dt/dX^2);
betaCoeffStar = (1/(1+r*dt))*(1 - vol^2*dt/dX^2);
gammaCoeffStar = (1/(1+r*dt))*...
    (0.5*(r-q-0.5*vol^2)*dt/dX + 0.5*vol^2*dt/dX^2);

for i=N-1:-1:1
% Time step for European (p) and American (P) put  
    c(j,i) = alphaCoeffStar.*c(j-1,i+1)...
                + betaCoeffStar.*c(j,i+1)...
                + gammaCoeffStar.*c(j+1,i+1);
    c(1,i) = 0; % OTM
    %ITM
    c(M,i) = Smax*exp(-q*((i-1)*dt))-K*exp(-r*(N-i)*dt); 
    % Compare European Value vs. Early Exercise
    C(:,i) = max ( c(:,i), (S-K) );
% Repeat Time step for European (p) and American (P) put    
    p(j,i) = alphaCoeffStar.*p(j-1,i+1)...
                + betaCoeffStar.*p(j,i+1)...
                + gammaCoeffStar.*p(j+1,i+1);
    p(M,i) = 0; % OTM
    %ITM %Smin~0
    p(1,i) = K*exp(-r*(N-i)*dt)-S(1)*exp(-q*((i-1)*dt)); 
    P(:,i) = max ( p(:,i), (K-S) );
end

% Option to plot numerical solution to PDE 
% as a function of time and underlying asset price
PlotFlag=1;
if (PlotFlag) 
    Smatrix=repmat(S,1,N);
    tmatrix=repmat(t,M,1);
      sizeSmatrix=size(Smatrix)
    sizetmatrix=size(tmatrix)
    sizec=size(c)
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
    alpha(.7); shading interp;
    hold off
    title ('Added Value of American vs. European Put') 
    xlabel ('Asset Price'); ylabel ('Time [Years]');
    zlabel ('Put Option Price');
end

% option value at current asset price should be at
% center of grid at time zero 
disp('European Put')
EuroPut = p(((M-1)/2),1); disp (EuroPut)
disp('American Put')
AmerPut = P(((M-1)/2),1); disp (AmerPut)
% Hull (2006) suggests control variate technique to  
% improve American Put price accuracy. Essential idea 
% is to compare PDE and analytical PDE. The difference
% is used to adjust American put value

% European call should equal american call value
% unless numerical error is present
disp('European Call')
EuroCall = c(((M-1)/2),1);  disp (EuroCall)
disp('American Call')
AmerCall = C(((M-1)/2),1); disp (AmerCall)
end