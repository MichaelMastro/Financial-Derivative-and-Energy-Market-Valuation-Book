function [ EuroCall AmerCall EuroPut AmerPut ] =...
    ImplicitPDE( S0, K, T, r, vol, Smax, M, q)
% ImplicitPDE employs PDE grid to calculate  price
% of European call and put options. PlotFlag Option
% turns on plotting of option value as a function
% of time and underlying asset price
% Implicit (vs. explicit) PDE is more stable and should
% require fewer steps in grid. The tradeoff is an 
% increase in computational time. Three techniques 
% (Inversion, Matrix Left Division, and 
% LU decomposition) are compared for accuracy and time
% By computational count, the LU decomposition should be
% the fastest. 

close all;
clc;
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);
    
if (nargin < 8), q = 0; end
if (nargin < 7), M = 101, end
if (nargin < 6), Smax = 100; end 
if (nargin < 5), vol = 0.2; end
if (nargin < 4), r = 0.03; end
if (nargin < 3), T = 5; end
if (nargin < 2), K = 50; end
if (nargin < 1), S0 = 50; end

dS=Smax/(M-1) % equivalent to M=Smax/dS + 1 
%Smax= (M-1)*dS % S=0,dS,2*dS,...,(M-1)*dS
S=0:dS:Smax;
dStest=S(2)-S(1)
% Use simple step size calibration of
tstepApprox = dS^2/sqrt(vol)/40
Napprox = round(T/tstepApprox) + 1 
% add 1 for zero point
%e.g., T=2yr, dt=1yr, N=3, t=0,1=0,dt,(N-1)*dt
t = linspace(0,T,Napprox);
N = length(t)
% Since Napprox rounded, need to find new dt
dt = t(2)-t(1)

c = zeros(M,N); 
maxSK = max(S-K,0); 
c(:,N) = maxSK; %time T expiration

p = zeros(M,N);
maxKS = max(K-S,0); 
p(:,N) = maxKS;  %time T expiration

j = [2:M-1]; % Stock Price Interior Nodes
jAll=[1:M]-1; % use for S=dS*jAll
aCoeff = 0.5*(r-q)*jAll*dt - 0.5*vol^2*jAll.^2*dt;
bCoeff = 1 + vol^2*jAll.^2*dt + r*dt;
cCoeff = -0.5*(r-q)*jAll*dt - 0.5*vol^2*jAll.^2*dt;
%%%%%%%

iAll=[1:N]; % use for t=(iAll-1)*dt  
%tau= T-t =(N-1)*dt - (iAll-1)*dt = (N-iAll)*dt
c(1,iAll) = 0; % OTM
c(M,iAll) = Smax*exp(-q*((iAll-1)*dt))...
                -K*exp(-r*(N-iAll)*dt); %ITM
p(M,iAll) = 0; % OTM
p(1,iAll) = K*exp(-r*(N-iAll)*dt); %ITM %Smin=0

% Implicit technique only calculates interior nodes.
% We throw out the first and last a,b,c Coefficient
CoeffMat=diag(aCoeff(3:M-1),-1)+...
    diag(bCoeff(2:M-1))+diag(cCoeff(2:M-2),1);

% As just mentioned, implicit technique only calculates 
% interior nodes. The top interior node needs information
% from the top edge node and the bottom interior node needs
% information from the bottom edge node. Since we can't 
% extend the coefficient matrix, the edge nodes 
% aCoeff(2) and cCoeff(M-1) have to be added in 
cEdge=zeros(M-2,1);
pEdge=zeros(M-2,1);

% Matrix Inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
invCoeffMat=inv(CoeffMat);
for i=N-1:-1:1
    cEdge(end)=-c(M,i)*cCoeff(M-1); %cEdge(1)=0;
    c(j,i)= (invCoeffMat * ( c(j,i+1) + cEdge) ) ;
    pEdge(1)=-p(1,i)*aCoeff(2); % pEdge(end)=0;
    p(j,i) = (invCoeffMat * ( p(j,i+1) + pEdge));   
end
MatrixInversionTime=toc;
sprintf('Matrix Inversion: Time = %e',...
    MatrixInversionTime)
EuroPut = interp1 (S, p(:,1), S0); 
EuroCall = interp1 (S, c(:,1), S0);
sprintf ('European Put = %f; European Call Put = %f',...
    EuroPut, EuroCall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix LU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[L,U]=lu(CoeffMat);
for i=N-1:-1:1
    cEdge(end)=-c(M,i)*cCoeff(M-1); %cEdge(1)=0;
    c(j,i)= (U \ (L \ ( c(j,i+1) + cEdge)) ); 
    pEdge(1)=-p(1,i)*aCoeff(2); % pEdge(end)=0;
    p(j,i) = (U \ (L \  ( p(j,i+1) + pEdge)) );   
end
MatrixLUtime=toc;
sprintf('Matrix LU: Time = %e',...
    MatrixLUtime)
EuroPut = interp1 (S, p(:,1), S0); 
EuroCall = interp1 (S, c(:,1), S0);
sprintf ('European Put = %f; European Call Put = %f',...
    EuroPut, EuroCall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix Division
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for i=N-1:-1:1  
    cEdge(end)=-c(M,i)*cCoeff(M-1); % cEdge(1)=0;
    c(j,i)= (CoeffMat \ ( c(j,i+1) + cEdge) );  
    pEdge(1)=-p(1,i)*aCoeff(2); % pEdge(end)=0;
    p(j,i) = (CoeffMat \ ( p(j,i+1) + pEdge));   
end
MatrixLeftTime=toc;
sprintf('Matrix Division: Time = %e', MatrixLeftTime)
EuroPut = interp1 (S, p(:,1), S0); 
EuroCall = interp1 (S, c(:,1), S0);
sprintf ('European Put = %f; European Call Put = %f',...
    EuroPut, EuroCall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analytical Black-Scholes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
bsCall = BlackScholesCall (K,S0,T,vol,r,0);
bsPut = BlackScholesPut (K,S0,T,vol,r,0);
bsTime=toc;
sprintf('Analytical Black Scholes: Time = %e', bsTime)
sprintf ('European Put = %f; European Call = %f',...
    bsPut, bsCall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Option to plot numerical solution to PDE 
% as a function of time and underlying asset price
PlotFlag=1;
if (PlotFlag) 
    Smatrix=repmat(S',1,N);
    tmatrix=repmat(t,M,1);
    surf (Smatrix,tmatrix,c)
    shading interp
    xlabel ('Asset Price'); ylabel ('Time [Years]');
    zlabel ('Call Option Price');
    title ('Implicit PDE: European Call') 
     
    figure
    surf (Smatrix,tmatrix,p)
    alpha(.7)
    shading interp
    title ('Implicit PDE: European Put') 
    xlabel ('Asset Price'); ylabel ('Time [Years]');
    zlabel ('Put Option Price');
end

end