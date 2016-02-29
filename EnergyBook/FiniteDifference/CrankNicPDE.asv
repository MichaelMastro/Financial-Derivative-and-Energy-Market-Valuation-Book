function [ EuroCall AmerCall EuroPut AmerPut ] =...
    CrankNicPDE( S0, K, T, r, vol, Smax, M, q)
% CrankNicPDE employs PDE grid to calculate  price
% of European call and put options. PlotFlag Option
% turns on plotting of option value as a function
% of time and underlying asset price
% Crank Nicolson improves convergence 
% (vs. explicit PDE) but is less stable 
% compared to Implicit PDE 
% Three techniques Inversion, Matrix Left Division,  
% LU decomposition) are compared for accuracy and time
% By computational count, the LU decomposition 
% should be the fastest but isn't noticeable
% This is attributable to the very efficient Matlab 
% built-in functions  
% A fast Direct Gaussian Elimination Approach is
% available and simple to code thanks to the 
% sparse, tridiagonal matrix

disp ('Crank Nicolson PDE')

close all;
clc;
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);
    
if (nargin < 8), q = 0; end
if (nargin < 7), M = 41, end % Accurate use M>101
if (nargin < 6), Smax = 100; end 
if (nargin < 5), vol = 0.2; end
if (nargin < 4), r = 0.03; end
if (nargin < 3), T = 5; end
if (nargin < 2), K = 50; end
if (nargin < 1), S0 = 50; end

dS=Smax/(M-1) % equivalent to M=Smax/dS + 1 
%Smax= (M-1)*dS % S=0,dS,2*dS,...,(M-1)*dS
S=0:dS:Smax;
% simple step size calibration of
tstepApprox = dS^2/sqrt(vol)/40
Napprox = round (T/tstepApprox+1);
t = linspace(0,T,Napprox);
N = length(t)
dt = t(2)-t(1)

c = zeros(M,N); %C = zeros(M,N);
maxSK = max(S-K,0); 
c(:,N) = maxSK; %C(:,N) = maxSK; %time T expiration

p = zeros(M,N); %P = zeros(M,N);
maxKS = max(K-S,0); 
p(:,N) = maxKS; %P(:,N) = maxKS; %time T expiration

j = [2:M-1]; % Stock Price Interior Nodes
jAll=[1:M]-1; % use for S=dS*jAll
alphaCoeff = 0.25*dt*(vol^2*jAll.^2 -(r-q)*jAll);
betaCoeff = -0.5*dt*(vol^2*jAll.^2+r);
gammaCoeff = 0.25*dt*(vol^2*jAll.^2 +(r-q)*jAll);
%%%%%%%
iAll=[1:N]; % use for t=(iAll-1)*dt  
%tau= T-t =(N-1)*dt - (iAll-1)*dt = (N-iAll)*dt
c(1,iAll) = 0; % OTM
c(M,iAll) = Smax*exp(-q*((iAll-1)*dt))...
                -K*exp(-r*(N-iAll)*dt); %ITM
p(M,iAll) = 0; % OTM
p(1,iAll) = K*exp(-r*(N-iAll)*dt); %ITM %Smin=0

% CN technique only calculates interior nodes.
% We throw out the first and last a,b,c Coefficient
CoeffMat1=diag(-alphaCoeff(3:M-1),-1)+...
    diag(1-betaCoeff(2:M-1))+diag(-gammaCoeff(2:M-2),1);
CoeffMat2=diag(alphaCoeff(3:M-1),-1)+...
    diag(1+betaCoeff(2:M-1))+diag(gammaCoeff(2:M-2),1);
% As just mentioned, CN technique only calculates 
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
invCoeffMat1=inv(CoeffMat1);
for i=N-1:-1:1
    cEdge(end)=(c(M,i)+c(M,i+1))*gammaCoeff(M-1); 
    %cEdge(1)=0;
    c(j,i)= (invCoeffMat1 *(CoeffMat2*c(j,i+1) + cEdge)); 
    pEdge(1)=(p(1,i)+p(1,i+1))*alphaCoeff(2); 
    %pEdge(end)=0;
    p(j,i) = (invCoeffMat1 *(CoeffMat2*p(j,i+1) + pEdge));   
end
MatrixInversionTime=toc;
fprintf (1,'\nMatrix Inversion: Time = %e',...
    MatrixInversionTime);
EuroPut = interp1 (S, p(:,1), S0); 
EuroCall = interp1 (S, c(:,1), S0);
fprintf (1,'\nEuropean Put = %f; European Call = %f\n',...
    EuroPut, EuroCall);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix LU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[L,U]=lu(CoeffMat1);
for i=N-1:-1:1
    cEdge(end)=(c(M,i)+c(M,i+1))*gammaCoeff(M-1); 
    % cEdge(1)=0;
    c(j,i)= (U \ (L \(CoeffMat2*c(j,i+1) + cEdge))); 
    pEdge(1)=(p(1,i)+p(1,i+1))*alphaCoeff(2); 
    % pEdge(end)=0;
    p(j,i) = (U \ (L \(CoeffMat2*p(j,i+1) + pEdge)));   
end
MatrixLUtime=toc;
fprintf (1,'\nMatrix LU: Time = %e',...
    MatrixLUtime)
EuroPut = interp1 (S, p(:,1), S0); 
EuroCall = interp1 (S, c(:,1), S0);
fprintf (1,'\nEuropean Put = %f; European Call = %f\n',...
    EuroPut, EuroCall);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix Division
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for i=N-1:-1:1  
    cEdge(end)=(c(M,i)+c(M,i+1))*gammaCoeff(M-1); 
    % cEdge(1)=0;
    c(j,i)= (CoeffMat1 \ (CoeffMat2*c(j,i+1) + cEdge));  
    pEdge(1)=(p(1,i)+p(1,i+1))*alphaCoeff(2); 
    % pEdge(end)=0;
    p(j,i) = (CoeffMat1 \ (CoeffMat2*p(j,i+1) + pEdge));   
end
MatrixLeftTime=toc;
fprintf (1,'\nMatrix Division: Time = %e', MatrixLeftTime);
EuroPut = interp1 (S, p(:,1), S0); 
EuroCall = interp1 (S, c(:,1), S0);
fprintf (1,'\nEuropean Put = %f, European Call = %f\n',...
    EuroPut, EuroCall);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gaussian Elimination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = zeros(M,N); c(:,N) = maxSK; 
p = zeros(M,N); p(:,N) = maxKS;  %time T expiration
c(1,iAll) = 0; % OTM
c(M,iAll) = Smax*exp(-q*((iAll-1)*dt))...
                -K*exp(-r*(N-iAll)*dt); %ITM
p(M,iAll) = 0; % OTM
p(1,iAll) = K*exp(-r*(N-iAll)*dt); %ITM %Smin=0
cEdge=zeros(M-2,1);
pEdge=zeros(M-2,1);

% place coeffecients into form to use standard
% tridiagonal matrix algorithm of Thomas
%aCoeff(j)c(j-1)+bCoeff(j)c(j)+aCoeff(j)c(j+1)=r(j)
aCoeff=-alphaCoeff;
bCoeff=1-betaCoeff;
cCoeff=-gammaCoeff;

tic
for i=N-1:-1:1  
    cEdge(end)=(c(M,i)+c(M,i+1))*gammaCoeff(M-1); 
    % cEdge(1)=0;
    rc(j)=  (CoeffMat2*c(j,i+1) + cEdge) ;
    pEdge(1)=(p(1,i)+p(1,i+1))*alphaCoeff(2);
    % pEdge(end)=0;
    rp(j) =  (CoeffMat2*p(j,i+1) + pEdge);
    
    %% modify coeffecients and right hand side
    cCoeffHat(2)=cCoeff(2)/bCoeff(2);
    rcHat(2)=rc(2)/bCoeff(2);
    rpHat(2)=rp(2)/bCoeff(2);
    for m=3:1:M-1 
        cCoeffHat(m)=cCoeff(m)/...
            (bCoeff(m)-cCoeffHat(m-1)*aCoeff(m));
        rcHat(m)=(rc(m)-rcHat(m-1)*aCoeff(m)) /...
            (bCoeff(m) - cCoeffHat(m-1)*aCoeff(m));
        
        rpHat(m)=(rp(m)-rpHat(m-1)*aCoeff(m)) /...
            (bCoeff(m) - cCoeffHat(m-1)*aCoeff(m));
    end
        c(M-1,i)=rcHat(M-1);
        p(M-1,i)=rpHat(M-1);
    for m=M-2:-1:2
        c(m,i)=rcHat(m)-cCoeffHat(m)*c(m+1,i);
        p(m,i)=rpHat(m)-cCoeffHat(m)*p(m+1,i);
    end
   
end
GaussElimTime=toc;
fprintf (1,'\nGaussian Elimination: Time = %e', GaussElimTime);
EuroPut = interp1 (S, p(:,1), S0); 
EuroCall = interp1 (S, c(:,1), S0);
fprintf (1,'\nEuropean Put = %f, European Call = %f\n',...
    EuroPut, EuroCall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analytical Black-Scholes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
bsCall = BlackScholesCall (K,S0,T,vol,r,0);
bsPut = BlackScholesPut (K,S0,T,vol,r,0);
bsTime=toc;
fprintf (1,'\nAnalytical Black Scholes: Time = %e', bsTime);
fprintf (1,'\nEuropean Put = %f; European Call = %f\n',...
    bsPut, bsCall);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Option to plot numerical solution to PDE 
% as a function of time and underlying asset price
PlotFlag=1;
if (PlotFlag) 
    Smatrix=repmat(S',1,N);
    tmatrix=repmat(t,M,1);
    figure
    surf (Smatrix,tmatrix,c)
    shading interp
    xlabel ('Asset Price'); ylabel ('Time [Years]');
    zlabel ('Call Option Price');
    title ('Crank-Nicolson European Call') 
     
    figure
    surf (Smatrix,tmatrix,p)
    alpha(.7)
    shading interp
    title ('Crank-Nicolson European Put') 
    xlabel ('Asset Price'); ylabel ('Time [Years]');
    zlabel ('Put Option Price');
end
end