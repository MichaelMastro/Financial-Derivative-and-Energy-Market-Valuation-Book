function [ EuroCall AmerCall EuroPut AmerPut ] =...
    CrankNicAmerPDE( S0, K, T, r, vol, Smax, M, q)
% CrankNicAmerPDE employs succesive over relation 
% (SOR) on a PDE grid to calculate  price
% of American call and put options. PlotFlag Option
% turns on plotting of option value as a function
% of time and underlying asset price
% Option Delta and Gamma can be easily extracted from
% grid. A 2D plot for Delta and Gamma at time zero is 
% generated as well as a 3D plot for Delta and Gamma
% vs time and asset price

disp ('Crank Nicolson PDE: American Techniques')

close all;
clc;
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);
    
if (nargin < 8), q = 0; end
if (nargin < 7), M = 41, end % accurate -> M>101
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
tstepApprox = dS^2/sqrt(vol)/200
Napprox = round (T/tstepApprox+1);
t = linspace(0,T,Napprox);
N = length(t)
dt = t(2)-t(1)
 
maxSK = (max(S-K,0))'; 
maxKS = (max(K-S,0))'; 

j = [2:M-1]; % Stock Price Interior Nodes
jAll=[1:M]'-1; % use for S=dS*jAll
alphaCoeff = 0.25*dt*(vol^2*jAll.^2 -(r-q)*jAll);
betaCoeff = -0.5*dt*(vol^2*jAll.^2+r);
gammaCoeff = 0.25*dt*(vol^2*jAll.^2 +(r-q)*jAll);

% Same Boundary Conditions for American
iAll=[1:N]; % use for t=(iAll-1)*dt  
%tau= T-t =(N-1)*dt - (iAll-1)*dt = (N-iAll)*dt
%%%%%%%%%%%%%%%%%%%%%%%%%%
% CN technique only calculates interior nodes.
% We ignore the first and last a,b,c Coefficient
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix Division
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
fprintf(1,'\nMatrix Division: Time = %e', MatrixLeftTime);
EuroPut = interp1 (S, p(:,1), S0); 
EuroCall = interp1 (S, c(:,1), S0);
fprintf (1,'\nEuropean Put = %f, European Call = %f\n\n',...
    EuroPut, EuroCall);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOR: American
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct Calculation of American Put (P) and Call (C)
% Successive  Over Relaxation Technique is required as
% each node is repeatedly checked for early payout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = zeros(M,N); C(:,N) = maxSK; 
P = zeros(M,N); P(:,N) = maxKS;  %time T expiration
C(1,iAll) = 0; % OTM
C(M,iAll) = Smax*exp(-q*((iAll-1)*dt))...
                -K*exp(-r*(N-iAll)*dt); %ITM
P(M,iAll) = 0; % OTM
P(1,iAll) = K;%K*exp(-r*(N-iAll)*dt); %ITM %Smin=0
Cchange = zeros(M,1); Pchange = zeros(M,1);
Cold = zeros(M,1); Pold=zeros(M,1); 
rC=zeros(M,1); rP=zeros(M,1);
CEdge=zeros(M-2,1); PEdge=zeros(M-2,1); 

tol = 0.0000001;
w=1.5;
maxcount=40;
tic
jNodes=[1:M];
for i=N-1:-1:1 
    %%%% American Call Calculation
    C(j,i)=C(j,i+1);
    CEdge(end)=(C(M,i+1))*gammaCoeff(M-1); %CEdge(1)=0; 
    rC(j)= ((CoeffMat2*C(j,i+1) + CEdge));
    CchangeNorm=1;
    counter =0;
    while ((CchangeNorm > tol) && (counter < maxcount)) 
        for m=2:M-1
            Cold(m)=C(m,i);
            C(m,i)=max(maxSK(m),...
                (C(m,i)+(w.*(1./(1-betaCoeff(m))).*...
                (rC(m) + alphaCoeff(m).*C((m-1),i)...
                - (1-betaCoeff(m)).*C(m,i)...
                + gammaCoeff(m).*C(m+1,i)))));  
        end
        Cchange(j) = C(j,i) - Cold(j);
            CchangeNorm=(norm(Cchange));
            counter = counter +1;      
    end %[i counter CchangeNorm ]
     
    %%%% American Put Calculation    
    P(j,i)=P(j,i+1);
    PEdge(1)=(P(1,i+1))*alphaCoeff(2); %pEdge(end)=0;
    rP(j) = ((CoeffMat2*P(j,i+1) + PEdge));  
    PchangeNorm=1;
    counter =0;
   while ((PchangeNorm > tol) && (counter < maxcount)) 
        for m=2:M-1
            Pold(m)=P(m,i);
            P(m,i)=max(maxKS(m),...
                (P(m,i)+(w.*(1./(1-betaCoeff(m))).*...
                (rP(m) + alphaCoeff(m).*P((m-1),i)...
                - (1-betaCoeff(m)).*P(m,i)...
                + gammaCoeff(m).*P(m+1,i)))));
        end
        Pchange(j) = P(j,i) - Pold(j);
        PchangeNorm=(norm(Pchange));
        counter = counter +1;
   end % [i counter PchangeNorm ]
end

sorTime=toc;
fprintf('\nSuccesive Over Relaxation: Time = %e', sorTime);
AmerPut = interp1 (S, P(:,1), S0); 
AmerCall = interp1 (S, C(:,1), S0);
fprintf (1,'\nAmerican Put = %f, American Call = %f\n\n',...
    AmerPut, AmerCall);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Option Greeks%%%%%%%%%%%%%%%%%%%%
DeltaCallAmer=(C(3:end,1)-C(1:end-2,1))/(2*dS);
DeltaPutAmer=(P(3:end,1)-P(1:end-2,1))/(2*dS);
DeltaCallEur=(c(3:end,1)-c(1:end-2,1))/(2*dS);
DeltaPutEur=(p(3:end,1)-p(1:end-2,1))/(2*dS);

GammaCallAmer=(C(3:end,1)-2*C(2:end-1,1)+C(1:end-2,1))/(dS^2);
GammaPutAmer=(P(3:end,1)-2*P(2:end-1,1)+P(1:end-2,1))/(dS^2);
GammaCallEur=(c(3:end,1)-2*c(2:end-1,1)+c(1:end-2,1))/(dS^2);
GammaPutEur=(p(3:end,1)-2*p(2:end-1,1)+p(1:end-2,1))/(dS^2);

DeltaCallAmerAll=(C(3:end,:)-C(1:end-2,:))/(2*dS);
GammaCallAmerAll=(C(3:end,:)-2*C(2:end-1,:)+C(1:end-2,:))/(dS^2);

DeltaPutAmerAll=(P(3:end,:)-P(1:end-2,:))/(2*dS);
GammaPutAmerAll=(P(3:end,:)-2*P(2:end-1,:)+P(1:end-2,:))/(dS^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option to plot numerical solution to PDE 
% as a function of time and underlying asset price   
    PlotFlag=1;
if (PlotFlag) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compare American and European call option price 
% evolution over time 
    figure
    Smatrix=repmat(S',1,N);
    tmatrix=repmat(t,M,1);
    surf (Smatrix,tmatrix,C)
    hold on
    surf (Smatrix,tmatrix,c)
    hold off
    shading interp
    xlabel ('Asset Price'); ylabel ('Time [Years]');
    zlabel ('Call Option Price');
    title ('Crank-Nicolson: European = American Call') 
    
% compare American and European put option price 
% evolution over time    
    figure
    surf (Smatrix,tmatrix,P)
    hold on
    surf (Smatrix,tmatrix,p)
    hold off
    alpha(.7)
    shading interp
    title ('Crank-Nicolson: American and European Put') 
    xlabel ('Asset Price'); ylabel ('Time [Years]');
    zlabel ('Put Option Price'); 
    
    figure
% compare american and european call and put option delta
    subplot (2,2,1); plot (S(2:(end-1))/S0,DeltaCallAmer,'-',...
        S(2:(end-1))/S0,DeltaCallEur,'--');
    title ('                Time zero: 5 years until Expiration')
    xlabel ('Asset Price / S_0'); 
    ylabel ('Call \Delta = \deltac/\deltaS');
    subplot (2,2,3); plot (S(2:(end-1))/S0,DeltaPutAmer,'-',...
        S(2:(end-1))/S0,DeltaPutEur,'--');
    xlabel ('Asset Price / S_0'); 
    ylabel ('Put \Delta = \deltap/\deltaS');

% compare american and european call and put option gamma
   subplot (2,2,2); plot (S(2:(end-1))/S0,GammaCallAmer,'-',...
        S(2:(end-1))/S0,GammaCallEur,'--');
    xlabel ('Asset Price / S_0'); 
    ylabel ('Call \Gamma = \delta^2c/\deltaS^2');
    subplot (2,2,4); plot (S(2:(end-1))/S0,GammaPutAmer,'-',...
        S(2:(end-1))/S0,GammaPutEur,'--');
    xlabel ('Asset Price / S_0'); 
    ylabel ('Put \Gamma = \delta^2p/\deltaS^2');
    legend('American','European')


    figure
% plot call delta vs time and vs asset price 
    subplot (2,2,1)
    Smatrix=repmat((S(2:end-1))',1,N)/ S0;
    tmatrix=repmat(t,M-2,1);
    colormap pink
    contourf (Smatrix,tmatrix,DeltaCallAmerAll)
    xlabel ('Asset Price / S_0'); ylabel ('Time [Years]');
    title ('American Call \Delta = \deltaP/\deltaS');
% plot put delta vs time and vs asset price 
    subplot (2,2,3)
    colormap pink
    contourf (Smatrix,tmatrix,DeltaPutAmerAll)
    xlabel ('Asset Price / S_0'); ylabel ('Time [Years]');
    title ('American Put \Delta = \deltaP/\deltaS');

% near expiration (last time few points) explosive increase 
% in gamma obscures plot 
% plot call gamma vs time and vs asset price
    subplot (2,2,2)
    contour (Smatrix,tmatrix,GammaCallAmerAll,50)
    xlabel ('Asset Price / S_0'); ylabel ('Time [Years]');
    title ('American Call \Gamma = \delta^2P/\deltaS^2');
% plot put gamma vs time and vs asset price
    subplot (2,2,4)
    contour (Smatrix,tmatrix,GammaPutAmerAll,50)
    xlabel ('Asset Price / S_0'); ylabel ('Time [Years]');
    title ('American Put \Gamma = \delta^2P/\deltaS^2');
end
end