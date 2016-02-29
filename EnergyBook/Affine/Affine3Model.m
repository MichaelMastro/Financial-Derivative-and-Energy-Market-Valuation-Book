function Affine3Model(data);
% Affine3Model Implements 3 term Affine Jump-Diffusion model
% using Kalman filter via fminsearch
% Jump-Diffusion Model and approach is based on 
% N.K. Nomikos, O.A. Soldatos(2010)Energy Policy 38,5671
% who applied model to Nordpool Electricity Prices

% Their model is unique in that the jump model (Y) is an
% independent third component in their affine model
%   F=f(t)+exp(Chi+Epsilon+Y) 
% where Chi is the short term deviation mean reversion 
%   dChi=-kappaChi*Chi*dt + sigmaChi*dWchi
% Epsilon is the long term drift arithmetic Brownian motion
%   dEps=(muEps-0.5sigmaEps^2)dt + sigmaEps*dWeps
% Chi is the jump process + mean reversion  
%   dY=-kappaY*Y*dt + J(muJ,sigmaJ^2)dq(lambda)

% Nomikos suggests pre-filtering the jumps 
% from the Spot Data to calculate the jump process
% parameters and to extract the jumps from the Spot Data
% The Kalman filter is feed a deseaonalied set of futures
% and jumpless spot data

% For the electricity market, Nomikos also added a seasonal
% dependency into the jump parameters. These featutes will 
% be ignored in this analysis of the natural gas market

% f(t) is a seasonality parameter: the graphs below
% show that seasonality is surprising small compared to the
% price of impact of supply changes (not shown) for the natural
% gas data of 1997-2012. 

close all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

global tau
global Obs
global knownParam
global latent
global time
global S
global MRvec

% Price Data 'ngSpotFut.dat' contains columns of time, spot price,
% 1, 2, 3, and 4 month Futures price data where the futures 
% settlement is at the end of the particular month
if (nargin == 0), 
    data=load('ngSpotFut.dat');
end

[Col,Rows] = size(data)
if (Col>Rows)
  disp('flip')
  data=data'; %flip S=(column) to S'=(rows)
  temp=Col; Col=Rows; Rows=temp
end
time = data(1,:)+1997;
S = data(2,:);
F = data(3:Col,:);
SMat=S'*ones (1,Col-2);
offset=4; % adjust if start date not at beginning of month
          % and end date not a beginning of month
Days=1:Rows;
% Simplify by assuming 21(=252/12) trading days per (every) month
tau=mod(-Days-offset,21)'; 
taustart=tau(1)
    %   S  1M-Fut  2M-Fut  3M-Fut     4M-Fut
tau = [tau*0, tau, tau+21, tau+2*21, tau+3*21]/252;

% Check Seasonality in Natural Gas Data from 1997 to 2012 
options=optimset('MaxIter',200);%'Display','iter',
SeasonPar =[-0.12 0.6816 -0.15 0 (mean(S))] % Initial Guess
% Nomikos suggests deterministic fit of historic spot date
% with two sinusoidal functions. Function Season.m contains
% SeasonVal = amp1*sin(2pi(t+phase1)) + amp2*sin(2pi(t+phase2))
SeasonPar=fminsearch('SineFit',SeasonPar,options);
amp1=SeasonPar(1)
phase1=SeasonPar(2)
amp2=SeasonPar(3)
phase2=SeasonPar(4)
LongRunMean=SeasonPar(5)
SeasonVal=Season(SeasonPar,time);
% Substrate Seasonal Oscillations from Spot and Futures
S=S-SeasonVal; 
F(1,:)=F(1,:)-Season(SeasonPar,time+tau(1)); 
F(2,:)=F(2,:)-Season(SeasonPar,time+tau(2)); 
F(3,:)=F(3,:)-Season(SeasonPar,time+tau(3)); 
F(4,:)=F(4,:)-Season(SeasonPar,time+tau(4)); 

%% Calculate Standard Deviation as if logarithmic returns 
% of log(St/S0) are normally distributed 
LogDelta=log(S(2:end))-log(S(1:end-1)); %
lenghtLD=length(LogDelta);
M1=mean(LogDelta);   %Technically this is the first raw moment 
StanDev = std(LogDelta)
neg3sd=M1-3*StanDev;
pos3sd=M1+3*StanDev;

% Returns +/- 3 standard deviations followed by strong reversion
% are considered jump 'Y spikes' by Nomikos
% Jumps without reversion are not considered 'Y spikes'
% Nomikos Suggests Repeating until all data +/- 3 sd 
NoJumpS=S; % Loop below will remove spikes
JumpS=zeros(1,Rows)*NaN; % Loop below will add spikes
LogDeltaNJ=zeros(1,Rows); % Return Vector with Spikes removed
JumpInd=[]; % Will Contain Spike Indices
JumpReversionVec = []; % Will Contain Spike Reversion
JumpSizeVec= []; % Will Contain Spike Size
for i= 2:Rows-1
    LogDeltaNJ(i)=log(NoJumpS(i))-log(NoJumpS(i-1)); %
    LogDeltaNJ(i+1)=log(NoJumpS(i+1))-log(NoJumpS(i)); %
if (((LogDeltaNJ(i)<neg3sd) && (LogDeltaNJ(i+1)>(pos3sd/2)))...
 ||(LogDeltaNJ(i) > pos3sd ) && (LogDeltaNJ(i+1) < (neg3sd/2)))
    JumpS(i) = S(i);
    NoJumpS(i)=NoJumpS(i-1);
    JumpInd=[JumpInd i];
    JumpReversionVec =[JumpReversionVec LogDeltaNJ(i+1)];
    JumpSizeVec= [JumpSizeVec LogDeltaNJ(i)];
%else % JumpS(i)=NaN;
    end
end

kappaY=sqrt(mean(JumpReversionVec.^2))
lambda=length(JumpInd)/ time(end)
muJ=mean(abs(JumpSizeVec))
sigmaJ=sqrt(mean(JumpSizeVec.^2))


figure

plot (time, SeasonVal+LongRunMean, time, S,'k--',...
    time, JumpS, 'bs',time, NoJumpS,'g-')
title ('Daily NYMEX Natural Gas Spot Price')
xlabel ('Time [Year]'); 
ylabel('Price [Dollars per Thousand Cubic Feet]');
axis tight
legend('Deterministic Seasonality', 'Deseasonalized Price',...
    'Spikes', 'w/o Spikes')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Obs=log([NoJumpS; F(1:end,:)]); % Obserations Data for Kalman
% Check for and replace missing Data (= NaN)
[a,b]=find(isnan(Obs));
for loop = 1:length(a)
    Obs(a(loop),b(loop))=Obs(a(loop),b(loop)-1);
end

R=diag(1e-3*ones(1,Col-1)); % Observation Noise Covariance

knownParam(1) = kappaY;   %Pre-Estimated Jump Process
knownParam(2) = lambda;
knownParam(3) = muJ;
knownParam(4) = sigmaJ;
dt = 1/252;
knownParam(5) = dt;   % e.g., 1/252 for daily
knownParam(6) = Col-1;  % Col-1= Spot + # of Futures
knownParam(7) = Rows;    % number of time points

%Supply initial guess values
muChiGuess = 0.95;
kappaChiGuess = 3;
sigChiGuess = 0.34;
muEGuess = 0;
kappaEGuess = 0.27;
sigEGuess = 0.32;
rhoGuess = 0.0;
phiChiGuess=0.65;
RvarGuess=1e-3;

param=zeros(7,1);
param(1)=muChiGuess;
param(2)=kappaChiGuess;
param(3)=sigChiGuess;
param(4)=muEGuess;
param(5)=kappaEGuess;
param(6)=sigEGuess;
param(7)=rhoGuess;
param(8)=phiChiGuess;
param(9:13)=RvarGuess;    % diagonal variance terms
% fminsearch for (negative) makimum likelihood of parameters 
% as calculated by Kalman filtration of all observed data
options=optimset('Display','iter','MaxIter',100);%options=[];%
pnew=fminsearch('A3Mloglikelihoodfn',param,options);
    muChi = pnew(1)
    kappaChi = pnew(2)
    sigChi = pnew(3)
    muE = pnew(4)
    kappaE = pnew(5)
    sigE = pnew(6)
    rho = pnew(7)
    phiChi = pnew(8)
fprintf(1, 'Measurement Error Variance Terms  \n');
disp (pnew(9:13));
%%% Plot Spot Data with Jumps removed Vs. 
% Calculated Equilibirum (Epsilon) 
% and Deviation (Chi) + Equilibirum (Epsilon) 
figure 
plot (time, log(NoJumpS), '--',...
    time, (latent(:,1)+latent(:,2)),'-.',...
    time, (latent(:,2)), time, (latent(:,1)),':') 
legend('Spot w/o Spikes', '(\chi+\epsilon)',...
    '(Equilibrium:\epsilon)','(Deviation:\chi)')
xlabel ('Time [Year]'); ylabel('Log Price');
title ('Fit to Daily Natural Gas Log Price')
axis tight;

%%%%%%%%% Predict Futures Curve
ChiCurrent=latent(end,1)
EpsCurrent=latent(end,2)
YCurrent=JumpS(end)
if (isnan(YCurrent)), YCurrent=0;, end

FutT=((0:251)/252);
djump=FutT*0;
maxj=5; jcount=linspace(0,maxj,maxj);
    deltaTau=FutT/maxj;
    for j=1:maxj
        TauChop=(j-0.5)*deltaTau;
        djump=djump+(exp(muJ*exp(-kappaY*TauChop)...
            +0.5*sigmaJ^2*exp(-2*kappaY*TauChop))-1);
    end
    djump=lambda*djump.*deltaTau;

AT= (sigChi^2/(4*kappaChi))*(1-exp(-2.*kappaChi.*FutT))...
 -((phiChi-rho*sigE*sigChi)/kappaChi)*(1-exp(-kappaChi.*FutT))...
    +muE*FutT;
exp(EpsCurrent + ChiCurrent)
FutF=exp(EpsCurrent + ChiCurrent*exp(-kappaChi*FutT)...
    + YCurrent*exp(-kappaY*FutT)+...
    + djump + AT);

timeS=time(end)
FutT=FutT+time(end);
FutFpSeason=FutF+Season(SeasonPar,FutT);

figure
plot(FutT, FutF, '--' ,FutT, FutFpSeason)
title ('Predicted Natural Gas Futures Curve')
xlabel ('Time [Year]'); 
ylabel('Price [Dollars per Thousand Cubic Feet]');
axis tight
legend ('Deseaonalized Futures', 'Futures Curve')

figure
ribbon (time, MRvec)
title ('Kalman Filter Measurement Residual')
ylabel ('Time [Year]'); axis tight
zlabel ('Log Price Residual'); 
legend ('Spot', '1M Futures', '2M Futures', '3M Futures',...
    '4M Futures') 
end