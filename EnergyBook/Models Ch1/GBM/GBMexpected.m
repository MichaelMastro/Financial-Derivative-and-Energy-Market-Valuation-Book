function GBMexpected ()
%GBMexpected.m simulates GBM (no input)
%of several (vectorized) price paths. 
%Goal is to show that the expected value of GBM 
%E[St]=S0exp(mu*t) is equal to the average of many price paths
%generated by loop St=S0*exp[(mu-0.5sig^2)t+sig*sqrt(t)*N(0,1) 
%Normally generated random movements N(0,1) are symmetric but
%exp(N(0,1)) movements asymmetric. 

close all
clear all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

TimeDelta=1/252; %assume daily prices
steps = 7.5*252; %assume 7.5 years
TimeLength=TimeDelta*steps; 
time = linspace(0,TimeLength,steps);  %years

loop=5000; %number of simulations 
mu=0.1; %percentage drift
Dsig =0.1;%annualized volatility 
Szero=30; %Initial Price

etadt=(mu-0.5*Dsig^2)*TimeDelta; %result of Ito
Expected=zeros(1,steps); SD=zeros(1,steps);
Expected = Szero.*exp((mu).*(time));     %E[St]=S0exp(mu*t)
SD = Expected.*sqrt(exp(Dsig^2.*time) -1); %Stand. Dev.

subplot (1,2,1)
plot(time,Expected,'-','color','blue','LineWidth',5);hold on;
plot(time,Expected+SD,'color','black');hold on;
plot(time,Expected-SD,'color',[0.5 0 0]); hold on;

%vector of random movements
sigma = Dsig*sqrt(TimeDelta).*randn(loop,steps);
S=zeros(loop,steps); 
S=Szero*cumprod(exp(etadt+sigma),2); 
%cumsum(logS0,etadt+sigma) followed by S=exp
%should be faster since addition is generally faster than mult.
aveS=mean(S);

plot(time,aveS,':','color','magenta','LineWidth',5);hold on;
iteration = 1:50:500; %Plot only a few price series
plot(time,S(iteration,:),':','color','cyan','LineWidth',0.1); 
hold on; 
    
axis tight
xlabel('Time [Years]'); ylabel('Asset Price')%%
legend('E[S_t]=e^\mu^t', 'E[S_t] + 1 Stand. Dev.',...
    'E[S_t] - 1 Stand. Dev.','Mean Sim. Price',...
    'S_t=S_0e^{(\mu-0.5\sigma^2)t+\sigmaN(0,1)t^{0.5}}',...
    'location','NorthWest');
title('Geometric Brownian Motion Simulated Price Process');
hold off;

Send=S(:,steps);
minS=min(Send); maxS=max(Send); 
subplot(1,2,2)
inc=minS:0.2:maxS; %define bin edges
hist(S(:,steps),inc);
xlabel('Final Price'); ylabel('Frequency'); 
title('Histogram of Final Simulated Price');
text(Expected(end),loop/180,'\downarrow E[S_T]');
axis tight
end
