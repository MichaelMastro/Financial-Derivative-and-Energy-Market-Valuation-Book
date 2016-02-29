%Geometric Brownian Motion with Jumps
%For a more detailed implementation see Floyd Hanson's
%"Stochastic Process and Control for Jump-Diffusions (2007)"
%deltaS/S=mu*deltat + sigma*randn*deltat^0.5 + nu*deltaq
%nu=average jump size - proportional to increase in asset price
%deltaq = 1 w/ probability lambda*deltat
%deltaq = 0 w/ probability 1-lambda*deltat
function SimJD ()
close all
clear all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
%set(0,'DefaultAxesColorOrder',[0 0 0;0 1 0;0 0 1; 1 0 0 ; 1 0 1; 0 1 1])
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

Szero=50;
mu=0.25; %drift
vol=0.1; %Volatility
musig2=mu-0.5*vol^2;
lambda=0.5; % rate of jumps per year = Intensity of Poisson Process
q1=-1.5; q2=0.5; %Max Jump Down, Up
%In constant jump case usually set nu to constant value relative to 
%current price, e.g., ~ -0.14 = downward jumps for market crash, default...
%For illustration  equate constant nu to numean of uniform jump process
%nu (1 Jump Size) = nuMean (distributed)
nuMean=(exp(q2)-exp(q1))/(q2-q1) - 1  
%nuMean = Average Jump size measured relative to previous stock price
logNuP1=log(nuMean+1) %Drift of ln(jumps)

steps=500; time=11;%years
dt=time/steps; SQRTdt=sqrt(dt);
Jump1Diff(1)=Szero;         Jump1(1)=0; %  Diffusion(1)=Szero;
JumpDistDiff(1)=Szero;      JumpDist(1)=0;
Jump1DiffExpect(1)=Szero;   JumpDistDiffExpect(1)=Szero;
    
for i=2:steps
Jump1DiffExpect(i)=Jump1DiffExpect(i-1)*exp((mu+nuMean*lambda)*dt);
JumpDistDiffExpect(i)=JumpDistDiffExpect(i-1)*...
    exp((mu+nuMean*lambda)*dt);
%we set nu(1Jump)=nuMean(distributed) -> E[1Jump]=E[Distributed Jump]
   if (lambda*dt>rand)
    Jump1Diff(i)=Jump1Diff(i-1)*...
        exp(musig2*dt+vol*randn*SQRTdt+logNuP1);
    Jump1(i)=Jump1Diff(i-1)*(exp(logNuP1)-1); %One-Size Jump
    Q=q1+(q2-q1)*rand;
    JumpDistDiff(i)= JumpDistDiff(i-1)*...
        exp(musig2*dt+vol*randn*SQRTdt+Q);
    JumpDist(i)=JumpDistDiff(i-1)*(exp(Q)-1); %Distributed-Size Jump
   else %No Jumps only Drift-Diffusion
    Jump1Diff(i)=Jump1Diff(i-1)*exp(musig2*dt+vol*randn*SQRTdt);
    JumpDistDiff(i)= JumpDistDiff(i-1)*exp(musig2*dt+vol*randn*SQRTdt);
    Jump1(i)=0; JumpDist(i)=0;
end
end

xax=(1:steps)*dt;
figure
subplot (2,1,1); plot(xax,Jump1DiffExpect,':',xax,Jump1Diff,xax,Jump1,'--')
legend('One-Size-J/D Expected','One-Size-J/D','One-Size-Jump ')
title('One-Size-Jump / Diffusion')
xlabel('Time [Years])'); ylabel('Price, S_t ')
axis tight

subplot (2,1,2); 
plot(xax,JumpDistDiffExpect,':',xax,JumpDistDiff,xax,JumpDist,'--')
legend('Distrbuted-J/D Expectected',...
    'Distrbuted-J/D','Distrbuted-Jump')
title('Distrbuted-Jump / Diffusion')
xlabel('Time [Years])'); ylabel('Price, S_t ')
axis tight
end
