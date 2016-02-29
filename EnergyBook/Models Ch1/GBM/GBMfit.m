function GBMfit (S)
%GBMfit().m Geometric Brownian Motion of Actual Data (S) 
%or Self-Simulation (no input)
%calculates drift mu and volatility from ln(St)-ln(St-1) data 
close all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

TimeDelta=1/252; 
%assume daily prices but could add as an input to function
if (nargin == 1), %Stock Price Data Input 
    %or read directly s=load('XOMprice.dat');
    steps = length(S);
    TimeLength=TimeDelta*steps; time = linspace(0,TimeLength,steps);
end
if (nargin == 0), %self-simulation 
    mu=0.1 %percentage drift
    Dsig =0.27 %annualized volatility 
    Szero=30; %Initial Price
    steps = 7.5*252; %assume 7.5 years
    TimeLength=TimeDelta*steps; 
    time = linspace(0,TimeLength,steps);  %years
    sigma = Dsig*sqrt(TimeDelta).*randn(1,steps); 
    %vector of random movements
    etadt=(mu-0.5*Dsig^2)*TimeDelta;
    S=zeros(1,steps); S(1)=Szero;
    for i = 2:steps
        S(i)=S(i-1)*exp(etadt+sigma(i));
    end
end
%logarithmic returns log(St/S0) is normally distributed 
LogDelta=log(S(2:end))-log(S(1:end-1)); %=log(s(i)/(s(i-1))
StanDev = std(LogDelta);
%variance vol^2*t 
EstVol=StanDev/sqrt(TimeDelta) %Estimated annualized volatility 
%log increments of GBM are normal
%mean(LogDelta)=(mu-0.5sig^2)t
EstMu=mean(LogDelta)/TimeDelta + 0.5*EstVol^2/2
Expected=zeros(1,steps); SD=zeros(1,steps);
Expected = S(1).*exp((EstMu).*(time));     %E[St]=S0exp(mu*t)
SD = Expected.*sqrt(exp(EstVol^2.*time) -1);
%SDV[St]=S0exp(mu*t)[exp(vol^2*t)-1]^0.5
time=time+2003; %start in Year 2003
figure
plot(time,Expected,'-',... 
    time,Expected+SD,'--',...
    time,S,':',...
    time,Expected-SD,'-.'), grid on
xlabel('Time '); ylabel('Price')%%
legend('Expected', '+1 Standard Deviation', 'Price Path',...
    '-1 Standard Deviation','location','NorthWest');
if (nargin == 1), %Stock Price Data Input 
    title('Geometric Brownian Motion fit to XOM');
end
if (nargin == 0), %self-simulated
    title('Geometric Brownian Motion of Simulated Price Process');
end
xlim ([2003 2011])
end
