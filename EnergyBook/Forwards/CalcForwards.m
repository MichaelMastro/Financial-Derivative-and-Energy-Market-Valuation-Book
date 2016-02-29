% CalcForwards simulates Geometric Brownian Motion of 
% spot price and Forward price
close all
clear all
clc

set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold');    
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);


TimeDelta=1/252; %assume daily prices but could add as an input to function
mu=0.1 %percentage drift
Dsig =0.1 %annualized volatility 
Szero=100; %Initial Price
steps = 7.5*252; %assume 7.5 years
TimeLength=TimeDelta*steps; time = linspace(0,TimeLength,steps);  %years
randn('state', 4);
sigma = Dsig*sqrt(TimeDelta).*randn(1,steps); %vector of random movements
alphadt=(mu-0.5*Dsig^2)*TimeDelta;
S=zeros(1,steps); S(1)=Szero;
for i = 2:steps
    S(i)=S(i-1)*exp(alphadt+sigma(i));
end

%logarithmic returns log(St/S0) is normally distributed 
LogDelta=log(S(2:end))-log(S(1:end-1)); %=log(s(i)/(s(i-1))
StanDev = std(LogDelta);
%variance sigma^2*t 
EstVol=StanDev/sqrt(TimeDelta) %Estimated annualized volatility 
%log increments of GBM are normal relative to price t
EstMu=mean(LogDelta)/TimeDelta + EstVol^2/2

Expected=zeros(1,steps); SD=zeros(1,steps);
Expected = S(1).*exp((EstMu).*(time));     %E[St]  =S0exp(mu*t)
time2003=time+2003; %start in Year 2003
ExpectedLine=Expected(end)*ones(1,steps);

r=0.03; %risk-free rate;
u=0.2; %storage cost
y1=0.05; y2=0.2; y3=0.35; %convenience yield
tau=time(end)-time;
Ftau1=S.*exp(tau.*(r+u-y1));%
Ftau2=S.*exp(tau.*(r+u-y2));%
Ftau3=S.*exp(tau.*(r+u-y3));%

figure
hl1 = line(time2003,S,'LineStyle',':','Color','b');
hl2 = line(time2003,Expected,'LineStyle','-.','Color','g');
hl2 = line(time2003,ExpectedLine,'LineStyle','-','Color','k');
ax1 = gca;
set(ax1,'XColor','b','YColor','k')
xmin1 = time2003(1);xmax1=time2003(end); % xlimits et 
%set(ax1,'xlim',xlimit)
ylimit = [0,max([Ftau1 S])*1.8]; % ylimits
%set(ax1,'ylim',ylimit)
axis([xmin1,xmax1,ylimit]);
xlabel('Time '); ylabel('Price')%%
legend('Spot Price', 'E[S_t]','E[S_T]','Location','NorthWest');
%title('Geometric Brownian Motion Spot and Future Price Process'); 
title('      \tau=T-t');

ax3 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');

hl3 = line(tau, Ftau1,'LineStyle','-','Color',[0.7 0 0],'Parent',ax3);
hl4 = line(tau, Ftau2,'LineStyle','-','Color',[0.5 0 0],'Parent',ax3);
hl5 = line(tau, Ftau3,'LineStyle','--','Color',[0.3 0 0],'Parent',ax3);
xmin=tau(1);xmax=tau(end);
axis([xmax,xmin,ylimit]);
set(ax3,'xdir','reverse')
legend('F_{\tau} y=0.05','F_{\tau} y=0.2','F_{\tau} y=0.35',...
    'Location','NorthEast')

str1(1) = {'Contango'};
str1(2) = {'F_{\tau}>E[S_T]'};
text(tau(end*.7),Ftau1(round(end*0.1)),str1,'HorizontalAlignment','right')
str2(1) = {'Normal Backwardation'};
str2(2) = {'F_{\tau}<E[S_T]'};
text(tau(end*.9),min(Ftau3),str2,'HorizontalAlignment','right')