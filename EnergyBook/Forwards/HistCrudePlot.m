% HistCrudePlot displays the WTI spot and 1, 2, 3, and 4 month
% futures data over the past approximate twenty years
% Also plots two characteristic forwards curves for a
% normal futures market and an inverted market
close all
clear all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold');    
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);


[d,S,f1,f2,f3, f4] = textread('CrudeOilHistT.txt','%f %f %f %f %f %f');
P=[S,f1,f2,f3, f4];
L=ones(1,length (d)); 
L0=0*L; L1=1*L; L2=2*L; L3=3*L;L4=4*L;

subplot(2,1,1)
plot3 (d,L0,S,d,L1,f1,d,L2,f2,d,L3,f3,d,L4,f4)

%surfc(P);
title('WTI Historical Spot and Futures');

ylabel('Spot/Futures Month');
xlabel('Date [Year]');
zlabel('Price US Dollars');
axis tight

SF4=S-f4;
SF2=S-f4;

[MaxSF,MaxIndex]=max(SF4);
[MinSf,MinIndex]=min(SF4);
fprintf (1, 'Severe Inversion in Year %4.0f Month %4.0f \n',...
    d(MaxIndex), 12*(-d(MaxIndex)+round(d(MaxIndex))) )
fprintf (1, 'Severe Normal in Year %4.0f Month %4.0f \n',...
    d(MinIndex), 12*(-d(MinIndex)+round(d(MinIndex))) )


subplot (2,1,2)
x=0:4;
plot(x,P(MinIndex,:),'-',x,P(MaxIndex,:),'--')
xlabel('Spot/Futures Month');
ylabel('Price US Dollars');
MaxDate=num2str(round(d(MaxIndex)));
MinDate=num2str (round(d(MinIndex)));
MinDate=['Normal ',MinDate];
MaxDate=['Inverted ',MaxDate];
legend (MinDate, MaxDate,'location','East')

