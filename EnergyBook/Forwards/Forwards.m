% Forwards calculates risk-free forward prices without storage cost and
% with Per Unit or Proportional Storage Cost
close all;
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold');    
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

S0=50; %Initial Stock Price
T=[0:.1:2]; %Time [Years]
r=0.05;         %Risk-Free Rate
U=1;            %Per Unit Cost of Storage
eu=0.02;        %Proportional Cost of Storage
F=(S0)*exp(r*T);        %Forward Price No Storage Cost
FU=(S0+U)*exp(r*T);     %Forward Price with Per Unit Cost of Storage
Feu=(S0)*exp((r+eu)*T); %Forward Price with Proportional Cost of Storage

subplot (2,1,1); plot(T,F); grid on
ylabel('Fair Value');
legend('S_0e^r^t','Location','NorthWest');
title('Risk-Free Forward Price');
midT=round(length(T)/2)+1; %Index of 1 Year Forward Price
text(T(midT),F(midT)+1,'F_0>(S_0)e^r^t\uparrow',...
    'HorizontalAlignment','right')
text(T(midT),F(midT)-1,'F_0<(S_0)e^r^t\downarrow',... 
     'HorizontalAlignment','right','VerticalAlignment','top')

subplot (2,1,2); plot(T,FU,'-',T,Feu,'--'); grid on
xlabel('Time '); ylabel('Value + Storage Cost');%%
legend('(S_0+U)e^r^t','S_0e^(^r^+^u^)^t','Location','NorthWest');

