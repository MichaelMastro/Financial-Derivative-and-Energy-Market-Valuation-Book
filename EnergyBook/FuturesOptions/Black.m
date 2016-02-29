function  Black ( );
%Black's model for option on futures contract
%forward price at maturity of option is log-normal distributed
%Generate various graphs to show dependence
%on T,K,F,vol as well as Greeks
close all;
clc;
set(0,'defaultaxeslinewidth',3); set(0,'defaultlinelinewidth',3);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);
set(0,'DefaultAxesLineStyleOrder',{'-',':','--','.-'})

r=0.05;vol=0.3;K=50;
%#1 Examine Approach to c=Max(0,e^{-rT}(F-K)
%as a function of (F) for three Time to Expiration
F=1:0.1:80; 
for i=1:3
    T=i-1;
    c(i,:) = BlackCall (K,F,T,vol,r );
end
figure
plot (F,c,K,0,'*');
xlabel('Futures Price'); ylabel('Call Value');
title('Black Model');
legend ('At Expiration','1 Year to Expiration',...
    '2 Years to Expiration','Location','NorthWest');
text(K-1,-1,'K')
text(K+10,10,'e^{-rT}(F-K)')
txstr1(1) = {'Increasing'}; 
txstr1(2) = {'Time to'};
txstr1(3) = {'Expiration \uparrow'}; 
text(20,5,txstr1);

%#2a Examine Delta as a function of futures price
%Greek Delta=dc/dF
clear F T c;
F=1:0.1:80;
for i=1:3
    T=i-1+1/365;
    DeltavF(i,:)=BlackDelta (K,F,T,vol,r );
end
figure
subplot (2,1,1); plot (F,DeltavF,K,0,'*');
xlabel('Futures Price');
ylabel('\Delta = \partialc/\partialF');
title('Delta Dependence in Black Model');
legend ('1 Day to Expiration','1 Year to Expiration',...
    '2 years to Expiration','Location','NorthWest');
text(K-1,-0.05,'K')
txstr1(1) = {'Increasing'}; 
txstr1(2) = {'Time to'};
txstr1(3) = {'Expiration \uparrow'}; 
text(20,0.25,txstr1);

%#2b Examine Delta as a function of time to expiration
clear F T c;
T=0:.05:5;
for i=1:3
    F=25*i;
    deltaSeries(i,:)= BlackDelta (K,F,T,vol,r );
end
subplot (2,1,2); plot (T,deltaSeries);
legend ('Out of the Money (F<K)','At the Money (F=K)',...
    'In the Money (F>K)','Location','North');
xlabel('Time to Expiration (\tau=T-t)'); 
ylabel('\Delta = \partialc/\partialF');

%#3 Examine Approach to time dependence of call value
%and the Greek theta=dc/d(T-t)
%as approaching expiration
clear F T c;
T=0:.05:5;
for i=1:3
    F=25*i;
    c(i,:) = BlackCall (K,F,T,vol,r);
    thetaSeries(i,:)= BlackTheta (K,F,T,vol,r );
end
figure
subplot (2,1,1);plot (T,c);
xlabel('\tau=T-t'); ; ylabel('Call Value');
title('Black Model');
legend ('Out of the Money (F<K)','At the Money (F=K)',...
    'In the Money (F>K)','Location','NorthEast');
subplot (2,1,2);plot (T,thetaSeries);
xlabel('Time to Expiration (\tau=T-t)'); 
ylabel('\theta = \partialc/\partial\tau');
title('Theta Dependence in Black Model');

%#5 Examine European put-call parity
figure
clear T F c;
F=1:0.1:80; 
T=1;
c = BlackCall (K,F,T,vol,r );
p = BlackPut (K,F,T,vol,r);
%FK=exp(-r*T)*(F-K)
plot (F,c,F,p,F,c-p,'--', F, p-c,'-.', K,0,'*');
xlabel('Futures Price'); ylabel('Value');
ylim([0 30]); xlim([20 80]);
title('European Put-Call Parity');
legend ('Call','Put',...
    'c-p=e^{-rT}(F-K)','p-c=e^{-rT}(K-F)',...
    'Location','North');
end

function c = BlackCall (K,F,T,vol,r);
    d1 = ( log(F./K)+ (vol^2/2).*T)./ (vol.*sqrt(T));
    d2 = d1-vol*sqrt(T);
    c = exp(-r*T).*(F.*myNormCDF(d1)-K.*myNormCDF(d2));
end

function p = BlackPut (K,F,T,vol,r);
    d1 = ( log(F./K)+ (vol^2/2).*T)./ (vol.*sqrt(T));
    d2 = d1-vol*sqrt(T);
    p = exp(-r*T).*(K.*myNormCDF(-d2)-F.*myNormCDF(-d1));
end

function theta = BlackTheta (K,F,T,vol,r);
    d1 = ( log(F./K)+ (vol^2/2).*T)./ (vol.*sqrt(T));
    d2 = d1-vol*sqrt(T);
    theta = -exp(-r*T).*F.*(K.*myNormPDF(d1)).*vol...
        -r.*K.*exp(-r*T).*myNormCDF(d2)...
        +r.*F.*exp(-r*T).*myNormCDF(d1);
end

function theta = BlackDelta (K,F,T,vol,r);
    d1 = ( log(F./K)+ (vol^2/2).*T)./ (vol.*sqrt(T));
    theta = exp(-r*T).*myNormCDF(d1);
end