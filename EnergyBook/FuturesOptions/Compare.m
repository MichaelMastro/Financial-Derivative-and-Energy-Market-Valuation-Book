function  compare ( );
% Black's model for option on futures contract
% forward price at maturity of option is log-normal distributed
% Generate various graphs to show dependence
% on T,K,F,vol as well as Greeks
close all;
clc;
set(0,'defaultaxeslinewidth',3); set(0,'defaultlinelinewidth',3);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

r=0.1; vol=0.5; K=50;
% #1 Examine Approach to c=Max(0,e^{-rT}(F-K)
% as a function of (F) for three Time to Expiration
F=20:0.1:80;
Expiration=[1/12,1,2];
for i=1:3
    T=Expiration(i);
    c(i,:) = BlackCall (K,F,T,vol,r );
end

figure
% ertFK=max(0,exp(-r.*T)*(F-K));
% FertK=max(0,(F-exp(-r.*T)*K));
FK=max(0,(F-K));

plot (F,c,F,FK,'--',K,0,'*');
xlabel('Futures Price'); ylabel('Call Value');
title('Black Model c=e^{-rT}(F\Phi(d_1)-K\Phi(d_2)');
legend ('1 Month to Expiration','1 Year to Expiration',...
    '2 Years to Expiration',...
    'F-K','Location','NorthWest');%'e^{-rT}(F-K)',
text(K-1,-1,'K')
text(K+10,10,'(F-K)')
txstr1(1) = {'Increasing'}; 
txstr1(2) = {'Time to'};
txstr1(3) = {'Expiration \uparrow'}; 
text(30,5,txstr1);

%%%%%%%%%%%%Second Graph
% Show impact of cost of carry b=r-d-y+u and 
% relation to dividend d, convenience yield y, 
% and storage u
figure
S=20:0.1:80; 
T=1;

% ertSK=max(0,exp(-r.*T)*(S-K));
% SertK=max(0,(S-exp(-r.*T)*K));

%%% %never exercise call early b<=r
% for example non-dividend paying stock 
% when d=0 -> b=r since b=r-d
d=0;%%% d=0 b=r-d=r
F=S.*exp((r-d).*T);
c = BlackCall (K,F,T,vol,r );
% Black model is interchangeable with Black-Scholes
SmX=max(0,exp(-d.*T).*S-exp(-r.*T).*K); 

SK=max(0,S-K);

subplot(3,1,1)
plot (S,c,S,SmX,'--',S,SK,'-.',K,0,'*'); %S,ertSK,
title('B-S Non-Dividend Stock d=0, b=r');
legend ('S\Phi(d_1)-Ke^{-rT}\Phi(d_2)',...
    'c_{min}=S-e^{-rT}K',...
        '(S-K)','Location','NorthWest');
txstr2(1) = {'       \downarrow '}; 
txstr2(2) = {'    c_{min}>S-K'}; 
txstr2(3) = {'Never Exercise'};
txstr2(4) = {'Early (b>=r)-> C=c'}; 
text(60,10,txstr2);

%%% d=0.75*r
%%%% potentially exercise early when d>0 -> b<r since b=r-d
d=0.75*r; 
F=S.*exp((r-d).*T);
c = BlackCall (K,F,T,vol,r );
SmX=max(0,exp(-d.*T).*S-exp(-r.*T).*K);

subplot(3,1,2)
plot (S,c,S,SmX,'--',S,SK,'-.',K,0,'*'); %S,ertSK,
title('B-S: Stock d=0.75r or Commodity y=0.75r');
xlabel('Asset Price'); 
legend ('e^{-(d=0.75r)T}S\Phi(d_1)-Ke^{-rT}\Phi(d_2)',...
    'c_{min}=e^{-(d=0.75r)T}S-e^{-rT}K',...
 '(S-K)', 'Location','NorthWest');
ylabel('Value');

text(S(end-20),SK(end-20),'\uparrow '); 
txstr3(1) = {'      c_{min}<S-K'}; 
txstr3(2) = {'(b<r)->Possible'};
txstr3(3) = {'Early Exercise'}; 
text(65,10,txstr3);

d=r;%%%d=r
F=S.*exp((r-d).*T);
c = BlackCall (K,F,T,vol,r );
SmX=max(0,exp(-d.*T).*S-exp(-r.*T).*K);
subplot(3,1,3)
plot (S,c,S,SmX,'--',S,SK,'-.',K,0,'*'); %S,ertSK
title('Black Futures b=0 or B-S Stock d=r');

legend...
('e^{-(d=r)T}S\Phi(d_1)-Ke^{-rT}\Phi(d_2)=e^{-rT}(F\Phi(d_1)-K\Phi(d_2))',...
    'c_{min}=e^{-(d=r)T}S-e^{-rT}K=e^{-rT}(F-K)',...
        '(S-K)', 'Location','NorthWest');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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