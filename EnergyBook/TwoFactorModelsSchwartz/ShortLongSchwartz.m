function ShortLongSchwartz(CL)
%ShortLongSchwartz fits Futures Prices to Schwartz-Smith Model
%A single set of futures prices at one instant in time
%is used to fit to 2-factor model of
%equilibrium price level (xi) and short-deviations (chi)
%ln(St)=chi+xi
%A time-series is need for correct fit but here we are
%interested in generating plot to illustrate the difference
%in futures-valuation in risk-neutral world(Q) vs.
%expected spot price under real-world measure (P).
%The difference principally manifests as short-term and 
%long-term risk premiums. Original derivation 
%in Schwartz and Smith, Short Term Variations 
%and Long-Term Dynamics, Management Science, 46, 893 (2000)
%Below use light sweet crude futures data from Sep 2010 
%

close all
clear all
clc

set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);
if (nargin == 0), %Check for Futures Data Input 
   CL=load('CL.dat'); %read directly 
end %Format should have 2 columns: time , price
t=(CL(:,1)-CL(1,1))/365+1/12; %Time in Years% t(spot)=zero
Data=log(CL(:,2)); %Prices of Futures Contracts

%find slope and intercept of long dated Futures 
mid=round(length(t)/2);
HalfT=t(mid:end);
HalfData=Data(mid:end);
%Use second half of Futures Data
AveT=mean(HalfT);
AveData=mean(HalfData);
slope=sum((HalfT-AveT).*(HalfData-AveData))...
    /sum((HalfT-AveT).^2);
%Slope of long dated futures is uXI-lambdaXI-0.5sigmaXI^2
LongTermIntercept=AveData-slope*AveData;
LongTermExtrap=slope*t+LongTermIntercept; %Line for Plot
%set Starting parameters, i.e., initial guess for fminsearch
%To Start assume closest futures contract approx spot
%Extrapolation long dated future curves back to time zero
%Guess Deviation chiStart is approximately difference
chiStart=Data(1)-LongTermIntercept;
XIstart=LongTermIntercept;
kappaStart=1.5; %Approximate Values found in Schwartz (2000)
sigmaKstart=.0029;
sigmaXIstart=.0015;
rhoStart=0.3;
lambdaChiStart=0.157;
%Slope of long dated futures is uXI-lambdaXI-0.5sigmaXI^2
%ignore lambdaXI and SigmaXI since small
uXIStarStart=slope; %mu(xi)star = mu(xi)-lambda(xi)

Guess=[chiStart,XIstart, kappaStart , sigmaKstart ,...
    sigmaXIstart, rhoStart, lambdaChiStart, uXIStarStart];
options=[];%options=optimset('Display','iter');
Estimates=fminsearch(@myfit,Guess,options,t,Data);
chi=Estimates(1) 
XI=Estimates(2)
kappa=Estimates(3)
sigmaC=Estimates(4)
sigmaXI=Estimates(5)
rho=Estimates(6)
lambdaChi=Estimates(7)
uXIStar=Estimates(8)
%Recalculate Risk-Neutral Futures curve
%with Same Formula in minimized function 'myfit'
A=uXIStar*t-(1-exp(-kappa*t))*lambdaChi/kappa...
    +0.5*((1-exp(-2*kappa*t))*0.5*sigmaC^2/kappa)...
    +0.5*sigmaXI^2*t+...
    (1-exp(-kappa*t))*rho*sigmaC*sigmaXI/kappa;
Fit=exp(-kappa*t)*chi+XI+A;

% Need time series data to estimate long run risk premium
% ,e.g., with a Kalman filter and thus real world mu(xi)
% Just for illustration assume short run risk premium is 
% same as long run risk premium
lambdaXI=lambdaChi;
uXI=uXIStar+lambdaXI %mu(xi)Start=mu(xi)-lambda(xi)
%Slope of Expected Spot Price far into the future
EstLongSlope=uXI+0.5*sigmaXI^2;
%In Long Run, Expected Spot prices acts as if started from
%effective Long run price, which is equilibrium price 
%XI plus premium for short term volatility
EffectiveLongRunPrice=...
    (XI+0.25*sigmaC^2/kappa+rho*sigmaC*sigmaXI/kappa)
LongSpotPrice=EstLongSlope*t...
    +EffectiveLongRunPrice;
%line for graph
%Above Calculated Asymptotic (Straight) lines for
%Long-Dated Futures Contracts (LongTermExtrap)
%and Long-Dated Expected Spot Price (LongSpotPrice)
%Difference in Slope is 
%lambda(xi)= Long-Term Risk Premium
%Extrapolate Both lines to Time Zero
%Difference in Intercepts 
%=(LongTermIntercept-EffectiveLongRunPrice)
%is lambda(chi)/kappa = Short-Term Risk Premium

%Calculate Expected Spot Price
A2=uXI*t+0.5*((1-exp(-2*kappa*t))*0.5*sigmaC^2/kappa)...
    +0.5*sigmaXI^2*t+...
    (1-exp(-kappa*t))*rho*sigmaC*sigmaXI/kappa;
ExpectSpot=exp(-kappa*t)*chi+XI+A2;

figure
tY=CL(:,1); %Calender time % Simplify notation
tNow=CL(1,1)-30;
plot(tY,Data,'*',tY,Fit,'black',tY,LongTermExtrap,'g',...
    tY,LongSpotPrice,tY,ExpectSpot,'.-',...
    tNow,XI,'x',tNow,XI+chi,'+'); 
datetick('x'); 
xlabel('Time'); ylabel('Ln(Price)')
title('Light Crude Futures: Schwartz and Smith 2-Factor')
%
legend('Futures Contracts', 'Risk-Neutral Fit',...
    'Futures Curve','Equilbrium Price', ...
    'Expected Spot','Location','SouthEast')

text(tNow,XI,'\xi ','HorizontalAlignment','right')
txstr3(1) = {' \xi+\chi'}; 
txstr3(2) = {'   Est. Spot'}; 
text(tNow,XI+chi,txstr3,...
    'VerticalAlignment','top','HorizontalAlignment','center')

text(tY(end-3),ExpectSpot(end-3),...
    'Slope=\mu_\xi+\sigma_\xi^2/2 ',...
    'VerticalAlignment','bottom','HorizontalAlignment','right')
text(tY(end-10),Data(end-10),...
    'Slope=\mu_\xi-\lambda_\xi+\sigma_\xi^2/2',...
    'VerticalAlignment','top','HorizontalAlignment','left')

txstr2(1) = {'\lambda_\xi Long Term'}; 
txstr2(2) = {'Risk Prem.'}; 
text(tY(end-4),ExpectSpot(end-6),txstr2)

txstr(1) = {'Short-Term'}; txstr(2) = {'Risk'}; 
txstr(3) = {'Premium'}; txstr(4) = {'\lambda_\chi/\kappa'};
text(tNow,LongTermExtrap(60),txstr)

text(tY(1),LongTermExtrap(6),'\downarrow')
text(tY(1),LongSpotPrice(1),'\uparrow')

eXI=exp(XI)
eXIplusChi=exp(XI+chi)

end

%%%Squared Difference of Experimental Futures Prices vs.
%%%Schwartz Smith 2-factor model calculated from params
function sse=myfit(params,T,lnF)
chi=params(1);
XI=params(2);
kappa=params(3);
sigmaC=params(4);
sigmaXI=params(5);
rho=params(6);
lambdaChi=params(7);
uXIStar=params(8);

A=uXIStar*T-(1-exp(-kappa*T))*lambdaChi/kappa...
    +0.5*((1-exp(-2*kappa*T))*0.5*sigmaC^2/kappa)...
    +0.5*sigmaXI^2*T+...
    (1-exp(-kappa*T))*rho*sigmaC*sigmaXI/kappa;
EstFit=exp(-kappa*T)*chi+XI+A;

Error_Vector=EstFit - lnF;
% % minimize the sum of squares error
sse=sum(Error_Vector.^2);
% % sse=Error_Vector(:)'*Error_Vector(:); 
end
