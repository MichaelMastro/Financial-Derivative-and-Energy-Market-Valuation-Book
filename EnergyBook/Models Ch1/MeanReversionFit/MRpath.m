function MRpath ()
%MRpath simulates Path of Mean Reversion Process
%then uses external function calls for unweighted 
%and weighted least squares fit as well as unweighted
%and weighted maximum likelihood fits. 
%The rapidly reverting data has less data points 
%so one can emphasize certain data to
%in theory provide a more accurate fit to lambda

%
close all;
clear all;
clc;
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesColorOrder',[0 0 0;0 1 0;0 0 1; 1 0 0 ; 1 0 1; 0 1 1])
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

steps = 201;
TimeLength=5; %years
lambda=3 %mean reversion rate
mu=2 %long-term mean
Dsig =1.5 

Szero=5; %Initial Price

S=zeros(1,steps); Expected=zeros(1,steps); SD=zeros(1,steps);

time = linspace(0,TimeLength,steps);
TimeDelta = time(2)-time(1)

Expected = Szero.*exp(-lambda.*time)+mu.*(1-exp(-lambda*time));
%SD is constant but calculated as vector 
SD=Dsig.*sqrt( (1-exp(-2.*lambda.*TimeDelta*ones(1,steps)))./(2.*lambda));
weightSD =SD.*linspace(0.1,1.9,steps);
%Place more weight (importance) into earlier points that are rapidly
%changing to hopefully better fit the mean reversion rate lambda
S= Expected + SD.*randn(1,steps);

x=S(1:end-1); y=S(2:end);

figure
plot(time,Expected,'-',... 
    time,Expected+SD,'--',...
    time,S,':',...
    time,Expected-SD,'-.'), grid on
xlabel('Time '); ylabel('Price')%%
legend('Expected', '+1 Standard Deviation',...
    'Price Path`','-1 Standard Deviation');
title(['Mean Reversion Process']);

TrueSlope = exp(-lambda*TimeDelta);
TrueIntercept = mu*(1-TrueSlope);
trueSD= Dsig.*sqrt( (1-exp(-2.*lambda.*time(steps))) ./ (2.*lambda)); 
TrueLine = polyval([TrueSlope TrueIntercept],x);

fprintf(1, '\t\t\t\t mu \t lambda \t sigma \t\t Slope');
fprintf(1, '\t\t Intercept \t Standard Deviation \n');
fprintf(1,...
'True  \t\t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \n',...
    mu, lambda, Dsig, TrueSlope, TrueIntercept, trueSD); 

%%%Calculate by Unweighted Least Squares
[LSMu, LSSigma, LSLambda, LSslope, LSintercept, LSstanDev]...
   = WeightedLeastSquaresOU (S,TimeDelta);
%LS fitting assumes equal weight
fprintf(1,...
'Standard LS\t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \n',...
   LSMu, LSSigma, LSLambda, LSslope, LSintercept, LSstanDev)
UnweightedLSline = polyval([LSslope, LSintercept],x);

%%%Calculate by Weighted Least Squares
[wLSMu, wLSSigma, wLSLambda, wLSslope, wLSintercept, wLSstanDev]...
   = WeightedLeastSquaresOU (S,TimeDelta, weightSD(2:end));
fprintf(1,...
'Weight LS\t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \n',...
   wLSMu, wLSSigma, wLSLambda, wLSslope, wLSintercept, wLSstanDev);  
WeightedlsLine = polyval([wLSslope, wLSintercept],x);

%%%Calculate by Unweighted Maximum likelihood
[MLmu, MLlambda, MLsigma] = weightedML(S,TimeDelta);
fprintf(1, 'Standard ML\t %6.2f \t %6.2f \t %6.2f \n',...
    MLmu, MLlambda, MLsigma)

%%%Calculate by Weighted Maximum likelihood
[wMLmu, wMLlambda, wMLsigma] = weightedML(S,TimeDelta, weightSD(2:end));
fprintf(1, 'Weight ML\t %6.2f  \t %6.2f \t %6.2f \n',...
    wMLmu, wMLlambda, wMLsigma)

figure
plot(x,y,'o',x,TrueLine,'-',x,UnweightedLSline,'--',x,WeightedlsLine,'.-')
xlabel('Previous Price, S_t_-_1 '); ylabel('Price, S_t ')%%
legend('Data','True Line','Unweigted LS Fit',...
    'Weighted LS Fit', 'location', 'NorthWest');

%Or use one of several Matlab fitting function 
%Matlabfit = polyfit(x,y,1);
%fprintf(1, 'Matlab calc. %6.2f Slope \t %6.2f Intercept \n', Matlabfit);
%MatlabfitLine = polyval(Matlabfit,x);
end
