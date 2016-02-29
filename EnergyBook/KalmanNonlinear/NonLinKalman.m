function  [] = NonLinKalman()
% NonLinKalman invokes Non-Linear (EKFbs) extended Kalman 
% Filter or (UTbs) unscented transform
% to calculate hidden stochastic volatility based on 
% observation of stock index price and at-the-money 
% 1-month duration call option price based on 
% Black-Scholes equation

close all
clc

set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

code='UTbs' %'EKFbs';%'UTbs' %%'EKFbs' or UTbs
 
% Global avoids passing dataset 
global obsS obsC obsRF TimeDelta 
% Retrieve Data from Filter for Plots
global xVec hxVec MRvec volSvec

obsS=load('PriceData.dat');
ObsSig=load('VolData.dat');
obsRF=load('RiskFree.dat');
obsC=load('CallData.dat');
len=length(obsS)

TimeDelta=1/12; % 1 Month Forward Option
Time=((1:len)/12)+1990; % Monthly Observations

%Supply initial guess values
    MuGuess=0.008; 
	RGuess=1e-3; 
    SigGuess=0.15;
    BetaGuess=0.93;
    AlphaGuess=0.02;

param=zeros(5,1);
	param(1)=MuGuess;
    param(2)=RGuess; 
    param(3)=SigGuess;
    param(4)=BetaGuess;
    param(5)=AlphaGuess; 
    options=optimset('MaxIter',20,'Display','iter');%
   
pnew=fminsearch(code,param,options)
MuNew = pnew(1);
RNew = pnew(2);
SigNew = pnew(3);
BetaNew = pnew(4);
AlphaNew = pnew(5);

fprintf(1, '\t\t\t\t mu \t Beta \t Alpha \t sigma \n');
fprintf(1,'Estimated \t %6.2f \t %6.2f \t %6.2f \t %6.2f\n',...
    MuNew, BetaNew, AlphaNew, SigNew); 
fprintf(1, 'Measurement Noise Covariance = R = %2.2e \n', RNew); 

figure; 
subplot (3,1,1)
plot (Time, obsS);
title(code)
legend ('Observed',...
    'Location','NorthWest' ); 
xlabel('Time'); ylabel('Index Price'); axis tight;
subplot (3,1,2)
plot (Time, MRvec,'--', Time, hxVec, '-.', Time, obsC);
legend ('Difference','Predicted', 'Observed',...
    'Location','NorthWest' ); 
xlabel('Time'); ylabel('Call Price'); axis tight;
subplot (3,1,3)
plot ( Time, ObsSig,Time, xVec, '-.', Time, volSvec,'--'); 
legend ( 'BS Imp. Vol.','Pred. Vol.',...
    'Stock Vol.','Location','NorthWest' ); axis tight;
xlabel('Time'); ylabel('Volatility'); 
end