% CYKalman fits 2-factor model parameters
% to a set of (1,2,3 and 4 month) futures contracts. 
% The hidden (factor) spot price and convenience yield 
% is extracted over the entire time series

% Innovation Covariance matrix for this particular model and data set 
% is close to singular, which creates ill-conditioned problem 
% for inversion of Innovation Cov.: inv(V) = inv(H*P*H' + R)  
% Four distinct techniques used to improve problem:
% 0:'Normal': Kalman filter in Vector form but
%  More Symmetric Joseph Covariance Correction
% 1:'scalar' Normal Kalman filter, Optimal Kalman Gain;
%   measurements are uncorrelated (=R is diagonal matrix)
%   -> Sequential Scalar Observation Update is used,
%   -> faster since Innovation Covariance is a scalar 
%   so do not need (slow) matrix inversion
% 2:'Schmidt-SR': error covariance matrix Pt|t and Pt|t-1 propagate  
%  in square root form -> maintain positive semi-definite nature of
%  error covariance: Kaminski, Bryson, Schmidt (1971) IEEE Trans.
%  Automatic Control 16, 727
% 3: 'Carraro-SR': Squaring / Square-Root Filter: error covariance   
%  matrix is not propagated in square root form.   
%  squaring the square root factors -> roundoff errors; 
%  however, this algorithm is simpler
%  and for this application, the covariance is needed to 
%  calculate the likelihood
%  Carraro (1988) Comp. Sci in Econ. and Management 1, 41

function CYKalman(flaginput)
% flag='Normal-Joseph';%'scalar';%'Schmidt-SR';%'Carraro-S/SR'
% flaginput 0 -> flag='Normal-Joseph' -> normal vector-
% based Kalman filter
% flaginput 1 -> flag='scalar' -> measurement cov mat
% R only has variance terms along its diagonal -> use
% division instead of matrix inversion
% flaginput 2 -> flag='Schmidt-SR' -> square root filter
% flaginput 3 -> flag='Carraro-S/SR'-> Squaring / square root 

tic
close all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold');    
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

global Obs CY   % avoid passing large arrays
global TimeDelta    %but this simplifies the code
global TimeMat MeasErr LnS %estimated Ln(Spot Price)
global flag

if (nargin==0)
    flaginput=0;
end

if (flaginput==1)
    flag='scalar'
elseif (flaginput==2) 
    flag='Schmidt-SR'
elseif (flaginput==3) 
    flag='Carraro-S/SR'
else
    flag='Normal-Joseph'
end

TimeMat=[1,2,3,4]'*(1/12); % Four series of futures contracts 
% with expiration of 1,2,3, and 4 months
TimeDelta=(1/12);

[d,S,f1,f2,f3, f4]=textread('CrudeOilHistT.txt','%f %f %f %f %f %f');
Obs=log([f1,f2,f3,f4]); 
samples=length(Obs(1,:));
LnS=ones(samples,1);
    RvarGuess=1e-7;
    AlphaGuess=1; kappaGuess=0.14; 
    sig1Guess=0.34; lambdaGuess=1.7;
    sig2Guess=0.43; rhoGuess=0.74;
    muGuess=0.58;
param=zeros(7,1);
	param(1)=AlphaGuess; param(2)=kappaGuess;  
    param(3)=sig1Guess; param(4)=lambdaGuess;
    param(5)=sig2Guess; param(6)=rhoGuess;
    param(7)=muGuess;
    param(8:11)=RvarGuess;    % diagonal variance terms
  
      
%options=optimset('TolX',1e-1,'TolFun',1e-1,'Display','iter');%
options=optimset('Display','iter','MaxIter',200);
%Set Tolerance to prevent unfeasible parameters
pnew=fminsearch('CYloglikelihoodfn',param,options);
   
fprintf(1, 'Estimation Type: ');
disp(flag);
fprintf(1, 'alphaStar \t kappa \t sigma1 \t lambda \t');
fprintf(1, 'sigma2 \t rho \t mu  \n');
fprintf(1, '%6.2f \t %6.2f \t %6.2f \t %6.2f \t', pnew(1:4)); 
fprintf(1, '%6.2f \t %6.2f \t %6.2f \n', pnew(5:7)); 
fprintf(1, 'Measurement Error Variance Terms  \n');
fprintf(1, '%6.2e\t%6.2e\t%6.2e\t%6.2e\n', pnew(8:11));
figure
Time=1988+[0:length(LnS)-1]/12;
subplot (2,1,1)
plot(Time,LnS,'-',Time,CY,'--')
title(flag)
xlabel('Time [Years]'); ylabel('ln (Dollars)'); axis tight
legend('Estimated Ln Spot','Convenience Yield','location','NorthWest')
subplot (2,1,2)
plot(Time,MeasErr)
xlabel('Time [Years]'); ylabel('Measurement Error'); axis tight
legend('1 Month Futures','2 Month Futures',...
    '3 Month Futures','4 Month Futures','location','NorthWest')
toc
end
%-------------------------------------------------------------