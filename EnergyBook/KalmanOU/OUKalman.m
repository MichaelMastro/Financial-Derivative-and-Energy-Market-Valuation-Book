% OUKalman simulats and then filters a mean reversion process. 
% Underlying parameters (mu, lambda, and sigma) of this OU process 
% are estimated by the matlab minization function fminsearch. 
% The beahvior of mu, lamda, and 
% sigma are tested by OUloglikelihoodfn, which returns the 
% negative of the maximum likelihood to fminsearch.
% SimulateOUkalman generates mean reverting data as well as a 
% corrsponding price noisey observation. This step can be replaced
% by importing real spot data (Obs) and the hidden spot (S) will be 
% intrinsically estimated by the Kalman filter
function OUKalman()
close all
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

% global avoids passing dataset 
global S; global Obs; 
global TimeDelta; %keeps fminsearch function call clearer
global xVec; % predicted price path

lambda=3 %mean reversion rate
mu=1 %long-term mean
Dsig =0.5 %Driving Sigma 
Szero=3; %Initial Price
stepsInput = 201; TimeLength=5; 
S=zeros(1,stepsInput); Obs=zeros(1,stepsInput);
time = linspace(0,TimeLength,stepsInput);  TimeDelta= time(2)-time(1);

SimulateOUkalman (mu, lambda,  Dsig, Szero, stepsInput); 
%generate S and Obs
%Can replace this step by importing observed spot data (obs)

%Supply initial guess values
    muGuess=9; 
	lambdaGuess=2;
    sigGuess=0.9;
    RGuess=0.01; 
param=zeros(3,1);
	param(1)=muGuess;
	param(2)=lambdaGuess;
    param(3)=sigGuess;
    param(4)=RGuess; 
pnew=fminsearch('OUloglikelihoodfn',param)
muNew = pnew(1);
lambdaNew = pnew(2);
DsigNew = pnew(3);
RNew = pnew(4);
    
fprintf(1, '\t\t\t\t mu \t lambda \t sigma \n');
fprintf(1, 'True  \t\t %6.2f \t %6.2f \t %6.2f \n',...
    mu, lambda, Dsig); 
fprintf(1, 'Estimated \t %6.2f \t %6.2f \t %6.2f \n',...
    muNew, lambdaNew, DsigNew); 
fprintf(1, 'Measurement Noise Covariance = R = %6.4f \n', RNew); 
figure; plot (time, S, time, Obs,'--', time, xVec, '-.');
legend ('Real Price', 'Observed Price','Predicted Price');
xlabel('Time'); ylabel('Asset Price')
end