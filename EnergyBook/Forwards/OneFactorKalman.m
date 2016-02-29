% OneFactorKalman fits the mean drift alpha, reversion rate kappa,
% sigma, and measurement error matrix R
% to a set of (1,2,3 and 4 month) futures contracts. 
% The hidden (factor) spot price is extracted 
% over the entire time series
function OneFactorKalman()

close all
clear all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold');    
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

global Obs          %global is poor programming technique 
global TimeDelta    %but this simplifies the code
global TimeMat MeasErr LnS %estimated Ln(Spot Price)

TimeMat=[1,2,3,4]*(1/12); %Four series of futures contracts 
%with expiration of 1,2,3, and 4 months are readily available
%from www.eia.doe.gov; a much larger maturation range is 
%used in practice and is more stable for model estimation
TimeDelta=(1/12);

[d,S,f1,f2,f3, f4]=textread('CrudeOilHistT.txt','%f %f %f %f %f %f');
Obs=log([f1,f2,f3,f4]); 
samples=length(Obs(1,:));
LnS=ones(samples,1);
    RvarGuess=0.001;
    AlphaGuess=11.6; AlphaStarGuess=2.43;
	kappaGuess=0.03; sigGuess=0.16;
param=zeros(4,1);
	param(1)=AlphaGuess; param(2)=kappaGuess;  
    param(3)=sigGuess; param(4)=AlphaStarGuess; 
    param(5:8)=RvarGuess;

options=optimset('Display','iter','MaxIter',100);%options=[];%
pnew=fminsearch('OneFactorloglikelihoodfn',param,options);
fprintf(1, '\t\t alpha(P) \t kappa \t sigma \t alpha(Q) \n');
fprintf(1, 'Estimated  %6.2f \t %6.3f\t %6.3f \t %6.2f \n', pnew(1:4)); 
fprintf(1, 'Measurement Error Variance Terms  %6.4f \t %6.4f \t %6.4f \t %6.4f\n', pnew(5:8));
figure
Time=1988+[0:length(LnS)-1]/12;
subplot (2,1,1)
plot(Time,exp(LnS),Time,exp(Obs))
xlabel('Time [Years]'); ylabel('Price US Dollars'); axis tight
legend('Estimated Spot','1 Month Futures','2 Month Futures',...
    '3 Month Futures','4 Month Futures','location','NorthWest')
subplot (2,1,2)
plot(Time,MeasErr)
xlabel('Time [Years]'); ylabel('Measurement Error'); axis tight
legend('1 Month Futures','2 Month Futures',...
    '3 Month Futures','4 Month Futures','location','NorthWest')

end
