% ModelJumpDiffusion fits data to Drift-Diffusion with Jumps
% deltaS/S=mu*deltat + sigma*randn*deltat^0.5 + nu*deltaq
% nu=average jump size - proportional to increase in asset price
% deltaq = 1 w/ probability lambda*deltat
% deltaq = 0 w/ probability 1-lambda*deltat
function ModelJumpDiffusion(S)
% Actual Data (S) or Self-Simulation (no input)

close all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesColorOrder',[0 0 0;0 0 1;0 1 0; 1 0 0 ; 1 0 1; 0 1 1])
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

global M1 M2
global fjd fData Centers lengthLD dt
global LSflag

LSflag=1; % 1 for LS; 2 for MME
dt=1/252; % assume daily prices but could add as an input to function
SQRTdt=sqrt(dt);  
% Stock Price Data Input or read directly s=load('XOMprice.dat');
if (nargin == 1), 
    steps = length(S);
    Years=steps*dt
    TimeLength=dt*steps; time = linspace(0,TimeLength,steps);
    [C,R] = size(S);
    if (C>R)
      S=S'; % flip S=(column) to S'=(rows)
    end
    LnS=log(S);
end
if (nargin == 0), % self-simulation 
    Szero=50; 
    mu=0.11; % drift
    vol=0.25; % Volatility
    musig2=mu-0.5*vol^2
    lambda=5; % rate of jumps per year = Intensity of Poisson Process
    q1=-0.14; q2=0.15; %Max Jump Down, Up
    nuMean=(exp(q2)-exp(q1))/(q2-q1) - 1  
    % nuMean = Average Jump size measured relative 
    % to previous stock price
    % logNuP1=log(nuMean+1) %Drift of ln(jumps)
    Years=7.5;  steps = Years*252
    lambdadt=lambda*dt;
    TimeLength=dt*steps; 
    time = linspace(0,TimeLength,steps);  %years
    S(1)=Szero;     %LnS(1)=log(Szero);
    rand('state', 0); randn('state', 0);
    UniDist=rand(1,steps);
    % Hanson suggests using center to avoid end bias in distribution
    % jumpleft=(1-lambda*dt)/2;jumpright=(1-jumpleft);
    for i=2:steps       % Calculate simulated price process
            if (lambdadt>UniDist(i))
           % if ((UniDist(i)>=jumpleft)&&(UniDist(i)<=jumpright))
                Q=q1+(q2-q1)*rand; % Distributd-Size Jump
                % LnS(i)= LnS(i-1)+(musig2*dt+vol*randn*SQRTdt+Q);
                S(i)=S(i-1)*exp(musig2*dt+vol*randn*SQRTdt+Q);
            else % No Jumps only Drift-Diffusion
                % LnS(i)= LnS(i-1)+(musig2*dt+vol*randn*SQRTdt);
                S(i)=S(i-1)*exp(musig2*dt+vol*randn*SQRTdt);
            end
    end
    % S=exp(LnS);     
end

LnS=log(S);
% calculate as if logarithmic returns log(St/S0) 
% are normally distributed 
% LogDelta=LnS(2:end)-LnS(1:end-1); %
LogDelta=log(S(2:end))-log(S(1:end-1)); %
lengthLD=length(LogDelta);
M1=mean(LogDelta);   % Technically this is the first raw moment 
StanDev = std(LogDelta); 
M2=StanDev^2;                     % 2nd Central Moment
M3=mean((LogDelta-M1).^3);        % 3rd Central Moment
M4=mean((LogDelta-M1).^4);        % 4th Central Moment
Skew=M3/(M2^1.5);
Kurtosis=M4/(M2^2)-3;

xmin =min(LogDelta);        estQ1=xmin; 
xmax=max(LogDelta);         estQ2=xmax; 
sorted=sort(LogDelta);
q25=sorted(floor(0.25*lengthLD)); q75=sorted(floor(0.75*lengthLD));
estMuJump=(estQ1+estQ2)/2
estNuMean=(exp(estQ2)-exp(estQ1))/(estQ2-estQ1) - 1  

% Count returns +/- 3 standard deviations as number of jumps and
% divide by Years of data to estimate lambda
% recalculate diffusion volatility without outliers
% Modification of approach presented in L. Clewlow, C. Strickland,
% V. Kaminski, "Extending Mean-Reversion Jump Diffusion"
outliersBottom=0; outliersTop=0;
neg3sd=M1-3*StanDev;
pos3sd=M1+3*StanDev;
bottom = 1; % step through sorted array to find outliers
while (sorted(bottom)<neg3sd) 
    outliersBottom=outliersBottom+1;  bottom=bottom+1;
end
top = lengthLD; % Could also use Matlab 'find' function
while (sorted(top)>pos3sd)
    outliersTop=outliersTop+1;   top=top-1;
end

StanDev = std(LogDelta(bottom:top));
estVol=StanDev/sqrt(dt);% Estimated annualized volatility 

estLambda=(outliersTop+outliersBottom)/Years;
estMuDsig2=(M1-estMuJump*estLambda*dt)/dt;
estMuD=estMuDsig2+0.5*estVol^2
% bin data based on process developed in D. Synowiez, Computers
% and Mathematics with applications, 56, 2120 (2008)
k = round ((((xmax-xmin)*lengthLD^(1/3)) / (2.64*(q75-q25)) )+1);
db=(xmax-xmin)/k;
Edges=zeros(1,k+2); Centers=zeros(1,k+1);
for i = 1:(k+2)
    Edges(i)=xmin+(i-1.5)*db; % -0.5*db...(k+0.5)*db
end
for i = 1:(k+1)
    Centers(i)=xmin+(i-1)*db; % +0*db...k*db
end
% fData= Actual frequency in each bin
[fData] = histc(LogDelta,Edges); fData=fData(1:end-1);

% Estimation more stable with proper pre-estimation (above) of 
% estLambda, estQ1, estQ2 
param=zeros(1,3);   % lambda, q1, q2 are independent variables    
param(1)=estLambda; param(2)=estQ1; param(3)=estQ2;

[pnew,likelihood]=fminsearch('LikeEval',param); 
% likelihood=-likelihood;
% fprintf(1, 'Likelihood = %6.4f \n', likelihood);

calcLam=pnew(1); a=pnew(2); b=pnew(3);
    muJump=(a+b)/2;
    sigJ2=(b-a)^2/12;
    calcMuDsig2=(M1-muJump*calcLam*dt)/dt;
    calcSig=sqrt((M2-(sigJ2+muJump^2)*calcLam*dt)/dt);
    calcMuD=calcMuDsig2+0.5*calcSig^2;
    calcNuMean=((exp(b)-exp(a))/(b-a)) - 1 ; 
    
fJD= ((1-calcLam*dt)/(calcSig*sqrt(dt)))...
 *myNormPDF((Centers-(calcMuDsig2)*dt)/(calcSig*sqrt(dt)))...
 +(calcLam*dt/(b-a))*...
  ( myNormCDF((Centers-a-(calcMuDsig2)*dt)/(calcSig*sqrt(dt)))...
    - myNormCDF((Centers-b-(calcMuDsig2)*dt)/(calcSig*sqrt(dt))) );
fJD=lengthLD*fJD/sum(fJD);

fprintf(1, '\t\t\t\t mu \t vol \t\t lambda \t q1 \t q2 \n');
if(nargin == 0),
fprintf(1,...
    'Simulated  \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \n',...
            mu, vol, lambda, q1, q2);
end
fprintf(1,...
    'Estimated  \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \n',...
            estMuD, estVol, estLambda, estQ1, estQ2); 
fprintf(1,...
    'Calculated \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \n',...
            calcMuD, calcSig, pnew); 
        
estM1=(calcMuDsig2+calcLam*muJump)*dt;
estM2=(calcSig^2+calcLam*(sigJ2+muJump^2))*dt;
estM3=(3*sigJ2+muJump^2)*muJump*calcLam*dt;
estM4=(muJump^4+3*sigJ2^2+6*muJump^2*sigJ2)*calcLam*dt...
        +3*(calcSig^2+calcLam*(sigJ2+muJump^2))^2*dt^2;
estSkew=estM3/(estM2^1.5);
estKurtosis=estM4/(estM2^2)-3;

fprintf(1, '\n\t\t\t  Skew \t Kurtosis \n');
fprintf(1,'Data  \t \t %6.3f \t %6.2f  \n', Skew, Kurtosis);
fprintf(1,'Estimated  \t %6.3f \t %6.2f  \n\n',estSkew,estKurtosis);
        
GoodFit(fJD,fData,lengthLD, Centers); % Goodness of fit statistics

fDataANDfJD=[fData; fJD]';                         
% figure; bar(Centers,fDataANDfJD,'group');

figure; 
subplot (1,2,1); plot(Centers,fDataANDfJD);
xlabel('Log-Return'); ylabel('Frequency');axis tight;
if (nargin == 0) % self-simulation 
   title('JD Self-Simulation');
else
   title('JD of Asset Data');
end
subplot (1,2,2); semilogy(Centers,fDataANDfJD);
xlabel('Log-Return'); ylabel('Log-Frequency'); axis tight;
legend('Data','J/D Fit','location','South')
if (LSflag == 1) % self-simulation 
   title('with Least-Squares fit');
else
   title('with Multinomial Estimation');
end

ES=zeros(1,steps); ES(1)=S(1); %ElnS(1)=log(S(1));
for i=2:steps       % Calculate expectation
   % ElnS(i)=ElnS(i-1)+...
  %      ((calcMuD+calcNuMean*calcLam)*dt);
    ES(i)=ES(i-1)*exp((calcMuD+calcNuMean*calcLam)*dt);
end
% ES=exp(ElnS);
time=time+2003; % start in Year 2003
figure
plot(time,ES,':',time,S)
legend('J/D Expectected','Data','location','NorthWest')
title('Uniform-Jump / Diffusion')
xlabel('Time [Years])'); ylabel('Price, S_t ')
axis tight
end