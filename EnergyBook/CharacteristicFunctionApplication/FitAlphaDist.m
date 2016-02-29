function FitAlphaDist (S)
% FitAlphaDist fits an alpha-stable dist to Actual Data (S) 
% or Self-Simulation (no input) converted to log return 
% i.e., ln(St)-ln(St-1).

close all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

global PDFexp
global xout
global TimeDelta
global PDFtheo

TimeDelta=1/252; 

%assume daily prices but could add as an input to function
if (nargin == 1), %Stock Price Data Input 
    %or read directly s=load('XOMprice.dat');
    N = length(S);
    TimeLength=TimeDelta*N; 
    time = linspace(0,TimeLength,N);
    %logarithmic returns log(St/S0) 
    LogDelta=log(S(2:end))-log(S(1:end-1));
end
if (nargin == 0), %self-simulation 

    Szero=30; %Initial Price
    N = 7.5*252 %assume 7.5 years
    sigma=0.4
    mu=0.1
    beta=0.1
    alpha=1.8
    %Calc vector of random movements drawn
   
    LogDelta= TimeDelta*GenAlphaDist (alpha, beta, mu, sigma,  N);

    S=zeros(1,N); S(1)=Szero;
    for j = 2:N
        S(j)=S(j-1)*exp(LogDelta(j));
    end
end

figure %Graph Asset Movement vs. time
subplot (3,1,1); plot (2003+(1:N)/252,S);
axis tight
if (nargin == 0)
   title ('Asset Price Self Simulation ')
elseif (nargin ==1)
   title ('XOM Asset Price')
end
ylabel('Price'); xlabel('Time [Years]');

%Bin experimental Data
[PDFraw,xout] = hist(LogDelta, round(sqrt(N)));
PDFexp=PDFraw/sum(PDFraw); %Normalize Total Distribution to One

estAlpha=1.7; estBeta=0.1; estMu=0.1; estSigma=0.3; 

param=zeros(1,4);    
param(1)=estAlpha; param(2)=estBeta; param(3)=estMu; param(4)=estSigma; 

Opt=optimset('Display','iter');
[pnew,likelihood]=fminsearch('LikeCF',param,Opt); 
likelihood=-likelihood;
fprintf(1, 'Likelihood = %6.4f \n', likelihood);

calcAlpha=pnew(1)
calcBeta=pnew(2)
calcMu=pnew(3)
calcSigma=pnew(4)

subplot (3,1,2);  % Linear PDF
plot (xout,PDFexp,xout,PDFtheo);legend ('Exp.', 'Theo.')
ylabel('PDF(Log Return)'); xlabel('Log Return'); axis tight

subplot (3,1,3);  % Semi-Log PDF
semilogy (xout,PDFexp,xout,abs(PDFtheo));legend ('Exp.', 'Theo.')
ylabel('PDF(Log Return)'); xlabel('Log Return'); axis tight

end




