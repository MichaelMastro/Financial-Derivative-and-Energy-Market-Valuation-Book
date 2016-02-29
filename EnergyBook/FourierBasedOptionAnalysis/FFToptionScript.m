% FFToptionScript compares FFT vs. analytical option
% value for Black-Scholes model and Merton model
% Lastly, the fractional FFT (FrFFT) is used to 
% decrease the k spacing to bunch the strike prices
% to more relevant values near the ATM option

close all; clear all; clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

T=1;        Ttext=['T = ' num2str(T)];    
r=0.05;     rtext=['r = ' num2str(r)];    
S0=100;     S0text=['S_0 = ' num2str(S0)];    

%%% Black-Scholes Parameters %%%
volBS = 0.2; 
Voltext=['\sigma_{BS} = ' num2str(volBS)];    
d=0;

CharFunc='phiBS'
Param.T=T;
Param.sigmaBS=volBS;
Param.S0=S0;
Param.r=r;
Param.mu=r-0.5*Param.sigmaBS^2; %risk-neutral BS

[cBSfft, k ] = FFToption ( CharFunc, Param);
K=S0*exp(k); % K is normalized to S0 in FFT calc
% use K returned from FFT option as input into BS 
cBSanalytical = BlackScholesCall(K,S0,T,volBS,r,d);
%cBSfftOTM = FFToptionOTM ( CharFunc, Param);

figure
plot (K, cBSanalytical, K, S0* cBSfft,  '+',...
    K, S0-K,'--')
legend ('Analytical BS', 'FFT BS')
xlim ([0.5, 1.5*S0]); ylim([0.1 100]);
bsStr(1)= {'Black-Scholes'};
bsStr(2)= {rtext};
bsStr(3)= {Voltext};
text(0.8*S0, 0.6*S0, bsStr)
genStr(1)= {S0text};
genStr(2)= {Ttext}; 
text(0.2*S0, 0.2*S0, genStr);

%%% Merton Series %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assume Merton  same underlying Black-Scholes 
% volatility augmented with a log-normal jump 
% process at a rate lambda
muJ = 0.1;       
sigmaJ = 0.1;    
lambda = 0.5;    
mmuJtext=['\mu_{Jump} = ' num2str(muJ)];    
msigmaJtext=['\sigma_{Jump} = ' num2str(sigmaJ)];    
lambdatext=['\lambda = ' num2str(lambda)];    

CharFunc='phiMerton'

Param.muJ=muJ;
Param.sigmaJ=sigmaJ;
Param.lambda=lambda;
Param.r=r;
% Merton risk-neutral drift
Param.mu=r-0.5*volBS^2 ...
            -lambda*(exp(muJ+0.5*sigmaJ^2)-1); 
        
cMERTONanalytical = MertonSeriesCall (K,S0,...
                T,volBS,r, d, muJ, sigmaJ, lambda);   
                    
[cMERTONfft, k ] = FFToption ( CharFunc, Param );
% Carr and Madan suggest alternative algorithm for 
% deep out of the money (OTM) options based on 
% Time Value of OTM options
% cBSfftOTM = FFToptionOTM ( CharFunc, Param );
K=S0*exp(k); %K is normalized to S0 in FFT calc
                    
figure
plot(K, cMERTONanalytical, K, S0*cMERTONfft,'+',...
    K,S0-K ,'--')
legend ('Merton Series Approx.', 'FFT Merton')
xlim ([0.5, 1.5*S0]); ylim([0.1 100]);
JumpStr(1)= {'Merton '};
JumpStr(2)= {mmuJtext};
JumpStr(3)= {msigmaJtext};
JumpStr(4)= {lambdatext};
text(0.8*S0, 0.6*S0, JumpStr)
genStr(1)= {S0text};
genStr(2)= {Ttext}; 
text(0.2*S0, 0.2*S0, genStr);    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option Price via Fractional Fourier Transform
% Normally FFT spacing is defined dw*dk =2pi/N
% or dw*dk/2pi = 1/N
% Thus fixed N and dw (or maximum w = N*dw ) 
% also sets dk by dk=2pi/(N*dw) 
% and maximum k by k(max) = N*dk
% -> Most k way too large or small
% FrFFT allows flexibility via the 'A Fraction'
% by dw*dk/2pi = A thus dk = A*2pi/dw
% and k spacing can be decreased
% In other words, max (and min k =- max k)
% are brought closer to the ATM and more option
% values are calculated at relevant strike prices
alpha=1.5; 
wEnd=1000;
figure
krange=0; % zero for just ATM; 
% Use larger number to look at range around ATM
% but this can be misleading as K spacing changes
nn=1:4
for n=nn;
    N=2^(n+6)
    ff=1:4;
    for f=ff
        A=1/(f*N); %set the 'A Fraction'
        [cBSFRfft k] =...
            FFToption(CharFunc,Param,alpha,wEnd,N);
        K=S0*exp(k); %K is normalized to S0 in FFT 
  cBSanalytical = BlackScholesCall(K,S0,T,volBS,r,d);
% Find ATM call when k~=0 ->K(ATM)=S0*exp(0)        
        [minKval,Ind] = min(k.^2) ;
        Error=S0*cBSFRfft(Ind-krange:Ind+krange)-...
                cBSanalytical(Ind-krange:Ind+krange);
        ErrorSum(n,f)=sum(Error.^2);
        subplot (max(ff),max(nn),n+(f-1)*(max(nn))); 
   plot (K, cBSanalytical, K, S0* cBSFRfft, '+',...
                K, S0-K,'--')
        ftext=['A=1/(' num2str(f) 'N) '];
        ntext=['N=' num2str(N)]; 
        nfText(1)={ftext}; nfText(2)={ntext};
        text(1.1*S0, 0.8*S0, nfText);    
        xlim ([0, 4*S0])
        ylim ([0 100])
    end
end
% Examine Error for change in N Spacing or A Fraction 
% Plot should show that FrFFT does does not change the 
% accuracy of the option calculation at a certain K
% Accuracy is determined primarily by N and dw as well
% as alpha and w-max
figure
semilogy (nn+6, ErrorSum(:,1),'+',nn+6, ErrorSum(:,2),'--',...
    nn+6, ErrorSum(:,3), '-.' ,nn+6, ErrorSum(:,4), ':')
legend('A=1/N', 'A=1/2N', 'A=1/3N' , 'A=1/4N')
xlabel ('2^n')
ylabel ('Square Error for ATM Call')
title ('Compare BS-FFT to BS-Analytical')
